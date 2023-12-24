#include "CPKalidindiBSCplUpdate.h"

registerMooseObject("baizeApp", CPKalidindiBSCplUpdate);
 
InputParameters
CPKalidindiBSCplUpdate::validParams()
{
  InputParameters params = CPKalidindiBackstressUpdate::validParams();
  params.addRequiredParam<UserObjectName>(
      "grain_tracker", "Name of GrainTracker user object that provides RankFourTensors");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  params.addRequiredParam<VectorPostprocessorName>(
      "vpp_name", "The name of the VectorPostprocessor that you want to use");

  return params;
}

CPKalidindiBSCplUpdate::CPKalidindiBSCplUpdate(
    const InputParameters & parameters)
  : CPKalidindiBackstressUpdate(parameters),
    _slip_resistance_copy(declareProperty<std::vector<Real>>(_base_name + "slip_resistance_copy")),
    _backstress_copy(declareProperty<std::vector<Real>>(_base_name + "backstress_copy")),
    _grain_tracker(getUserObject<GrainTrackerMatProp>("grain_tracker")),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _number_active_grains(declareProperty<unsigned int>("number_active_grains")),
    _number_active_grains_old(getMaterialPropertyOld<unsigned int>("number_active_grains")),
    _fe_problem(*parameters.get<FEProblem *>("_fe_problem")),
    _t_step(_fe_problem.timeStep()),
    _vpp_name(getParam<VectorPostprocessorName>("vpp_name"))
{
}

void
CPKalidindiBSCplUpdate::initQpStatefulProperties()
{
  CPKalidindiBackstressUpdate::initQpStatefulProperties();

  _slip_resistance_copy[_qp] = _slip_resistance[_qp];
  _backstress_copy[_qp] = _backstress[_qp];

  _number_active_grains[_qp] = 0;
}

void
CPKalidindiBSCplUpdate::setInitialConstitutiveVariableValues()
{
  CPKalidindiBackstressUpdate::setInitialConstitutiveVariableValues();

  if (_t_step > 0)
    convertStateVariablesFromPFtoCP();
}

void
CPKalidindiBSCplUpdate::convertStateVariablesFromPFtoCP()
{
  // Get list of active order parameters from grain tracker
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

  std::vector<unsigned int> active_op_indexs;
  unsigned int max_op_index = 0;
  Real max_op_value = (*_vals[max_op_index])[_qp];
  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue; 

    active_op_indexs.push_back(op_index);
    if (max_op_value < (*_vals[op_index])[_qp])
    {
      max_op_index = op_index;
      max_op_value = (*_vals[max_op_index])[_qp];
    }
  }

  _number_active_grains[_qp] = active_op_indexs.size();

  if (_number_active_grains[_qp] == 1) // inside grain
  {
    if (_number_active_grains[_qp] == _number_active_grains_old[_qp])
    {
      _slip_resistance[_qp] = _slip_resistance_old[_qp];
      _backstress[_qp] = _backstress_old[_qp];
      }
    else
    {
      _slip_resistance[_qp].resize(_number_slip_systems, _gss_initial);
      _backstress[_qp].resize(_number_slip_systems, 0.0);
    }
  }
  else if (_number_active_grains[_qp] > 1) // at GB or junction
  {
    std::vector<std::vector<Real>> grain_to_resistance_slips(active_op_indexs.size());
    std::vector<std::vector<Real>> grain_to_backstresses(active_op_indexs.size());

    // TODO： want need set it to be _vpp_object_ptr
    FeatureMatPropCplVectorPostprocessor * vpp_object_ptr = dynamic_cast<FeatureMatPropCplVectorPostprocessor *>(const_cast<VectorPostprocessor *>(&_fe_problem.getVectorPostprocessorObjectByName(_vpp_name)));

    if (!vpp_object_ptr)
      mooseError("Pointer cast failed! Object is not of expected type");    
    
    if (_number_active_grains[_qp] > _number_active_grains_old[_qp]) // have some active grains
    {
      std::vector<Real> vector_active_op_val;
      for (auto op_index : active_op_indexs)
        vector_active_op_val.push_back((*_vals[op_index])[_qp]);

      std::sort(vector_active_op_val.rbegin(), vector_active_op_val.rend()); // descending sort

      unsigned int new_active_number = _number_active_grains[_qp] - _number_active_grains_old[_qp];
      Real critical_op_val = vector_active_op_val[new_active_number];

      for (unsigned int i = 0; i < active_op_indexs.size(); i++)
      {
        auto & op_index = active_op_indexs[i];
        Real op_val = (*_vals[op_index])[_qp];
        auto grain_id = op_to_grains[op_index];

        if (op_val >= critical_op_val && op_val == max_op_index)
        {
          grain_to_resistance_slips[i] = _slip_resistance_old[_qp];
          grain_to_backstresses[i] = _backstress_old[_qp];
        }
        else if (op_val >= critical_op_val)
        {
          grain_to_resistance_slips[i] = vpp_object_ptr->getStateVariable(grain_id, state_variable::slip_resistance); // 
          grain_to_backstresses[i] = vpp_object_ptr->getStateVariable(grain_id, state_variable::backstress); // 
        }
        else
        {
          grain_to_resistance_slips[i].resize(_number_slip_systems, _gss_initial);
          grain_to_backstresses[i].resize(_number_slip_systems, 0.0);
        }
      }
    }
    else // No new activated grains or reduce the number of activated grains
    {
      for (unsigned int i = 0; i < active_op_indexs.size(); i++)
      {
        auto & op_index = active_op_indexs[i];
        if (op_index != max_op_index) // second: _slip_resistance_ave
        {
          auto grain_id = op_to_grains[op_index];
          grain_to_resistance_slips[i] = vpp_object_ptr->getStateVariable(grain_id, state_variable::slip_resistance); // 
          grain_to_backstresses[i] = vpp_object_ptr->getStateVariable(grain_id, state_variable::backstress); // 
        }
        else // primary
        {
          grain_to_resistance_slips[i] = _slip_resistance_old[_qp];
          grain_to_backstresses[i] = _backstress_old[_qp];
        }
      }
    }

    // Interpolate _slip_resistance[_qp]
    Real sum_h = 0.0;
    _slip_resistance[_qp].clear();
    std::vector<Real> sum_slip_resistance(_number_slip_systems, 0.0);
    std::vector<Real> sum_backstress(_number_slip_systems, 0.0);

    for (unsigned int i = 0; i < active_op_indexs.size(); i++)
    {
      auto & op_index = active_op_indexs[i];
      
      Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5))) / 2.0;
      sum_h += h;

      for (unsigned int sr_index = 0; sr_index < _number_slip_systems; sr_index++)
      {
        sum_slip_resistance[sr_index] += grain_to_resistance_slips[i][sr_index] * h;
        sum_backstress[sr_index] += grain_to_backstresses[i][sr_index] * h;
      }
    }

    const Real tol = 1.0e-10;
    sum_h = std::max(sum_h, tol);

    for (unsigned int sr_index = 0; sr_index < _number_slip_systems; sr_index++)
    {
      _slip_resistance[_qp][sr_index] = sum_slip_resistance[sr_index] / sum_h;
      _backstress[_qp][sr_index] = sum_backstress[sr_index] / sum_h;
    }
  }

  // _previous_substep_slip_resistance = _slip_resistance[_qp];
  // _previous_substep_backstress = _backstress[_qp];
}