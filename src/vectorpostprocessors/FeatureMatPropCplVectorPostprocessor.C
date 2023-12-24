#include "FeatureMatPropCplVectorPostprocessor.h"

registerMooseObject("baizeApp", FeatureMatPropCplVectorPostprocessor);

InputParameters
FeatureMatPropCplVectorPostprocessor::validParams()
{
  InputParameters params = FeatureDataVectorPostprocessor::validParams();

  params.addParam<std::string>("base_name", "Name to append to reporters."); 
  params.addParam<unsigned int>("number_slip_systems", 12,
                                "The total number of possible active slip systems for the crystalline material");  

  return params;                                
}

FeatureMatPropCplVectorPostprocessor::FeatureMatPropCplVectorPostprocessor(
    const InputParameters & parameters)
  : FeatureDataVectorPostprocessor(parameters),
    _number_slip_systems(getParam<unsigned int>("number_slip_systems")),
    _slip_resistances(declareRestartableData<std::vector<std::vector<Real>>>("slip_resistances")),    
    _backstresses(declareRestartableData<std::vector<std::vector<Real>>>("backstresses")),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _slip_resistance_copy(getMaterialProperty<std::vector<Real>>(_base_name + "slip_resistance_copy")),    
    _backstress_copy(getMaterialProperty<std::vector<Real>>(_base_name + "backstress_copy"))
{
}

void
FeatureMatPropCplVectorPostprocessor::execute()
{
  FeatureDataVectorPostprocessor::execute();

  const auto num_features = _feature_counter.getTotalFeatureCount();

  // Reset the state variables vector
  _slip_resistances.resize(num_features);
  _backstresses.resize(num_features);
  for (auto grain_index : make_range(num_features))
  {
    _slip_resistances[grain_index].assign(_number_slip_systems, 0.0);
    _backstresses[grain_index].assign(_number_slip_systems, 0.0);
  }

  for (const auto & elem : _mesh.getMesh().active_local_element_ptr_range())
  {
    _fe_problem.setCurrentSubdomainID(elem, 0);
    _fe_problem.prepare(elem, 0);
    _fe_problem.reinitElem(elem, 0);

    /**
     * Here we retrieve the var to features vector on the current element.
     * We'll use that information to figure out which variables are non-zero
     * (from a threshold perspective) then we can sum those values into
     * appropriate grain index locations.
     */
    const auto & var_to_features = _feature_counter.getVarToFeatureVector(elem->id());

    accumulateStateVariables(elem, var_to_features, num_features);
  }
}

void
FeatureMatPropCplVectorPostprocessor::accumulateStateVariables(
    const Elem * elem,
    const std::vector<unsigned int> & var_to_features,
    std::size_t libmesh_dbg_var(num_features))
{
  for (MooseIndex(var_to_features) var_index = 0; var_index < var_to_features.size(); ++var_index)
  {
    // Only sample "active" variables
    if (var_to_features[var_index] != FeatureFloodCount::invalid_id)
    {
      auto feature_id = var_to_features[var_index];
      mooseAssert(feature_id < num_features, "Feature ID out of range");

      std::vector<Real> slip_resistance_integral_value = computeStateVarivalesIntegral(var_index, _slip_resistance_copy);

      std::vector<Real> backstress_integral_value = computeStateVarivalesIntegral(var_index, _backstress_copy);

      for (auto sr_index : make_range(_number_slip_systems))
      {
        _slip_resistances[feature_id][sr_index] += slip_resistance_integral_value[sr_index];
        _backstresses[feature_id][sr_index] += backstress_integral_value[sr_index];
      } 
    }  
  }  
}

std::vector<Real> 
FeatureMatPropCplVectorPostprocessor::computeStateVarivalesIntegral(std::size_t var_index, const MaterialProperty<std::vector<Real>> & state_variable_copy) const
{
  std::vector<Real> sum(_number_slip_systems, 0.0);

  for (auto sr_index : make_range(_number_slip_systems))
    for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
      sum[sr_index] += _JxW[qp] * _coord[qp] * (*_coupled_sln[var_index])[qp] * state_variable_copy[qp][sr_index];

  return sum;
}

void
FeatureMatPropCplVectorPostprocessor::finalize()
{
  // Do the parallel sum
  _communicator.sum(_feature_volumes);

  auto num_features = _feature_volumes.size();
  sum_state_variables(_slip_resistances, num_features);
  sum_state_variables(_backstresses, num_features);
}

void 
FeatureMatPropCplVectorPostprocessor::sum_state_variables(std::vector<std::vector<Real>> & stat_variables, const unsigned int & num_features)
{
  for (auto grain_index : make_range(num_features))
    _communicator.sum(stat_variables[grain_index]);

  for (auto grain_index : make_range(num_features))
  {
    if (_feature_volumes[grain_index] > 0.0)
      for (auto sr_index : make_range(_number_slip_systems))
        stat_variables[grain_index][sr_index] /= _feature_volumes[grain_index];
    else
      stat_variables[grain_index].assign(_number_slip_systems, 0.0);
  }
}

std::vector<Real> 
FeatureMatPropCplVectorPostprocessor::getStateVariable(unsigned int feature_id, const state_variable & state_variable_name)
{
  mooseAssert(feature_id < _slip_resistances.size(), "feature_id is out of range");

  switch (state_variable_name) 
  {
    case state_variable::slip_resistance:
        return _slip_resistances[feature_id];
    case state_variable::backstress:
        return _backstresses[feature_id];
  }

  return _backstresses[feature_id];
}


