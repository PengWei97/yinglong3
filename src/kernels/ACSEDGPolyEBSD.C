//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ACSEDGPolyEBSD.h"

registerMooseObject("baizeApp", ACSEDGPolyEBSD);

InputParameters
ACSEDGPolyEBSD::validParams()
{
  InputParameters params = ACBulk<Real>::validParams();
  params.addClassDescription("Stored Energy contribution to grain growth");
  params.addRequiredCoupledVar("v", "Array of coupled variable names");
  params.addParam<Real>("length_scale", 1.0e-9, "Length scale in m, where default is nm");
  params.addParam<bool>("concurrent_recovery", false, "The concurrent recovery would be considered if true");

  params.addParam<Real>("rho_end_l2", 2.10e12, "Dislocation density of grains after long-term concurrent recovery in level 2");
  params.addParam<Real>("rho_end_l3", 1.8e13, "Dislocation density of grains after long-term concurrent recovery in level 3");
  params.addParam<Real>("rho_critical", 3.9e13, "Critual dislocation density of grains");
  params.addParam<Real>("a_rho_l2", 4.6e-4, "Evolution coefficient of dislocation during medium time recovery at level 2"); 
  params.addParam<Real>("a_rho_l3", 6.0e-4, "Evolution coefficient of dislocation during medium time recovery at level 3");

  params.addRequiredParam<UserObjectName>("grain_tracker",
                                          "The GrainTracker UserObject to get values from.");
  params.addRequiredParam<UserObjectName>("GNDs_provider", "GNDs provider for EBSD reader");                                    
  params.addRequiredParam<unsigned int>("op_index", "The index for the current order parameter");
  return params;
}

ACSEDGPolyEBSD::ACSEDGPolyEBSD(const InputParameters & parameters)
  : ACBulk<Real>(parameters),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _vals_var(coupledIndices("v")),
    _length_scale(getParam<Real>("length_scale")),
    _concurrent_recovery(getParam<bool>("concurrent_recovery")),  

    _rho_end_l2(getParam<Real>("rho_end_l2")),
    _rho_end_l3(getParam<Real>("rho_end_l3")),
    _rho_critical(getParam<Real>("rho_critical")),
    _a_rho_l2(getParam<Real>("a_rho_l2")),
    _a_rho_l3(getParam<Real>("a_rho_l3")), 
       
    _beta(getMaterialProperty<Real>("beta")),
    _rho_eff(getMaterialProperty<Real>("rho_eff")),
    _grain_tracker(getUserObject<GrainTrackerInterface>("grain_tracker")),
    _GNDs_provider(getUserObject<EBSDReader>("GNDs_provider")),
    _op_index(getParam<unsigned int>("op_index"))
{
}

Real
ACSEDGPolyEBSD::computeDFDOP(PFFunctionType type)
{
  Real SumEtaj = 0.0;
  for (unsigned int i = 0; i < _op_num; ++i)
    SumEtaj += (*_vals[i])[_qp] * (*_vals[i])[_qp];

  // Add the current OP to the sum
  Real SumEtai2 = SumEtaj + _u[_qp] * _u[_qp];

  const auto & op_to_grain = _grain_tracker.getVarToFeatureVector(_current_elem->id());
  const auto grn_index = op_to_grain[_op_index]; // Grain ID

  unsigned int num_vaild_grain = 0;
  for (MooseIndex(op_to_grain) op_index = 0; op_index < op_to_grain.size(); ++op_index)
  {

    if (op_to_grain[op_index] == FeatureFloodCount::invalid_id)
      continue;

    num_vaild_grain += 1;
  }
  
  Real rho_i = _rho_eff[_qp];
  
  if (num_vaild_grain >= 2)
    rho_i = getGNDsFromEBSD(grn_index);

  // Calculate the contributions of the deformation energy to the residual and Jacobian
  Real drho_eff_detai = 2.0 * _u[_qp] * (rho_i - _rho_eff[_qp]) / SumEtai2;

  // Calculate the Stored Energy contribution to either the residual or Jacobian of the grain growth
  // free energy
  switch (type)
  {
    case Residual:
    {
      if (_fe_problem.time() > 0.01)
        return _beta[_qp] * drho_eff_detai;
      else
        return 0.0;
    }

    case Jacobian:
      return 0.0; // _beta[_qp] * _phi[_j][_qp] *
            //  (2.0 * SumEtai2 * ((rho_i - _rho_eff[_qp]) - _u[_qp] * drho_eff_detai) -
            //   4.0 * _u[_qp] * _u[_qp] * (rho_i - _rho_eff[_qp])) /
            //  (SumEtai2 * SumEtai2);
  }
  mooseError("Invalid type passed in");
}

Real
ACSEDGPolyEBSD::getGNDsFromEBSD(const unsigned int & grain_id)
{
  Real rho_init = 2.0e15 * (_length_scale * _length_scale); // TODO - need to be more physically based
  Real rho_i = rho_init;

  auto & time_current = _fe_problem.time(); // current simulation time s

  if (grain_id < _GNDs_provider.getGrainNum())
    rho_init = _GNDs_provider.getAvgData(grain_id)._custom[0] * (_length_scale * _length_scale); 

  const Real rho_end_l2 = _rho_end_l2 * (_length_scale * _length_scale);
  const Real rho_end_l3 = _rho_end_l3 * (_length_scale * _length_scale);
  const Real rho_critical = _rho_critical * (_length_scale * _length_scale);
  
  if (_concurrent_recovery)
  {
    if (rho_init >= rho_end_l3) // yellow grain
      rho_i = (rho_init - rho_end_l3) * std::exp(-_a_rho_l3 * time_current) + rho_end_l3; // func1 - level 3
    else if ((rho_init > rho_end_l2) && (rho_init < rho_critical)) // blue grain
      rho_i = (rho_init - rho_end_l2) * std::exp(-_a_rho_l2 * time_current) + rho_end_l2; // func2 - level 2
    else
      rho_i = rho_init;
  }
  else
    rho_i = rho_init;

  return rho_i; 
}
