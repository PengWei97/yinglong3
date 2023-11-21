//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CPKalidindiBackstressUpdate.h"

registerMooseObject("baizeApp", CPKalidindiBackstressUpdate);
 
InputParameters
CPKalidindiBackstressUpdate::validParams()
{
  InputParameters params = CrystalPlasticityKalidindiUpdate::validParams();
  params.addClassDescription("Kalidindi version of homogeneous crystal plasticity considering backstress term.");
  params.addRequiredRangeCheckedParam<Real>("c_bs","c_bs>=0.0","Direct hardening coefficient for backstress");
  params.addRequiredRangeCheckedParam<Real>("d_bs","d_bs>=0.0","Dynamic recovery coefficient for backstress");

  return params;
}

CPKalidindiBackstressUpdate::CPKalidindiBackstressUpdate(
    const InputParameters & parameters)
  : CrystalPlasticityKalidindiUpdate(parameters),

    // Backstress parameters
    _c_bs(getParam<Real>("c_bs")),
    _d_bs(getParam<Real>("d_bs")),

    // Backstress variable
    _backstress(declareProperty<std::vector<Real>>(_base_name + "backstress")),
    _backstress_old(getMaterialPropertyOld<std::vector<Real>>(_base_name + "backstress")),
    _backstress_increment(declareProperty<std::vector<Real>>(_base_name + "backstress_increment")),

    // resize local caching vectors used for substepping
    _previous_substep_backstress(_number_slip_systems, 0.0),
    _backstress_before_update(_number_slip_systems, 0.0)
{
}

void
CPKalidindiBackstressUpdate::initQpStatefulProperties()
{
  // Initialize _tua, schmid tensor, g^\alpha, gamma^\alpha
  CrystalPlasticityKalidindiUpdate::initQpStatefulProperties();

  // Initialize the backstress size
  _backstress[_qp].resize(_number_slip_systems);
  for (const auto i : make_range(_number_slip_systems))
    _backstress[_qp][i] = 0.0;
}

void
CPKalidindiBackstressUpdate::setInitialConstitutiveVariableValues()
{
  CrystalPlasticityKalidindiUpdate::setInitialConstitutiveVariableValues();

  _backstress[_qp] = _backstress_old[_qp];
  _previous_substep_backstress = _backstress_old[_qp];
}

void
CPKalidindiBackstressUpdate::setSubstepConstitutiveVariableValues()
{
  CrystalPlasticityKalidindiUpdate::setSubstepConstitutiveVariableValues();

  _backstress[_qp] = _previous_substep_backstress;
}

bool
CPKalidindiBackstressUpdate::calculateSlipRate()
{
  // Difference between RSS and backstress temporary variable for each slip system
  Real effective_stress;

  for (const auto i : make_range(_number_slip_systems))
  {
    effective_stress = _tau[_qp][i] - _backstress[_qp][i];

    _slip_increment[_qp][i] =
        _ao * std::pow(std::abs(effective_stress / _slip_resistance[_qp][i]), 1.0 / _xm);
    if (effective_stress < 0.0)
      _slip_increment[_qp][i] *= -1.0;

    if (std::abs(_slip_increment[_qp][i]) * _substep_dt > _slip_incr_tol)
    {
      if (_print_convergence_message)
        mooseWarning("Maximum allowable slip increment exceeded ",
                     std::abs(_slip_increment[_qp][i]) * _substep_dt);

      return false;
    }
  }
  return true;
}

void
CPKalidindiBackstressUpdate::calculateConstitutiveSlipDerivative(
    std::vector<Real> & dslip_dtau)
{
  Real effective_stress;

  for (const auto i : make_range(_number_slip_systems))
  {
    effective_stress = _tau[_qp][i] - _backstress[_qp][i];

    if (MooseUtils::absoluteFuzzyEqual(effective_stress, 0.0))
      dslip_dtau[i] = 0.0;
    else
      dslip_dtau[i] = _ao / _xm *
                      std::pow(std::abs(effective_stress / _slip_resistance[_qp][i]), 1.0 / _xm - 1.0) /
                      _slip_resistance[_qp][i];
  }
}

bool
CPKalidindiBackstressUpdate::areConstitutiveStateVariablesConverged()
{
  if (isConstitutiveStateVariableConverged(_slip_resistance[_qp],
                                           _slip_resistance_before_update,
                                           _previous_substep_slip_resistance,
                                           _resistance_tol) &&
      isConstitutiveStateVariableConverged(_backstress[_qp],
                                           _backstress_before_update,
                                           _previous_substep_backstress,
                                           _rel_state_var_tol))
    return true;
  return false;                      
}

void
CPKalidindiBackstressUpdate::updateSubstepConstitutiveVariableValues()
{
  _previous_substep_slip_resistance = _slip_resistance[_qp];
  _previous_substep_backstress = _backstress[_qp];
}

void
CPKalidindiBackstressUpdate::cacheStateVariablesBeforeUpdate()
{
  _slip_resistance_before_update = _slip_resistance[_qp];
  _backstress_before_update = _backstress[_qp];
}

void
CPKalidindiBackstressUpdate::calculateStateVariableEvolutionRateComponent()
{
  CrystalPlasticityKalidindiUpdate::calculateStateVariableEvolutionRateComponent();
  
  ArmstrongFrederickBackstressUpdate();
}

bool
CPKalidindiBackstressUpdate::updateStateVariables()
{
  // Now perform the check to see if the slip system should be updated
  if (CrystalPlasticityKalidindiUpdate::updateStateVariables())
  {
    // now update backstress variable, which can be positive or negative
    for (const auto i : make_range(_number_slip_systems))
    { 
      _backstress_increment[_qp][i] *= _substep_dt;
      _backstress[_qp][i] = _previous_substep_backstress[i] + _backstress_increment[_qp][i];
    }
    return true;
  }

  return false;
}

void
CPKalidindiBackstressUpdate::ArmstrongFrederickBackstressUpdate()
{
  _backstress_increment[_qp].resize(_number_slip_systems, 0.0);
  for (const auto i : make_range(_number_slip_systems)) 
  {
    _backstress_increment[_qp][i] = _c_bs * _slip_increment[_qp][i];
    _backstress_increment[_qp][i] -= _d_bs * _backstress[_qp][i] * std::abs(_slip_increment[_qp][i]);  
  }
}