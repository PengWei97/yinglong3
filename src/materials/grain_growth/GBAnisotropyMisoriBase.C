//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GBAnisotropyMisoriBase.h"
#include "MooseMesh.h"

registerMooseObject("baizeApp", GBAnisotropyMisoriBase);

InputParameters
GBAnisotropyMisoriBase::validParams()
{
  InputParameters params = Material::validParams();
  params.addCoupledVar("T", 300.0, "Temperature in Kelvin");
  params.addParam<Real>("length_scale", 1.0e-9, "Length scale in m, where default is nm");
  params.addParam<Real>("time_scale", 1.0e-9, "Time scale in s, where default is ns");
  params.addParam<Real>("GBsigma_HAGB", 0.9, "GB energy of a high angle GB");
  params.addParam<Real>("GBmob_HAGB", 2.5e-6, "GB mobility of a high angle GB");
  params.addParam<Real>("Q_HAGB", 0.23, "initial value of Q, where default is eV");
  params.addRequiredParam<Real>("wGB", "Diffuse GB width in nm");
  params.addParam<Real>("scale_factor_matrix", 1.0e-2,"scale factor for matrix sigma or mob");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  return params;
}

GBAnisotropyMisoriBase::GBAnisotropyMisoriBase(const InputParameters & parameters)
  : Material(parameters),
    _mesh_dimension(_mesh.dimension()),
    _length_scale(getParam<Real>("length_scale")),
    _time_scale(getParam<Real>("time_scale")),
    _GBsigma_HAGB(getParam<Real>("GBsigma_HAGB")),
    _GBmob_HAGB(getParam<Real>("GBmob_HAGB")),
    _Q_HAGB(getParam<Real>("Q_HAGB")),
    _wGB(getParam<Real>("wGB")),
    _scale_factor_matrix(getParam<Real>("scale_factor_matrix")),
    _T(coupledValue("T")),
    _kappa(declareProperty<Real>("kappa_op")),
    _gamma(declareProperty<Real>("gamma_asymm")),
    _L(declareProperty<Real>("L")),
    _mu(declareProperty<Real>("mu")),
    _kb(8.617343e-5),      // Boltzmann constant in eV/K
    _JtoeV(6.24150974e18), // Joule to eV conversion
    _mu_qp(0.0),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _grad_vals(coupledGradients("v"))
{
  // reshape vectors
  _sigma.resize(_op_num);
  _mob.resize(_op_num);
  _Q.resize(_op_num);
  _kappa_gamma.resize(_op_num);
  _a_g2.resize(_op_num);

  for (unsigned int op = 0; op < _op_num; ++op)
  {
    _sigma[op].resize(_op_num);
    _mob[op].resize(_op_num);
    _Q[op].resize(_op_num);
    _kappa_gamma[op].resize(_op_num);
    _a_g2[op].resize(_op_num);
  }
}

void
GBAnisotropyMisoriBase::computeQpProperties()
{
  // Initialize sigma_ij, mob_ij, Q_ij
  std::fill(_sigma.begin(), _sigma.end(), std::vector<Real>(_op_num, _GBsigma_HAGB));
  std::fill(_mob.begin(), _mob.end(), std::vector<Real>(_op_num, _GBmob_HAGB));
  std::fill(_Q.begin(), _Q.end(), std::vector<Real>(_op_num, _Q_HAGB));

  computeGBProperties();

  // convert unit
  for (unsigned int m = 0; m < _op_num - 1; ++m)
    for (unsigned int n = m + 1; n < _op_num; ++n)
    {
      _sigma[m][n] *= _JtoeV * (_length_scale * _length_scale); // eV/nm^2

      _mob[m][n] *= _time_scale / (_JtoeV * (_length_scale * _length_scale * _length_scale *
                                             _length_scale)); // Convert to nm^4/(eV*ns);
    } 

    computerPFParameters();

    interpolatePFParams();
}

void
GBAnisotropyMisoriBase::computeGBProperties()
{
}

void
GBAnisotropyMisoriBase::computerPFParameters()
{
  Real sigma_init;
  Real g2 = 0.0;
  Real f_interf = 0.0;
  Real a_0 = 0.75;
  Real a_star = 0.0;
  Real kappa_star = 0.0;
  Real gamma_star = 0.0;
  Real y = 0.0; // 1/gamma
  Real yyy = 0.0;

  Real sigma_big = 0.0;
  Real sigma_small = 0.0;

  for (unsigned int m = 0; m < _op_num - 1; ++m)
    for (unsigned int n = m + 1; n < _op_num; ++n)
    {
      if (m == 0 && n == 1)
      {
        sigma_big = _sigma[m][n];
        sigma_small = sigma_big;
      }

      else if (_sigma[m][n] > sigma_big)
        sigma_big = _sigma[m][n];

      else if (_sigma[m][n] < sigma_small)
        sigma_small = _sigma[m][n];
    }

  sigma_init = (sigma_big + sigma_small) / 2.0;
  _mu_qp = 6.0 * sigma_init / _wGB;

  for (unsigned int m = 0; m < _op_num - 1; ++m)
    for (unsigned int n = m + 1; n < _op_num; ++n) // m<n
    {
      a_star = a_0;
      a_0 = 0.0;

      while (std::abs(a_0 - a_star) > 1.0e-9)
      {
        a_0 = a_star;
        kappa_star = a_0 * _wGB * _sigma[m][n];
        g2 = _sigma[m][n] * _sigma[m][n] / (kappa_star * _mu_qp);
        y = -5.288 * g2 * g2 * g2 * g2 - 0.09364 * g2 * g2 * g2 + 9.965 * g2 * g2 - 8.183 * g2 +
            2.007;
        gamma_star = 1 / y;
        yyy = y * y * y;
        f_interf = 0.05676 * yyy * yyy - 0.2924 * yyy * y * y + 0.6367 * yyy * y - 0.7749 * yyy +
                   0.6107 * y * y - 0.4324 * y + 0.2792;
        a_star = std::sqrt(f_interf / g2);
      }

      _kappa_gamma[m][n] = kappa_star; // upper triangle stores the discrete set of kappa values
      _kappa_gamma[n][m] = gamma_star; // lower triangle stores the discrete set of gamma values

      _a_g2[m][n] = a_star; // upper triangle stores "a" data.
      _a_g2[n][m] = g2;     // lower triangle stores "g2" data.
    }
}

void
GBAnisotropyMisoriBase::interpolatePFParams()
{
  Real sum_kappa = 0.0;
  Real sum_gamma = 0.0;
  Real sum_L = 0.0;
  Real Val = 0.0;
  Real sum_val = 0.0;
  Real f_sigma = 1.0;
  Real f_mob = 1.0;
  Real gamma_value = 0.0;

  for (unsigned int m = 0; m < _op_num - 1; ++m)
  {
    for (unsigned int n = m + 1; n < _op_num; ++n) // m<n
    {
      gamma_value = _kappa_gamma[n][m];

      Val = (100000.0 * ((*_vals[m])[_qp]) * ((*_vals[m])[_qp]) + 0.01) *
            (100000.0 * ((*_vals[n])[_qp]) * ((*_vals[n])[_qp]) + 0.01);

      sum_val += Val;
      sum_kappa += _kappa_gamma[m][n] * f_sigma * Val;
      sum_gamma += gamma_value * Val;
      
      // Following comes from substituting Eq. (36c) from the paper into (36b)
      sum_L += Val * _mob[m][n] * std::exp(-_Q[m][n] / (_kb * _T[_qp])) * f_mob * _mu_qp * _a_g2[n][m] / _sigma[m][n];
    }
  }

  _kappa[_qp] = sum_kappa / sum_val;
  _gamma[_qp] = sum_gamma / sum_val;
  _L[_qp] = sum_L / sum_val;
  _mu[_qp] = _mu_qp;
}

void
GBAnisotropyMisoriBase::computerPFGBIsotropy()
{
  const Real length_scale4 = _length_scale * _length_scale * _length_scale * _length_scale;
  Real sigma_iso = _GBsigma_HAGB * _JtoeV * _length_scale * _length_scale;
  Real mob_iso = _GBmob_HAGB * _time_scale / (_JtoeV * length_scale4) * std::exp(-_Q_HAGB / (_kb * _T[_qp]));

  _kappa[_qp] = 3.0 / 4.0 * sigma_iso * _wGB;
  _gamma[_qp] = 1.5;
  _L[_qp] = 4.0 / 3.0 * mob_iso / _wGB;
  _mu[_qp] = 3.0 / 4.0 * 1.0 / 0.125 * sigma_iso / _wGB;
}