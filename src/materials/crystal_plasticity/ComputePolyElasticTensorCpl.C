#include "ComputePolyElasticTensorCpl.h"
#include "RotationTensor.h"
#include "EulerAngles.h"

registerMooseObject("baizeApp", ComputePolyElasticTensorCpl);

InputParameters
ComputePolyElasticTensorCpl::validParams()
{
  InputParameters params = ComputeElasticityTensorBase::validParams();
  params.addClassDescription(
      "Compute an evolving elasticity tensor coupled to a grain growth phase-field and crystal plasticity model.");
  params.addRequiredParam<UserObjectName>(
      "grain_tracker", "Name of GrainTracker user object that provides RankFourTensors");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  params.addParam<Real>("length_scale", 1.0e-9, "Length scale of the problem, in meters");
  params.addParam<Real>("pressure_scale", 1.0e6, "Pressure scale of the problem, in pa");
  return params;
}

ComputePolyElasticTensorCpl::ComputePolyElasticTensorCpl(
    const InputParameters & parameters)
  : ComputeElasticityTensorBase(parameters),
    _D_elastic_tensor(_op_num),
    _grain_tracker(getUserObject<GrainTrackerMatProp>("grain_tracker")),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _crysrot(declareProperty<RankTwoTensor>(_base_name + "crysrot")),
    _length_scale(getParam<Real>("length_scale")),
    _pressure_scale(getParam<Real>("pressure_scale")),
    _JtoeV(6.24150974e18)
{
  _D_elastic_tensor.resize(_op_num);

  // Loop over variables (ops)
  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
    _D_elastic_tensor[op_index] = &declarePropertyDerivative<RankFourTensor>(
        _elasticity_tensor_name, coupledName("v", op_index)); // declare elasticity tensor derivative properties
        
  // the base class guarantees constant in time, but in this derived class the tensor will rotate over time once plastic deformation sets in
  // TODO - Need to consider the issue of grain rotation
  revokeGuarantee(_elasticity_tensor_name, Guarantee::CONSTANT_IN_TIME);  
}

void
ComputePolyElasticTensorCpl::computeQpElasticityTensor()
{
  // Get list of active order parameters from grain tracker
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

  // Calculate elasticity tensor
  _elasticity_tensor[_qp].zero();

  Real sum_h = 0.0;

  // 找到最大的 _id，用于设定[_qp]处的欧拉角
  unsigned int max_id = 0;
  Real max_var = (*_vals[max_id])[_qp];

  // 对 _elasticity_tensor 执行插值处理
  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    // Interpolation factor for elasticity tensors
    Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5))) / 2.0;

    // Sum all rotated elasticity tensors
    _elasticity_tensor[_qp] += _grain_tracker.getElasticTensor(grain_id) * h;
    sum_h += h;

    if ((*_vals[op_index])[_qp] > max_var)
    {
      max_id = op_index;
      max_var = (*_vals[max_id])[_qp];
    }
  }

  const Real tol = 1.0e-10;
  sum_h = std::max(sum_h, tol);
  _elasticity_tensor[_qp] /= sum_h;

  // Calculate elasticity tensor derivative: Cderiv = dhdopi/sum_h * (Cop - _Cijkl)
  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
    (*_D_elastic_tensor[op_index])[_qp].zero();

  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    Real dhdopi = libMesh::pi * std::cos(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5)) / 2.0;
    RankFourTensor & C_deriv = (*_D_elastic_tensor[op_index])[_qp];

    C_deriv = (_grain_tracker.getElasticTensor(grain_id) - _elasticity_tensor[_qp]) * dhdopi / sum_h;

    // Convert from XPa to eV/(xm)^3, where X is pressure scale and x is length scale;
    C_deriv *= _JtoeV * (_length_scale * _length_scale * _length_scale) * _pressure_scale;
  }

  // TODO - 最大的晶粒索引对应的欧拉角 ~ 是否合理？
  auto grain_id = op_to_grains[max_id];
  auto & angles = _grain_tracker.getEulerAngles(max_id);

  RotationTensor R(angles);
  _crysrot[_qp] = R.transpose();
}