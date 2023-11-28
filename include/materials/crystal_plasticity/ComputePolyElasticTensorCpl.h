//Purpose: Based on crystal plasticity finite element model + grainTracker custom class to manage elastic modulus
// Reference 1: poly_grain_growth_2D_eldrforce.i + ComputePolycrystalElasticityTensor + GrainTrackerElasticity
// Reference 2: GrainPropertyReadFileCP + ComputeElasticityTensorCP

// TODO: Consider the rotation of grains during loading

#pragma once

// ComputePolyElasticTensorCpl ~ Compute Polycrystal Elastic Tensor in Coupled Model

#include "ComputeElasticityTensorBase.h"
#include "GrainTrackerMatProp.h"

class EulerAngleProvider;

class ComputePolyElasticTensorCpl : public ComputeElasticityTensorBase
{
public:
  static InputParameters validParams();

  ComputePolyElasticTensorCpl(const InputParameters & parameters);

protected:
  virtual void computeQpElasticityTensor();

  /// vector of elasticity tensor material properties
  std::vector<MaterialProperty<RankFourTensor> *> _D_elastic_tensor;

  const GrainTrackerMatProp & _grain_tracker;

  const unsigned int _op_num; /// Number of order parameters
  const std::vector<const VariableValue *> _vals; /// Order parameters

  /// Crystal Rotation Matrix used to rotate the slip system direction and normal
  MaterialProperty<RankTwoTensor> & _crysrot;

  Real _length_scale;
  Real _pressure_scale;  
  const Real _JtoeV; /// Conversion factor from J to eV3
};