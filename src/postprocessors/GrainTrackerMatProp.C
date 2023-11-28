//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GrainTrackerMatProp.h"
#include "EulerAngleProvider.h"
#include "RotationTensor.h"

registerMooseObject("baizeApp", GrainTrackerMatProp);

InputParameters
GrainTrackerMatProp::validParams()
{
  InputParameters params = GrainTracker::validParams();
  params.addClassDescription(
      "Used to manage the material parameters corresponding to each grain id, such as elastic tensor, euler angles");

  params.addRequiredParam<std::vector<Real>>("C_ijkl", "Unrotated stiffness tensor");
  params.addParam<MooseEnum>(
      "fill_method", RankFourTensor::fillMethodEnum() = "symmetric9", "The fill method");
  params.addRequiredParam<UserObjectName>("euler_angle_provider",
                                          "Name of Euler angle provider user object");
  params.addParam<bool>("random_rotations",
                        true,
                        "Generate random rotations when the Euler Angle "
                        "provider runs out of data (otherwise error "
                        "out)");
  return params;
}

GrainTrackerMatProp::GrainTrackerMatProp(const InputParameters & parameters)
  : GrainTracker(parameters),
    _euler_angles_vec(declareRestartableData<std::vector<EulerAngles>>("euler_angles_vec")),
    _C_ijkl_vec(declareRestartableData<std::vector<RankFourTensor>>("C_ijkl_vec")),
    _random_rotations(getParam<bool>("random_rotations")),
    _C_ijkl(getParam<std::vector<Real>>("C_ijkl"),
            getParam<MooseEnum>("fill_method").getEnum<RankFourTensor::FillMethod>()),
    _euler(getUserObject<EulerAngleProvider>("euler_angle_provider"))
{
}

const EulerAngles & 
GrainTrackerMatProp::getEulerAngles(unsigned int grain_id) const
{
  mooseAssert(grain_id < _euler_angles_vec.size(), "Requested data for invalid grain index.");
  return _euler_angles_vec[grain_id];
}

const RankFourTensor & 
GrainTrackerMatProp::getElasticTensor(unsigned int grain_id) const
{
  mooseAssert(grain_id < _C_ijkl_vec.size(), "Requested data for invalid grain index.");
  return _C_ijkl_vec[grain_id];
}

void
GrainTrackerMatProp::newGrainCreated(unsigned int new_grain_id)
{
  if (_euler_angles_vec.size() <= new_grain_id)
  {
    _euler_angles_vec.resize(new_grain_id + 1);
    _C_ijkl_vec.resize(new_grain_id + 1);
  }

  _euler_angles_vec[new_grain_id] = newGrainEuler(new_grain_id);
  _C_ijkl_vec[new_grain_id] = newGrainElastic(new_grain_id);
}

EulerAngles
GrainTrackerMatProp::newGrainEuler(unsigned int new_grain_id)
{
  EulerAngles angles;

  if (new_grain_id < _euler.getGrainNum())
    angles = _euler.getEulerAngles(new_grain_id);
  else
  {
    if (_random_rotations)
      angles.random();
    else
      mooseError("GrainTrackerMatProp has run out of grain rotation data.");
  }
  return angles;
}

RankFourTensor
GrainTrackerMatProp::newGrainElastic(unsigned int new_grain_id)
{
  EulerAngles angles = _euler_angles_vec[new_grain_id];

  RankFourTensor C_ijkl = _C_ijkl;
  C_ijkl.rotate(RotationTensor(RealVectorValue(angles)));

  return C_ijkl;  
}