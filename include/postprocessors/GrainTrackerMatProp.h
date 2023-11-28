#pragma once

// Purpose: manage the map from [grain id -> material properties] based on GrainTracker
// material properties: euler angles, elastic tensor

#include "GrainTracker.h"
#include "RankFourTensor.h"
#include "EulerAngles.h"

class EulerAngleProvider;

class GrainTrackerMatProp : public GrainTracker
{
public:
  static InputParameters validParams();

  GrainTrackerMatProp(const InputParameters & parameters);

  const EulerAngles & getEulerAngles(unsigned int grain_id) const;
  const RankFourTensor & getElasticTensor(unsigned int grain_id) const;

protected:
  virtual void newGrainCreated(unsigned int new_grain_id);

  virtual EulerAngles newGrainEuler(unsigned int new_grain_id);

  virtual RankFourTensor newGrainElastic(unsigned int new_grain_id);

  std::vector<EulerAngles> & _euler_angles_vec;

  std::vector<RankFourTensor> & _C_ijkl_vec;

  /// generate random rotations when the Euler Angle provider runs out of data (otherwise error out)
  const bool _random_rotations;

  /// unrotated elasticity tensor
  RankFourTensor _C_ijkl;

  /// object providing the Euler angles
  const EulerAngleProvider & _euler;
};
