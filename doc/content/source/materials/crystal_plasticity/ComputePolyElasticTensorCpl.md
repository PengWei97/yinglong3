# ComputePolyElasticTensorCpl

!syntax description /Materials/crystal_plasticity/ComputePolyElasticTensorCpl

This material generates an elastic tensor in grain growth phase-field and crystal 
plasticity coupled model.

## Description

Drawing on the interpolation processing method of `ComputePolycrystalElasticityTensor`, obtain the ordinal parameter value corresponding to each [_qp] through `GrainTrackerMatProp`. Then use the interpolation function h to perform interpolation processing to calculate `_elasticity_tensor`.

The rotation matrix used in this class,`ComputePolyElasticTensorCpl`, is the transpose of the rotation
matrix created from the Bunge Euler angles in the base class. This difference in the rotation matrix is because of the active rotation convention used in crystal plasticity simulations.

## Rotation Tensor Conventions

The [Euler angle convention](http://mathworld.wolfram.com/EulerAngles.html) used in
`ComputePolyElasticTensorCpl` is the $z$-$x'$-$z'$ (3-1-3) convention.  The Euler angles arguments are
expected in degrees, not radians, and are denoted as $\phi_1$, $\Phi$, and $\phi_2$, corresponding to
the axis rotations.  

The rotation matrix is obtained by finding the euler angles corresponding to the grain_id corresponding to the maximum order parameter at [_qp], which does not consider the issue of grain rotation during plastic deformation.

## Example Input File Syntax

<!-- !listing modules/tensor_mechanics/test/tests/crystal_plasticity/stress_update_material_based/update_method_test.i block=Materials/elasticity_tensor

!syntax parameters /Materials/ComputeElasticityTensorCP

!syntax inputs /Materials/ComputeElasticityTensorCP

!syntax children /Materials/ComputeElasticityTensorCP

!bibtex bibliography -->
