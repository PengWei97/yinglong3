# chaoSuan job

```bash
#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 2
#SBATCH -n 256

mpiexec -n 256 ~/projects/qinglong/qinglong-opt -i GNSoTi_AGG_level2a3_700du.i

#  --recover
# ./case4_recovery_v3/out_case4_recovery_v3_cp/2201
```

# Linux

## 打包
```bash
tar -cvf - ex_case4_recovery_v92/*.e-s0??[2].* ex_case4_recovery_v92/*.e-s0115.* | pigz -9 -p 20 > ex_case4_recovery_v92.tgz 
tar -cvf - ex_case4_recovery_v92/*.e-s0070.* ex_case4_recovery_v92/*.e-s0101.* | pigz -9 -p 20 > ex_case4_recovery_v92.tgz 
ll ex_case4_recovery_v92/*.e-s*.056
```

## 创建

code /home/pw-moose/projects/baize/include/materials/crystal_plasticity/CPKalidindiBackstressUpdate.h
code /home/pw-moose/projects/baize/src/materials/crystal_plasticity/CPKalidindiBackstressUpdate.C
code /home/pw-moose/projects/baize/doc/content/source/materials/crystal_plasticity/CPKalidindiBackstressUpdate.md

mkdir /home/pw-moose/projects/baize/include/userobjects/
mkdir /home/pw-moose/projects/baize/src/userobjects/
mkdir -p /home/pw-moose/projects/baize/doc/content/source/materials/crystal_plasticity/

cp /home/pw-moose/projects/c_pfor_am/src/materials/ComputeElasticityTensorCPGrain.C /home/pw-moose/projects/baize/src/materials/ComputeElasticityTensorCPGrain.C
cp /home/pw-moose/projects/c_pfor_am/include/materials/ComputeElasticityTensorCPGrain.h /home/pw-moose/projects/baize/include/materials/ComputeElasticityTensorCPGrain.h

# git 开发
```bash
mamba install git-flow
```

## Git flow流程开发
```bash
git checkout develop
git flow feature start add-BackstressCPFEM
...
git checkout develop
git pull
git flow release start '3.0.0'
git flow release finish '3.0.0'
git tag # develop
git push origin develop
git chechout main
git push origin main
git push --tags
```

