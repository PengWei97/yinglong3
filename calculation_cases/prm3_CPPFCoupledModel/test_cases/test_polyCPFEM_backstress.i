my_d_bs = 1.0e2
my_c_bs = 1.0

my_filename = "t2_n200_tensile_v1"
my_filename2 = "t2_n200_tensile_v1"

[Functions]
  [./cycle_load]
    type = ParsedFunction
    expression = '0.1*t' # 0.1*t 0.1*sin(2*pi*t) 0.2*sin(2*pi*t)
  [../]
[]

[Mesh]
  [./nepermesh]
    type = FileMeshGenerator
    file = n200-id1_v2.msh
  [../]
  [./left_modifier]
    type = BoundingBoxNodeSetGenerator
    input = nepermesh
    new_boundary = Left
    top_right = '0.1 100.1 0'
    bottom_left = '-0.1 -0.1 0'
  [../]
  [./bottom_modifier]
    type = BoundingBoxNodeSetGenerator
    input = left_modifier
    new_boundary = Bottom
    top_right = '100.1 0.1 0'
    bottom_left = '-0.1 -0.1 0'
  [../]
  [./top_modifier]
    type = BoundingBoxNodeSetGenerator
    input = bottom_modifier
    new_boundary = Top
    top_right = '100.1 100.1 0.0'
    bottom_left = '-0.1 99.9 0.0'
  [../]    
  parallel_type = distributed
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[UserObjects]
  [./prop_read]
    type = GrainPropertyReadFileCP
    prop_file_name = euler_ang_test_200.inp
    # Enter file data as prop#1, prop#2, .., prop#nprop
    nprop = 3
	  ngrain = 200
    read_type = indexgrain
  [../]
  [./term]
    type = Terminator
    expression = 'e_yy > 0.20'
  [../]
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Modules/TensorMechanics/Master/all]
  strain = FINITE
[]

[AuxVariables]
  [./pk2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fp_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_increment]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./backstress]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./euler1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./euler2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./euler3]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./crysrot11]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./pk2]
   type = RankTwoAux
   variable = pk2
   rank_two_tensor = second_piola_kirchhoff_stress
   index_j = 1
   index_i = 1
   execute_on = timestep_end
  [../]
  [./fp_yy]
    type = RankTwoAux
    variable = fp_yy
    rank_two_tensor = plastic_deformation_gradient
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]
  [./e_yy]
    type = RankTwoAux
    variable = e_yy
    rank_two_tensor = total_lagrangian_strain
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]
  [./gss]
   type = MaterialStdVectorAux
   variable = gss
   property = slip_resistance
   index = 0
   execute_on = timestep_end
  [../]
  [./slip_inc]
   type = MaterialStdVectorAux
   variable = slip_increment
   property = slip_increment
   index = 0
   execute_on = timestep_end
  [../]
  [./backstress]
   type = MaterialStdVectorAux
   variable = backstress
   property = backstress
   index = 0
   execute_on = timestep_end
  [../]
  [./stress_yy]
    type = RankTwoAux
    variable = stress_yy
    rank_two_tensor = stress
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]  
  [./stress_xx]
    type = RankTwoAux
    variable = stress_xx
    rank_two_tensor = stress
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  [../]
  [./euler1]
    type = MaterialRealVectorValueAux
    variable = euler1
    property = Euler_angles
    component = 0
    execute_on = timestep_end
  [../]
  [./euler2]
    type = MaterialRealVectorValueAux
    variable = euler2
    property = Euler_angles
    component = 1
    execute_on = timestep_end
  [../]
  [./euler3]
    type = MaterialRealVectorValueAux
    variable = euler3
    property = Euler_angles
    component = 2
    execute_on = timestep_end
  [../]
  [./crysrot11]
    type = RankTwoAux
    variable = crysrot11
    rank_two_tensor = crysrot
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  [../]
[]

[BCs]
  [./symmy]
    type = DirichletBC
    variable = disp_y
    boundary = Bottom
    value = 0
  [../]
  [./symmx]
    type = DirichletBC
    variable = disp_x
    boundary = Left
    value = 0
  [../]
  [./tdisp]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = Top
    function = cycle_load
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCPGrain
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.754e5 0.754e5 0.754e5'
    fill_method = symmetric9
    read_prop_user_object = prop_read
  [../]
  [./stress]
    type = ComputeMultipleCrystalPlasticityStress
    crystal_plasticity_models = 'trial_xtalpl'
    tan_mod_type = exact
    rtol = 1e-6 # Constitutive stress residual relative tolerance
    maxiter_state_variable = 50 # Maximum number of iterations for stress update

    maximum_substep_iteration = 25 # Maximum number of substep iteration

    use_line_search = true
  [../]
  [./trial_xtalpl]
    type = CPKalidindiBackstressUpdate
    crystal_lattice_type = FCC # BCC FCC HCP
    number_slip_systems = 12 
    slip_sys_file_name = input_slip_sys.txt
    number_cross_slip_directions = 0
    number_cross_slip_planes = 0

    slip_increment_tolerance = 0.1 # 2e-2, Maximum allowable slip in an increment for each individual constitutive model
    stol = 0.1 # Constitutive internal state variable relative change tolerance
    resistance_tol = 0.1 # 2e-2, Con
    zero_tol = 1e-12 # Tolerance for residual check when variable value is zero for each individual constitutive model
    print_state_variable_convergence_error_messages = false

    # https://www.sciencedirect.com/science/article/pii/S0921509317300898
    c_bs = ${my_c_bs} # 1.0e1 # 7.56e4 1e4
    d_bs = ${my_d_bs} # 1.0e-1 # 7.20e2 1e2
  [../]
[]

[Postprocessors]
  [./dofs]
    type = NumDOFs
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
  [./stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  [../]
  [./pk2]
   type = ElementAverageValue
   variable = pk2
  [../]
  [./fp_yy]
    type = ElementAverageValue
    variable = fp_yy
  [../]
  [./e_yy]
    type = ElementAverageValue
    variable = e_yy
  [../]
  [./gss]
    type = ElementAverageValue
    variable = gss
  [../]
  [./slip_increment]
   type = ElementAverageValue
   variable = slip_increment
  [../]
  [./backstress]
    type = ElementAverageValue
    variable = backstress
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -ksp_type -ksp_gmres_restart'
  petsc_options_value = ' asm      2              lu            gmres     200'
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-10
  nl_abs_step_tol = 1e-10

  # l_max_its = 30 # Max number of linear iterations
  nl_max_its = 20 # Max number of nonlinear iterations

  start_time = 0.0
  end_time = 300
  # num_steps = 5
  # dt = 0.1
  dtmin = 0.1e-6
  dtmax = 5.0

  [./TimeStepper] 
    type = IterationAdaptiveDT
    dt = 0.01 # Initial time step.  In this simulation it changes.
    optimal_iterations = 30 # Time step will adapt to maintain this number of nonlinear iterations
    iteration_window = 5
  [../]
  # [./Adaptivity] # 加上网格自适应，计算很慢
  #   # Block that turns on mesh adaptivity. Note that mesh will never coarsen beyond initial mesh (before uniform refinement)
  #   initial_adaptivity = 2 # Number of times mesh is adapted to initial condition
  #   refine_fraction = 0.7 # Fraction of high error that will be refined
  #   coarsen_fraction = 0.1 # Fraction of low error that will coarsened
  #   max_h_level = 4 # Max number of refinements used, starting from initial mesh (before uniform refinement)
  # [../]
[]

[Outputs]
  [./my_checkpoint]
    file_base = ./${my_filename}/out_${my_filename}
    interval = 5
    type = Checkpoint
    additional_execute_on = 'FINAL'
  [../]
  [my_exodus]
    file_base = ./ex_${my_filename2}/out_${my_filename} 
    interval = 10
    type = Nemesis
    additional_execute_on = 'FINAL'
  [../]
  [./csv]
    file_base = ./csv_${my_filename}/out_${my_filename}
    # interval = 5
    type = CSV
  [../]
  print_linear_residuals = false
[]