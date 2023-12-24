# 这是创建耦合模型的基本测试算例

my_filename = "t0_gg_elastic_fatigue_bi"
my_filename2 = "t0_gg_elastic_fatigue_bi"

[Materials]
  [./Copper]
    type = GBEvolution
    block = 0
    T = 500 # K
    wGB = 2 # nm
    GBmob0 = 2.5e-6 # m^4/(Js) from Schoenfelder 1997
    Q = 0.23 # Migration energy in eV
    # GBMobility = 0.0
    GBenergy = 0.708 # GB energy in J/m^2

    # length_scale = 1.0e-9
    # time_scale = 1.0
  [../]
  [./ElasticityTensor]
    type = ComputePolyElasticTensorCpl
    grain_tracker = grain_tracker
  [../]
  [./stress]
    type = ComputeMultCPStressCpl
    crystal_plasticity_models = 'trial_xtalpl'
    tan_mod_type = exact
    use_line_search = true

    # rtol = 1e-6 # Constitutive stress residual relative tolerance
    # maxiter_state_variable = 50 # Maximum number of iterations for stress update
    # maximum_substep_iteration = 25 # Maximum number of substep iteration
  [../]
  [./trial_xtalpl]
    type = CPKalidindiBSCplUpdate
    crystal_lattice_type = FCC
    number_slip_systems = 12 
    slip_sys_file_name = input_slip_sys.txt

    ao = 0.001 # 
    gss_initial = 30.8
    t_sat = 148
    # h = 180
    slip_increment_tolerance = 0.1 # Maximum allowable slip in an increment
    stol = 0.1 # Constitutive internal state variable relative change tolerance
    resistance_tol = 0.1

    grain_tracker = grain_tracker
    vpp_name = grain_volumes

    c_bs = 0.0 # 1.0e4
    d_bs = 0.0 # 1.0e2

    output_properties = 'slip_resistance'
    outputs = my_exodus
  [../]
  [./strain]
    type = ComputeFiniteStrain
    block = 0
    displacements = 'disp_x disp_y'
  [../]
[]

[Kernels]
  [./PolycrystalKernel]
    # ACGrGrPoly ACInterface TimeDerivative
  [../]
  [./PolycrystalElasticDrivingForce]
    # ACGrGrElasticDrivingForces
  [../]
  [./TensorMechanics]
    use_displaced_mesh = true
    displacements = 'disp_x disp_y'
  [../]
[]

[UserObjects]
  [./euler_angle_file]
    type = EulerAngleFileReader
    file_name = grn_36_rand_2D.tex
  [../]
  [./grain_tracker]
    type = GrainTrackerMatProp # GrainTrackerElasticity
    threshold = 0.2
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_begin'

    C_ijkl = '1.27e5 0.708e5 0.708e5 1.27e5 0.708e5 1.27e5 0.7355e5 0.7355e5 0.7355e5'
    fill_method = symmetric9

    flood_entity_type = ELEMENTAL
    euler_angle_provider = euler_angle_file
  [../]
  # [./term]
  #   type = Terminator
  #   expression = 'run_time > 7200'
  # [../]    
[]

[VectorPostprocessors]
  [./grain_volumes]
    type = FeatureMatPropCplVectorPostprocessor
    flood_counter = grain_tracker

    number_slip_systems = 12
    # mat_prop = elastic_energy
    execute_on = 'INITIAL TIMESTEP_BEGIN' 
    # 'INITIAL TIMESTEP_BEGIN LINEAR NONLINEAR'
    
    # TIMESTEP_BEGIN ~ Must exist, otherwise an error will be reported
  [../]
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 20
  xmax = 100
  ymax = 100
  elem_type = QUAD4
  uniform_refine = 2
[]

[GlobalParams]
  op_num = 2
  var_name_base = gr

  displacements = 'disp_x disp_y'

  length_scale = 1.0e-9
  pressure_scale = 1.0e6
[]

[Variables]
  [./PolycrystalVariables] # action
  [../]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./BicrystalBoundingBoxIC]
      x1 = 0
      y1 = 0
      x2 = 50
      y2 = 100
    [../]
  [../]
[]

[AuxVariables]
  [./pk2]
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
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  [./elastic_strain22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./unique_grains]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./var_indices]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./C1111]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./euler_angle]
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
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
  [./elastic_strain22]
    type = RankTwoAux
    variable = elastic_strain22
    rank_two_tensor = total_lagrangian_strain
    index_i = 1
    index_j = 1
    execute_on = timestep_end
  [../]
  [./stress22]
    type = RankTwoAux
    variable = stress22
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    execute_on = timestep_end
  [../]
  [./unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    execute_on = timestep_end
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
  [../]
  [./var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    execute_on = timestep_end
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
  [../]
  [./C1111]
    type = RankFourAux
    variable = C1111
    rank_four_tensor = elasticity_tensor
    index_l = 0
    index_j = 0
    index_k = 0
    index_i = 0
    execute_on = timestep_end
  [../]
  [./vonmises_stress]
    type = RankTwoScalarAux
    variable = vonmises_stress
    rank_two_tensor = stress
    scalar_type = VonMisesStress
    execute_on = timestep_end
  [../]
  [./euler_angle]
    type = OutputEulerAngles
    variable = euler_angle
    euler_angle_provider = euler_angle_file
    grain_tracker = grain_tracker
    output_euler_angle = 'phi1'
    execute_on = 'initial timestep_end'
  [../]
[]

[Functions]
  [./tdisp]
    type = ParsedFunction
    expression = '0.5*sin(2*pi*t)' # fatigue
    # expression = 'if (t<10,0.1*t,1.0)' # 'if (t<90,0.1*t,9.0)'
  [../]
[]

[BCs]
  # [./Periodic]
  #   [./All]
  #     auto_direction = 'x'
  #     # variable = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8'
  #     variable = 'gr0 gr1 gr2 gr3 gr4'
  #   [../]
  # [../]
  [./top_displacement]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = tdisp
  [../]
  [./x_anchor]
    type = DirichletBC
    variable = disp_x
    boundary = left # 'left right'
    value = 0.0
  [../]
  [./y_anchor]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
[]

[Postprocessors]
  [./elastic_strain22]
    type = ElementAverageValue
    variable = elastic_strain22
  [../]
  [./stress22]
    type = ElementAverageValue
    variable = stress22
  [../]
  [./pk2]
   type = ElementAverageValue
   variable = pk2
  [../]
  [./ngrains]
    type = FeatureFloodCount
    variable = bnds
    threshold = 0.7
  [../]
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
  [./gss]
    type = ElementAverageValue
    variable = gss
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    coupled_groups = 'disp_x,disp_y'
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 31 0.7'
  l_tol = 1.0e-4
  l_max_its = 30
  nl_max_its = 25
  nl_rel_tol = 1.0e-7
  start_time = 0.0
  # dt = 0.001
  dtmax = 0.2
  end_time = 100
  # num_steps = 10
  
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.05
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    initial_adaptivity = 3
    refine_fraction = 0.8
    coarsen_fraction = 0.05
    max_h_level = 3
  [../]
[]

[Outputs]
  [./my_checkpoint]
    file_base = ./${my_filename}/out_${my_filename}
    interval = 10
    type = Checkpoint
    additional_execute_on = 'INITIAL FINAL'
  [../]
  [my_exodus]
    file_base = ./ex_${my_filename2}/out_${my_filename} 
    # interval = 100
    type = Nemesis
    additional_execute_on = 'FINAL'
    sync_only = true
    sync_times = '1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	43	44	45	46	47	48	49	50	51	52	53	54	55	56	57	58	59	60	61	62	63	64	65	66	67	68	69	70	71	72	73	74	75	76	77	78	79	80	81	82	83	84	85	86	87	88	89	90	91	92	93	94	95	96	97	98	99	100'
  [../]
  [./csv]
    file_base = ./csv_${my_filename}/out_${my_filename}
    interval = 10
    type = CSV
  [../]
  [pgraph]
    type = PerfGraphOutput
    interval = 10
    level = 2                     # Default is 1
    heaviest_branch = false        # Default is false
    heaviest_sections = 7         # Default is 0
    execute_on = 'TIMESTEP_END'
    additional_execute_on = 'FINAL'  # Default is "final"
  []
  print_linear_residuals = false
[]
