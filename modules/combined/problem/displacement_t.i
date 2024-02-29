#use miug mius mm
[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = fracture.i
    execute_on = TIMESTEP_END
  []
[]

[Transfers]
  [from_c]
    type = MultiAppCopyTransfer
    from_multi_app = fracture
    variable = c
    source_variable = c
  []
  [to_shearband_driving_force]
    type = MultiAppCopyTransfer
    to_multi_app = fracture
    variable = shearband_driving_force
    source_variable = shearband_driving_force
  []
[]

[Mesh]
  allow_renumbering = false
  [TWC]
    type = AnnularMeshGenerator
    rmax = 5
    rmin = 3.5
    nt = 320
    nr = 21
  []
[]

#[Problem]
  #Note that the suffix is left off in the parameter below.
  #restart_file_base = displacement_t_out_cp/LATEST  # You may also use a specific number here
#[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
[]

[AuxVariables]
  [c]
  []
  [shearband_driving_force]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = MaterialRealAux
      property = local_energy_dev
    []
  []
  [vonmises]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = RankTwoScalarAux
      rank_two_tensor = stress
      variable = vonmises
      scalar_type = VonMisesStress
    []
  []
  [pressure]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = RankTwoScalarAux
      rank_two_tensor = stress
      variable = pressure
      scalar_type = Hydrostatic
    []
  []
  [EPS]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = MaterialRealAux
      property = eqv_plastic_strain
      variable = EPS
    []
  []
  [yield_factor]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = MaterialRealAux
      property = yield_factor
      variable = yield_factor
    []
  []
  [accel_x]
  []
  [vel_x]
  []
  [accel_y]
  []
  [vel_y]
  []
[]


[AuxKernels]
  [accel_x]
      type = TestNewmarkTI
      variable = accel_x
      displacement = disp_x
      first = false
  []
  [vel_x]
      type = TestNewmarkTI
      variable = vel_x
      displacement = disp_x
  []
  [accel_y]
      type = TestNewmarkTI
      variable = accel_y
      displacement = disp_y
      first = false
  []
  [vel_y]
      type = TestNewmarkTI
      variable = vel_y
      displacement = disp_y
  []
[]

[Kernels]
  [DynamicTensorMechanics]
    displacements = 'disp_x disp_y'
    strain = FINITE
    incremental = true
    decomposition_method = HughesWinget
  []
  [inertia_x]
    type = InertialForce
    variable = disp_x
  []
  [inertia_y]
    type = InertialForce
    variable = disp_y
  []   
[]

[BCs]
  [Pressure]
      [outter_surface]
          boundary = 'rmax'
          function = triangle_pulse
          displacements = 'disp_x disp_y'
      []
  []
[]

[Functions]
  [triangle_pulse]
      #4Gpa=4e9n/m^2=4e9 kg/(m*s^2)=4e9 kg/(e3e12mm*mius^2)=4e3 miug/(mm*mius^2)
      type = PiecewiseLinear
      x = '0	0.1	0.2	0.3	0.4	0.5	0.6	0.7	0.8	0.9	1	1.1	1.2	1.3	1.4	1.5	1.6	1.7	1.8	1.9	2	2.1	2.2	2.3	2.4	2.5	2.6	2.7	2.8	2.9	3	3.1	3.2	3.3	3.4	3.5	3.6	3.7	3.8	3.9	4	4.1	4.2	4.3	4.4	4.5	4.6	4.7	4.8	4.9	5	5.1	5.2	5.3	5.4	5.5	5.6	5.7	5.8	5.9	6	6.1	6.2	6.3	6.4	6.5	6.6	6.7	6.8	6.9	7	7.1	7.2	7.3	7.4	7.5	7.6	7.7	7.8	7.9	8	8.1	8.2	8.3	8.4	8.5	8.6	8.7	8.8	8.9	9	9.1	9.2	9.3	9.4	9.5	9.6	9.7	9.8	9.9	10	10.1	10.2	10.3	10.4	10.5	10.6	10.7	10.8	10.9	11	11.1	11.2	11.3	11.4	11.5	11.6	11.7	11.8	11.9	12	12.1	12.2	12.3	12.4	12.5	12.6	12.7	12.8	12.9	13	13.1	13.2	13.3	13.4	13.5	13.6	13.7	13.8	13.9	14	14.1	14.2	14.3	14.4	14.5	14.6	14.7	14.8	14.9	15	15.1	15.2	15.3	15.4	15.5	15.6	15.7	15.8	15.9	16	16.1	16.2	16.3	16.4	16.5	16.6	16.7	16.8	16.9	17	17.1	17.2	17.3	17.4	17.5	17.6	17.7	17.8	17.9	18	18.1	18.2	18.3	18.4	18.5	18.6	18.7	18.8	18.9	19	19.1	19.2	19.3	19.4	19.5	19.6	19.7	19.8	19.9	20
'
      y = '

4.02E-02	4.03E-02	3.96E-02	3.93E-02	2.76E-02	-8.19E-03	-8.96E-02	-2.25E-01	-4.23E-01	-6.42E-01	-8.92E-01	-1.06E+00	-1.00E+00	-5.29E-01	4.22E-01	2.00E+00	4.27E+00	7.34E+00	1.12E+01	1.59E+01	2.15E+01	2.80E+01	3.54E+01	4.36E+01	5.25E+01	6.20E+01	7.22E+01	8.31E+01	9.45E+01	1.07E+02	1.19E+02	1.33E+02	1.47E+02	1.61E+02	1.76E+02	1.91E+02	2.07E+02	2.23E+02	2.40E+02	2.56E+02	2.73E+02	2.90E+02	3.07E+02	3.24E+02	3.41E+02	3.59E+02	3.76E+02	3.94E+02	4.11E+02	4.28E+02	4.46E+02	4.63E+02	4.81E+02	4.99E+02	5.17E+02	5.34E+02	5.52E+02	5.69E+02	5.86E+02	6.02E+02	6.19E+02	6.35E+02	6.51E+02	6.67E+02	6.83E+02	6.98E+02	7.13E+02	7.27E+02	7.41E+02	7.55E+02	7.69E+02	7.82E+02	7.95E+02	8.07E+02	8.18E+02	8.29E+02	8.39E+02	8.48E+02	8.56E+02	8.64E+02	8.70E+02	8.75E+02	8.80E+02	8.84E+02	8.87E+02	8.89E+02	8.90E+02	8.90E+02	8.89E+02	8.87E+02	8.84E+02	8.80E+02	8.75E+02	8.68E+02	8.60E+02	8.52E+02	8.43E+02	8.32E+02	8.21E+02	8.09E+02	7.96E+02	7.83E+02	7.69E+02	7.55E+02	7.40E+02	7.26E+02	7.12E+02	6.97E+02	6.82E+02	6.67E+02	6.50E+02	6.33E+02	6.15E+02	5.97E+02	5.79E+02	5.61E+02	5.42E+02	5.23E+02	5.03E+02	4.83E+02	4.64E+02	4.45E+02	4.26E+02	4.07E+02	3.88E+02	3.70E+02	3.52E+02	3.35E+02	3.17E+02	3.00E+02	2.82E+02	2.65E+02	2.49E+02	2.33E+02	2.17E+02	2.02E+02	1.87E+02	1.72E+02	1.59E+02	1.45E+02	1.32E+02	1.20E+02	1.08E+02	9.71E+01	8.63E+01	7.61E+01	6.63E+01	5.69E+01	4.83E+01	4.03E+01	3.31E+01	2.66E+01	2.09E+01	1.59E+01	1.15E+01	7.91E+00	4.95E+00	2.66E+00	9.90E-01	4.47E-02	-2.36E-01	2.15E-01	1.34E+00	3.22E+00	5.72E+00	8.89E+00	1.27E+01	1.72E+01	2.23E+01	2.81E+01	3.47E+01	4.18E+01	4.97E+01	5.81E+01	6.73E+01	7.71E+01	8.72E+01	9.79E+01	1.09E+02	1.20E+02	1.32E+02	1.44E+02	1.57E+02	1.70E+02	1.83E+02	1.96E+02	2.10E+02	2.24E+02	2.38E+02	2.52E+02	2.67E+02	2.82E+02	2.97E+02	3.12E+02	3.27E+02	3.42E+02	3.56E+02	3.71E+02	3.85E+02	3.99E+02	4.13E+02

'
  []
[]

[Materials]
  [density]
    type = Density
    displacements = 'disp_x disp_y'
    density = 7.9e3 #miug/mm^3
  []
  [strain]
    type = ComputeFiniteStrain
    displacements = 'disp_x disp_y'
    implicit = false
    decomposition_method = HughesWinget
  []
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    shear_modulus = 77.5e3  #77.5Gpa
    poissons_ratio = 0.3
  []
  [damage_stress]
    type = JCEOSnl
    A = 110#Mpa
    B = 1500
    n = 0.36
    C = 0.025
    reference_strain_rate = 1
    m = 0.81
    T_m = 1800#kelvin
    T_r = 298
    T_i = 298
    SpecHeat = 500e-6 #1e-6mm^2/(mius^2*kelvin)=1*m*m/(s*s*kelvin)
    D_i = 7.9e3
    ShearModule = 77.5e3
    C0 = 4569e-3 #1e-3 mm/mius = 1m/s
    S1 = 1.49
    gamma0 = 2.17
    a = 0
    read_prop_user_object = prop_read
    c = c
    D_name = 'degradation'
    shearband_driving_force_name = local_energy_dev
  []
  [degradation]#a1
    type = DerivativeParsedMaterial
    property_name = degradation
    coupled_variables = 'c'
    expression = '((1.0-c)^2)/(((1.0-c)^2)+a1*c*(1-0.5*c))'
    constant_names       = 'a1'
    constant_expressions = '0.8'
    derivative_order = 2
  []
[]

[UserObjects]
  [prop_read]
    type = PropertyReadFile
    prop_file_name = '0.02722-1.txt'
    nprop = 1
    read_type = element
  []
[]

[Executioner]
  type = Transient
  start_time = 0
  end_time = 20
  dt = 0.005 # local speed is 6000m/s=6000mium/mius——1e-4mius for 0.6mium
   petsc_options_iname = '-pc_type -pc_factor_shift  -snes_type'
   petsc_options_value = 'lu NONZERO vinewtonrsls'
  [TimeIntegrator]
    type = CentralDifference
    solve_type = lumped
  []
[]

[Outputs]
  exodus = true
  interval = 20
  checkpoint = true
[]