#use miug mius mm
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

[GlobalParams]
  implicit = false
[]

#[Problem]
  #Note that the suffix is left off in the parameter below.
  #restart_file_base = displacement_t_out_cp/LATEST  # You may also use a specific number here
#[]

[AuxVariables]
  [shearband_driving_force]
    order = CONSTANT
    family = MONOMIAL
  []
  [bounds_dummy]
  []
[]

[Bounds]
  [irreversibility]
    type = VariableOldValueBoundsAux
    variable = bounds_dummy
    bounded_variable = c
    bound_type = lower
  []
  [upper]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = c
    bound_type = upper
    bound_value = 1.0
  []
[]

# [ICs]
#   [pre]
#     type = BoundingBoxIC
#     x1 = -0.0044
#     y1 = 1.75
#     x2 = 0.0044
#     y2 = 1.7715
#     variable = c
#     inside = 0.999
#     outside = 0
#   []
# []

#[Modules]
  #[PhaseField]
   # [Nonconserved]
     # [c]
      #  free_energy = F#refers to floc
       # kappa = kappa_op
       # mobility = L
    #  []
   # []
 # []
#[]

[Variables]
  [c]
  []
[]

[Kernels]
  [trans]
    type = TimeDerivative
    lumping = true
    variable = c
    implicit = false
  []
  [Bulk]
    type = AllenCahn
    variable = c
    f_name = F
    mob_name = L
  []
  [ACInterface]
    type = ACInterface
    variable = c
    kappa_name = kappa_op
    mob_name = L
  []
[]

[Materials]
  [pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l visco'
    prop_values = '230.90706	0.137444679	0.036378564' # *2/12 before:/8
    # prop_values = '14.448 0.86e-2 0.5814'
    # prop_values = '3e-3 0.25 2.2'
    #16.8kj/m^2 = 16.8 1e3 N*m/m^2 = 16.8 1e3 N/m = 16.8 10e3 kg/s^2 = 16.8 miug/mius^2
    #10 mium = 1e-2 mm
    #yita~(zhang) = yita* c_a * l/gc = 0.01 mius, c_a = 2
    #yita(origin) = 0.01*gc/(c_a*l) miug/(mm * mius)
    #visco(moose) = yita(origin)/gc = yita~(zhang)/(c_a*l) = 0.5 mius/mm
  []
  [define_mobility]
    type = ParsedMaterial
    material_property_names = 'gc_prop visco'
    property_name = L
    expression = '1.0/(gc_prop * visco)'
  []
  [define_kappa]
    type = ParsedMaterial
    material_property_names = 'gc_prop l'
    property_name = kappa_op
    expression = 'gc_prop * l'#kappa = 2gc * l/c_a), c_a = 2
  []
  [local_fracture_energy]#c_a = 2, a(d) = d^2, f_frac = (1/c_a) * (gc/l) * a(d)
    type = DerivativeParsedMaterial
    property_name = local_fracture_energy
    coupled_variables = 'c'
    material_property_names = 'gc_prop l'
    expression = '(c^2 * gc_prop)/(2 * l)'
    derivative_order = 2
  []
  [./fracture_driving_energy]
    type = DerivativeParsedMaterial
    coupled_variables = 'c shearband_driving_force'
    material_property_names = 'local_fracture_energy(c) degradation(c)'
    expression = 'local_fracture_energy + degradation * shearband_driving_force'
    derivative_order = 2
    property_name = F
  [../]
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

#[Preconditioning]
  #[spm]
    #type = SMP
  #  full=true
#  []
#[]

[Executioner]
  type = Transient

  #solve_type = PJFNK
  #petsc_options_iname = '-pc_type  -snes_type'
  #petsc_options_value = 'lu vinewtonrsls'

   petsc_options_iname = '-pc_type -pc_factor_shift  -snes_type'
   petsc_options_value = 'lu NONZERO vinewtonrsls'
  [TimeIntegrator]
    type = CentralDifference
    solve_type = lumped
  []
 #nl_rel_tol = 1e-8
 #nl_abs_tol = 1e-10
 #l_max_its = 100
 #nl_max_its = 100
[]