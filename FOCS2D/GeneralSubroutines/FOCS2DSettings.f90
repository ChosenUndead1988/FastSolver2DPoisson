module FOCS2DSettings

  use FOCS2DSpecifyPrecision

  logical, parameter :: MULTICOLOR_SMOOTHER = .true.
  integer, parameter :: NUMBER_AUXILIARY_GRIDS = 5
  integer, parameter :: LEADING_FACTOR = 1
  integer, parameter :: ORDER_RESTRICTION_OPERATOR = 2
  integer, parameter :: ORDER_INTERPOLATION_OPERATOR = 4
  integer, parameter :: ORDER_INTERPOLATION_NESTED_ITERATION = 4
  integer, parameter :: NUM_SWEEPS_DESCENT_MULTIGRID = 1
  integer, parameter :: NUM_SWEEPS_ASCENT_MULTIGRID = 1
  integer, parameter :: NUM_CYCLES_NESTED_ITERATION = 3
  integer, parameter :: MAXIMUM_CYCLES_MULTIGRID = 20
  character(len=1), parameter :: CYCLE_TYPE_MULTIGRID = 'V'
  real(dp), parameter :: ABSOLUTE_TOLERANCE_MULTIGRID = 10.0_dp**(-12)
  real(dp), parameter :: RELATIVE_TOLERANCE_MULTIGRID = 10.0_dp**(-12)

end module FOCS2DSettings
