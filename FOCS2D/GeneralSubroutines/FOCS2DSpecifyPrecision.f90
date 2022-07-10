module FOCS2DSpecifyPrecision

  implicit none

  !> Single precision real numbers, 6 digits, range 10⁻³⁷ to 10³⁷-1; 32 bits
  integer, parameter :: sp = selected_real_kind(6, 37)
  !> Double precision real numbers, 15 digits, range 10⁻³⁰⁷ to 10³⁰⁷-1; 64 bits
  integer, parameter :: dp = selected_real_kind(15, 307)
  !> Quadruple precision real numbers, 33 digits, range 10⁻⁴⁹³¹ to 10⁴⁹³¹-1; 128 bits
  integer, parameter :: qp = selected_real_kind(33, 4931)

end module FOCS2DSpecifyPrecision
