module FOCS2DConstants

  use FOCS2DSpecifyPrecision

  implicit none

  !> The real number pi. Since arctan( 1.0 ) = pi/4.0
  real(kind=dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
  !> The smallest number x such tha 1 + x > 1
  real(kind=dp), parameter :: machineEps = epsilon(1.0_dp)
  !> zero prerecorded
  real(kind=dp), parameter :: ZERO = 0.0_dp

end module FOCS2DConstants
