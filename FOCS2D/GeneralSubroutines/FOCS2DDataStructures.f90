module FOCS2DDataStructures

  use FOCS2DSpecifyPrecision

  type MultigridStructure

    integer :: gridPointsPerAxisX
    integer :: gridPointsPerAxisY
    real(dp) :: stepSizeX
    real(dp) :: stepSizeY
    real(dp), dimension(:,:), allocatable :: RHS
    real(dp), dimension(:,:), allocatable :: approximateSolution

   end type MultigridStructure

end module FOCS2DDataStructures
