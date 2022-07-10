module FOCS2DMatvecs

  use FOCS2DSpecifyPrecision
  use FOCS2DConstants
  use FOCS2DMiscellaneous

  private
  public :: getResidualFOCS, applyLHSOperatorFOCS, applyRHSOperatorFOCS

  contains

  subroutine getResidualFOCS(Nx, Ny, stepSizeX, stepSizeY, &
    RHS, gridFunction, residual)
    integer, intent(in) :: Nx, Ny
    real(dp), intent(in) :: stepSizeX, stepSizeY, RHS(0:,0:), gridFunction(0:,0:)
    real(dp), intent(out) :: residual(0:,0:)
    real(dp) :: augmentedGridFunction(0:(Nx+1),0:(Ny+1))

    residual = 0.0_dp

    call applyLHSOperatorFOCS(Nx, Ny, stepSizeX, stepSizeY, &
      gridFunction, augmentedGridFunction)

    residual(1:Nx, 1:Ny) = RHS(1:Nx, 1:Ny) - augmentedGridFunction(1:Nx, 1:Ny)

  end subroutine getResidualFOCS

  subroutine applyLHSOperatorFOCS(Nx, Ny, stepSizeX, stepSizeY, &
    gridFunction, augmentedGridFunction)
    integer, intent(in) :: Nx, Ny
    integer :: i, j
    real(dp), intent(in) :: stepSizeX, stepSizeY, gridFunction(0:,0:)
    real(dp), intent(out) :: augmentedGridFunction(0:,0:)
    real(dp) :: m1, m2, m3, m4, lambda

    lambda = stepSizeX/stepSizeY

    m1 = (1.0_dp + lambda*lambda )/2.0_dp
    m2 = - 1.0_dp + 5.0_dp*lambda*lambda
    m3 = - lambda*lambda + 5.0_dp
    m4 = - 10.0_dp*(1.0_dp + lambda**2)

    augmentedGridFunction = 0.0_dp

    do j = 1, Ny
      do i = 1, Nx

        augmentedGridFunction(i,j) = m4*gridFunction(i,j) &
          + m2*(gridFunction(i,j-1) + gridFunction(i,j+1)) &
          + m3*(gridFunction(i-1,j) + gridFunction(i+1,j)) &
          + m1*(gridFunction(i-1,j-1) + gridFunction(i+1,j-1) &
          + gridFunction(i-1,j+1) + gridFunction(i+1,j+1))

      end do
    end do

  end subroutine applyLHSOperatorFOCS

  subroutine applyRHSOperatorFOCS(Nx, Ny, stepSizeX, RHS, augmentedRHS)
    integer, intent(in) :: Nx, Ny
    integer :: i, j
    real(dp), intent(in) :: stepSizeX, RHS(0:,0:)
    real(dp), intent(out) :: augmentedRHS(0:,0:)

    augmentedRHS = 0.0_dp

    do j = 1, Ny
      do i = 1, Nx

        augmentedRHS(i,j) = (stepSizeX*stepSizeX/2.0_dp)*( &
          8.0_dp*RHS(i,j) + RHS(i,j-1) + RHS(i-1,j) &
          + RHS(i+1,j) + RHS(i,j+1))

      end do
    end do

  end subroutine applyRHSOperatorFOCS

end module FOCS2DMatvecs
