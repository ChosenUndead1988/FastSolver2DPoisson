module FOCS2DSmoother

  use FOCS2DSpecifyPrecision
  use FOCS2DConstants

  implicit none

  private
  public :: SmootherGaussSeidel

  contains

  subroutine SmootherGaussSeidel(flagMulticolor, Nx, Ny, numSweeps, &
    stepSizeX, stepSizeY, RHS, gridFunction)
    logical, intent(in) :: flagMulticolor
    integer, intent(in) :: Nx, Ny, numSweeps
    integer :: loop
    real(dp), intent(in) :: stepSizeX, stepSizeY, RHS(0:,0:)
    real(dp), intent(inout) :: gridFunction(0:,0:)

    if (flagMulticolor .eqv. .true.) then
      do loop = 1, numSweeps
        call RedBlackGaussSeidel(Nx, Ny, stepSizeX, stepSizeY, RHS, gridFunction)
      end do
    else
      do loop = 1, numSweeps
        call lexicographicGaussSeidel(Nx, Ny, stepSizeX, stepSizeY, RHS, gridFunction)
      end do
    end if

  end subroutine SmootherGaussSeidel

  subroutine RedBlackGaussSeidel(Nx, Ny, stepSizeX, stepSizeY, RHS, gridFunction)
    integer, intent(in) :: Nx, Ny
    integer :: i, j
    real(dp), intent(in) :: stepSizeX, stepSizeY, RHS(0:,0:)
    real(dp), intent(inout) :: gridFunction(0:,0:)

    ! red nodes - (odd, odd)
    do j = 1, Ny, 2
      do i = 1, Nx, 2

        call localGaussSeidel(i, j, stepSizeX, stepSizeY, RHS, gridFunction)

      end do
    end do

    ! blue nodes - (even, even)
    do j = 2, Ny - 1, 2
      do i = 2, Nx - 1, 2

        call localGaussSeidel(i, j, stepSizeX, stepSizeY, RHS, gridFunction)

      end do
    end do

    ! green nodes - (even, odd)
    do j = 1, Ny, 2
      do i = 2, Nx - 1, 2

        call localGaussSeidel(i, j, stepSizeX, stepSizeY, RHS, gridFunction)

      end do
    end do

    ! black nodes - (odd, even)
    do j = 2, Ny - 1, 2
      do i = 1, Nx, 2

        call localGaussSeidel(i, j, stepSizeX, stepSizeY, RHS, gridFunction)

      end do
    end do

  end subroutine RedBlackGaussSeidel

  subroutine lexicographicGaussSeidel(Nx, Ny, stepSizeX, stepSizeY, RHS, gridFunction)
    integer, intent(in) :: Nx, Ny
    integer :: i, j
    real(dp), intent(in) :: stepSizeX, stepSizeY, RHS(0:,0:)
    real(dp), intent(inout) :: gridFunction(0:,0:)

    do j = 1, Ny
      do i = 1, Nx

        call localGaussSeidel(i, j, stepSizeX, stepSizeY, RHS, gridFunction)

      end do
    end do

  end subroutine lexicographicGaussSeidel

  subroutine localGaussSeidel(i, j, stepSizeX, stepSizeY, RHS, gridFunction)
    integer, intent(in) :: i, j
    real(dp), intent(in) :: stepSizeX, stepSizeY, RHS(0:,0:)
    real(dp), intent(inout) :: gridFunction(0:,0:)
    real(dp) :: m1, m2, m3, m4, lambda, corner, sideX, sideY

    lambda = stepSizeX/stepSizeY

    m1 = (1.0_dp + lambda*lambda )/2.0_dp
    m2 = - 1.0_dp + 5.0_dp*lambda*lambda
    m3 = - lambda*lambda + 5.0_dp
    m4 = - 10.0_dp*(1.0_dp + lambda**2)

    sideY = m2*(gridFunction(i,j+1) + gridFunction(i,j-1))
    sideX = m3*(gridFunction(i+1,j) + gridFunction(i-1,j))

    corner = m1*(gridFunction(i-1,j-1) + gridFunction(i+1,j-1) &
      + gridFunction(i-1,j+1) + gridFunction(i+1,j+1))

    gridFunction(i,j) = (RHS(i,j) - sideX - sideY - corner)/m4

  end subroutine localGaussSeidel

end module FOCS2DSmoother
