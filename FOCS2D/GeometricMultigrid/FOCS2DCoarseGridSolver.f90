module FOCS2DCoarseGridSolver

  use FOCS2DSpecifyPrecision
  use FOCS2DConstants

  implicit none

  private
  public :: getCoarseGridOperator, getLUFactorizationCoarseGridOperator, solveCoarseSystem


  contains

  subroutine solveCoarseSystem(gridPointsPerAxisX, gridPointsPerAxisY, pivots, coarseGridOp, RHS, solution)
    integer, intent(in) :: gridPointsPerAxisX, gridPointsPerAxisY, pivots(:)
    integer :: i, j, loop, numVars, info
    real(dp), intent(in) :: coarseGridOp(:,:)
    real(dp), intent(in) :: RHS(0:,0:)
    real(dp), intent(out) :: solution(0:,0:)
    real(dp), allocatable :: x(:)

    numVars = gridPointsPerAxisX*gridPointsPerAxisY
    allocate( x(numVars) )

    do j = 1, gridPointsPerAxisY
      do i = 1, gridPointsPerAxisX

        call grid2Vec(i, j, gridPointsPerAxisX, gridPointsPerAxisY, loop)
        x(loop) = RHS(i,j)

      end do
    end do

    call DGETRS( 'N', numVars, 1, coarseGridOp, numVars, pivots, &
      x, numVars, info )

    solution = 0.0_dp

    do j = 1, gridPointsPerAxisY
      do i = 1, gridPointsPerAxisX

        call grid2Vec(i, j, gridPointsPerAxisX, gridPointsPerAxisY, loop)
        solution(i,j) = x(loop)

      end do
    end do


  end subroutine solveCoarseSystem

  subroutine getLUFactorizationCoarseGridOperator(numVarsCoarse, pivots, &
    coarseGridOp )
    integer, intent(in) :: numVarsCoarse
    integer, intent(out) :: pivots(:)
    integer :: info
    real(dp), intent(inout) :: coarseGridOp(:,:)

    call DGETRF( numVarsCoarse, numVarsCoarse, coarseGridOp, &
      numVarsCoarse, pivots, info )

  end subroutine getLUFactorizationCoarseGridOperator

  subroutine getCoarseGridOperator(gridPointsPerAxisX, gridPointsPerAxisY, &
    stepSizeX, stepSizeY, coarseGridOp)
    integer, intent(in) :: gridPointsPerAxisX, gridPointsPerAxisY
    integer :: i, j, row, col
    real(dp), intent(in) :: stepSizeX, stepSizeY
    real(dp), intent(out) :: coarseGridOp(:,:)
    real(dp) :: lambda, m1, m2, m3, m4

    lambda = stepSizeX/stepSizeY

    m1 = (1.0_dp + lambda*lambda )/2.0_dp
    m2 = - 1.0_dp + 5.0_dp*lambda*lambda
    m3 = - lambda*lambda + 5.0_dp
    m4 = - 10.0_dp*(1.0_dp + lambda**2)

    coarseGridOp = 0.0_dp

    row = 0

    do j = 1, gridPointsPerAxisY
      do i = 1, gridPointsPerAxisX

        row = row + 1

        call grid2Vec(i, j, gridPointsPerAxisX, gridPointsPerAxisY, col)

        coarseGridOp(row, col) = m4

        call grid2Vec(i, j - 1, gridPointsPerAxisX, gridPointsPerAxisY, col)

        if (col .gt. 0) then
          coarseGridOp(row, col) = m2
        end if

        call grid2Vec(i, j + 1, gridPointsPerAxisX, gridPointsPerAxisY, col)

        if (col .gt. 0) then
          coarseGridOp(row, col) = m2
        end if

        call grid2Vec(i - 1, j, gridPointsPerAxisX, gridPointsPerAxisY, col)

        if (col .gt. 0) then
          coarseGridOp(row, col) = m3
        end if

        call grid2Vec(i + 1, j, gridPointsPerAxisX, gridPointsPerAxisY, col)

        if (col .gt. 0) then
          coarseGridOp(row, col) = m3
        end if

        call grid2Vec(i - 1, j - 1, gridPointsPerAxisX, gridPointsPerAxisY, col)

        if (col .gt. 0) then
          coarseGridOp(row, col) = m1
        end if

        call grid2Vec(i + 1, j - 1, gridPointsPerAxisX, gridPointsPerAxisY, col)

        if (col .gt. 0) then
          coarseGridOp(row, col) = m1
        end if

        call grid2Vec(i - 1, j + 1, gridPointsPerAxisX, gridPointsPerAxisY, col)

        if (col .gt. 0) then
          coarseGridOp(row, col) = m1
        end if

        call grid2Vec(i + 1, j + 1, gridPointsPerAxisX, gridPointsPerAxisY, col)

        if (col .gt. 0) then
          coarseGridOp(row, col) = m1
        end if

      end do
    end do

  end subroutine getCoarseGridOperator

  subroutine grid2Vec(i, j, gridPointsPerAxisX, gridPointsPerAxisY, index)
    integer, intent(in) :: i, j, gridPointsPerAxisX, gridPointsPerAxisY
    integer, intent(out) :: index

    if ((i .lt. 1) .or. (i .gt. gridPointsPerAxisX) &
      .or. (j .lt. 1) .or. (j .gt. gridPointsPerAxisY)) then
      index = - 1
    else
      index = (j - 1)*gridPointsPerAxisX + i
    end if

  end subroutine grid2Vec

end module FOCS2DCoarseGridSolver
