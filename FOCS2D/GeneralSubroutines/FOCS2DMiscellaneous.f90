module FOCS2DMiscellaneous

  use FOCS2DSpecifyPrecision
  use FOCS2DConstants

  implicit none

  contains

  function KroneckerDelta(i,j) result(delta)
    integer, intent(in) :: i,j
    real(dp) :: delta

    if (i .eq. j) then
      delta = 1.0_dp
    else
      delta = 0.0_dp
    end if

  end function

  subroutine getGridPointPerAxis(leadingFactor, numGrids, N)
    integer, intent(in) :: leadingFactor, numGrids
    integer, intent(out) :: N

    N = leadingFactor*(2**numGrids) - 1

  end subroutine getGridPointPerAxis

  subroutine getCoordinates(i, j, gridPointsPerAxisX, gridPointsPerAxisY, &
    leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY, x, y)
    integer, intent(in) :: i, j, gridPointsPerAxisX, gridPointsPerAxisY
    real(dp), intent(in) :: leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY
    real(dp), intent(out) :: x, y
    real(dp) :: stepSizeX, stepSizeY

    call getStepSize(gridPointsPerAxisX, gridPointsPerAxisY, &
      leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY, &
      stepSizeX, stepSizeY)

    x = leftEndpointX + real(i, dp)*stepSizeX
    y = leftEndpointY + real(j, dp)*stepSizeY

  end subroutine

  subroutine getStepSize(gridPointsPerAxisX, gridPointsPerAxisY, &
    leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY, &
    stepSizeX, stepSizeY)
    integer, intent(in) :: gridPointsPerAxisX, gridPointsPerAxisY
    real(dp), intent(in) :: leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY
    real(dp), intent(out) :: stepSizeX, stepSizeY

    stepSizeX = (rightEndpointX - leftEndpointX)/real(gridPointsPerAxisX + 1, dp)
    stepSizeY = (rightEndpointY - leftEndpointY)/real(gridPointsPerAxisY + 1, dp)

  end subroutine getStepSize

  subroutine getInfinityNormError(analyticSolution, &
    gridPointsPerAxisX, gridPointsPerAxisY, &
    leftEndpointX, rightEndpointX, &
    leftEndpointY, rightEndpointY, &
    approximateSolution, error)
    interface
      subroutine analyticSolution(x, y, solution)
        use FOCS2DSpecifyPrecision
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: solution
      end subroutine analyticSolution
    end interface

    integer, intent(in) :: gridPointsPerAxisX, gridPointsPerAxisY
    integer :: i, j
    real(dp), intent(in) :: leftEndpointX, rightEndpointX, &
      leftEndpointY, rightEndpointY, approximateSolution(0:,0:)
    real(dp) :: error
    real(dp) :: deviation, norm, true, x, y

    norm = 0.0_dp

    do j = 1, gridPointsPerAxisY
      do i = 1, gridPointsPerAxisX

        call getCoordinates(i, j, gridPointsPerAxisX, gridPointsPerAxisY, &
          leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY, x, y)

        call analyticSolution(x, y, true)

        deviation = abs( approximateSolution(i,j) - true )

        if (deviation .gt. norm) then
          norm = deviation
        end if

      end do
    end do

    error = norm

  end subroutine getInfinityNormError

  subroutine checkTerminationSolver(normCurrent, normInitial, &
    relativeToleranceResidual, absoluteToleranceResidual, stopIteration)
    real(dp), intent(in) :: relativeToleranceResidual, &
      absoluteToleranceResidual, normCurrent, normInitial
    real(kind=dp) :: bound
    logical, intent(out) :: stopIteration

    bound = relativeToleranceResidual*normInitial + absoluteToleranceResidual

    if (normCurrent .LE. bound) then
      stopIteration = .TRUE.
    else
      stopIteration = .FALSE.
    end if

  end subroutine checkTerminationSolver

end module FOCS2DMiscellaneous
