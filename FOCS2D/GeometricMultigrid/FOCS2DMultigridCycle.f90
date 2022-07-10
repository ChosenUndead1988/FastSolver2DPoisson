module FOCS2DMultigridCycle

  use FOCS2DSpecifyPrecision
  use FOCS2DConstants
  use FOCS2DDataStructures
  use FOCS2DMiscellaneous
  use FOCS2DCoarseGridSolver, only: solveCoarseSystem
  use FOCS2DMatvecs, only: getResidualFOCS
  use FOCS2DRestrictionOperator, only: RestrictionOperator
  use FOCS2DSmoother, only: SmootherGaussSeidel
  use FOCS2DInterpolationOperator, only: interpolationOperator

  implicit none

  integer, parameter :: ADDED_SMOOTHING_STEPS_COARSE = 1

  private
  public :: MGCYC, subMultigridCycleV

  contains

  subroutine MGCYC(MGStruct, flagMulticolor, cycleType, &
    orderRestriction, orderInterpolation, preSweeps, postSweeps, &
    numGrids, maximumIterations, numMGCYC, pivots, relativeTolerance, &
    absoluteTolerance, coarseGridOp, RHS, solution)
    type(MultigridStructure), intent(inout) :: MGStruct(:)
    logical, intent(in) :: flagMulticolor
    character(len=*), intent(in) :: cycleType
    integer, intent(in) :: orderRestriction, orderInterpolation, &
      preSweeps, postSweeps, numGrids, maximumIterations, pivots(:)
    integer, intent(inout) :: numMGCYC
    integer :: Nx, Ny, i, j
    real(dp), intent(inout) :: solution(0:,0:)
    real(dp), intent(in) :: coarseGridOp(:,:)
    real(dp), intent(in) :: RHS(0:,0:)
    real(dp), intent(in) :: relativeTolerance, absoluteTolerance
    real(dp) :: InfNorm, InfNorm0
    real(dp), allocatable:: residual(:,:), approx(:,:)
    logical :: stopIteration

    Nx = MGStruct(1)%gridPointsPerAxisX
    Ny = MGStruct(1)%gridPointsPerAxisY

    allocate( approx(0:(Nx + 1), 0:(Ny + 1)) )
    allocate( residual(0:(Nx + 1), 0:(Ny + 1)) )

    call getResidualFOCS(Nx, Ny, &
      MGStruct(1)%stepSizeX, &
      MGStruct(1)%stepSizeY, &
      RHS(0:,0:), solution(0:,0:), residual(0:,0:))

    InfNorm0 = maxval( abs( reshape( residual(1:Nx,1:Ny), (/ Nx*Ny /)) ) )

    numMGCYC = 0

    do

      call performMGCYC(MGStruct, flagMulticolor, cycleType, &
        orderRestriction, orderInterpolation, preSweeps, postSweeps, &
        numGrids, pivots, coarseGridOp, solution, RHS, approx)

      numMGCYC = numMGCYC + 1

      call getResidualFOCS(Nx, Ny, &
        MGStruct(1)%stepSizeX, &
        MGStruct(1)%stepSizeY, &
        RHS(0:,0:), approx(0:,0:), residual(0:,0:))

      InfNorm = maxval( abs( reshape( residual(1:Nx,1:Ny), (/ Nx*Ny /)) ) )

      solution(0:,0:) = approx(0:,0:)

      call checkTerminationSolver(InfNorm, InfNorm0, &
        relativeTolerance, absoluteTolerance, stopIteration)

      if (stopIteration .EQV. .TRUE.) then
        exit
      else if (numMGCYC .GE. maximumIterations) then
        exit
      end if

    end do

    deallocate( approx )
    deallocate( residual )

  end subroutine MGCYC

  ! This is the main subroutine to perform multigrid cycles. We may
  ! perform the V or F cycles.
  subroutine performMGCYC(MGStruct, flagMulticolor, cycleType, &
    orderRestriction, orderInterpolation, preSweeps, postSweeps, &
    numGrids, pivots, coarseGridOp, initialGuess, RHS, approximateSolution)
    logical, intent(in) :: flagMulticolor
    character(len=*), intent(in) :: cycleType
    integer, intent(in) :: orderRestriction, orderInterpolation, &
      preSweeps, postSweeps, numGrids, pivots(:)
    integer :: i
    type(MultigridStructure), intent(inout) :: MGStruct(numGrids)
    real(dp), intent(in) :: coarseGridOp(:,:)
    real(dp), intent(inout) :: initialGuess(0:,0:)
    real(dp), intent(in) :: RHS(0:,0:)
    real(dp), intent(inout) :: approximateSolution(0:,0:)

    select case (cycleType)

      case ('V','v')

        call MultigridCycleV(MGStruct, flagMulticolor, orderRestriction, &
          orderInterpolation, preSweeps, postSweeps, numGrids, pivots, &
          coarseGridOp, initialGuess, RHS, approximateSolution )

      case ('F','f')

        call MultigridCycleF(MGStruct, flagMulticolor, orderRestriction, &
          orderInterpolation, preSweeps, postSweeps, numGrids, pivots, &
          coarseGridOp, initialGuess, RHS, approximateSolution )

      case default

        call MultigridCycleV(MGStruct, flagMulticolor, orderRestriction, &
          orderInterpolation, preSweeps, postSweeps, numGrids, pivots, &
          coarseGridOp, initialGuess, RHS, approximateSolution )

    end select

  end subroutine performMGCYC

  subroutine MultigridCycleF(MGStruct, flagMulticolor, orderRestriction, &
    orderInterpolation, preSweeps, postSweeps, numGrids, pivots, &
    coarseGridOp, initialGuess, RHS, approximateSolution )
    logical, intent(in) :: flagMulticolor
    integer, intent(in) :: orderRestriction, orderInterpolation,  &
      preSweeps, postSweeps, numGrids, pivots(:)
    integer :: i
    type(MultigridStructure), intent(inout) :: MGStruct(numGrids)
    real(dp), intent(in) :: coarseGridOp(:,:)
    real(kind=dp), intent(in) :: initialGuess(0:,0:)
    real(kind=dp), intent(in) :: RHS(0:,0:)
    real(kind=dp), intent(out) :: approximateSolution(0:,0:)

    MGStruct(1)%approximateSolution(0:,0:) = initialGuess(0:,0:)
    MGStruct(1)%RHS(0:,0:) = RHS(0:,0:)

    if ( numGrids .EQ. 2 ) then

      call MultigridCycleV(MGStruct, flagMulticolor, &
        orderRestriction, orderInterpolation, &
        preSweeps, postSweeps, numGrids, pivots, coarseGridOp, &
        initialGuess, RHS, approximateSolution )

    else

      call MultigridDescentPortion(MGStruct, flagMulticolor,  &
        orderRestriction, 1, preSweeps, numGrids, pivots, coarseGridOp)

      call MultigridAscentPortion(MGStruct, flagMulticolor, &
            orderInterpolation, numgrids - 1, postSweeps, numGrids)

      do i = numGrids - 1, 2, -1

        call MultigridDescentPortion(MGStruct, flagMulticolor,  &
          orderRestriction, i, preSweeps, numGrids, pivots, coarseGridOp)

        call MultigridAscentPortion(MGStruct, flagMulticolor, &
          orderInterpolation, i - 1, postSweeps, numGrids)

      end do

      approximateSolution(0:,0:) = MGStruct(1)%approximateSolution(0:,0:)

    end if

  end subroutine MultigridCycleF

  subroutine MultigridCycleV(MGStruct, flagMulticolor, &
    orderRestriction, orderInterpolation, &
    preSweeps, postSweeps, numGrids, pivots, coarseGridOp, &
    initialGuess, RHS, approximateSolution )
    logical, intent(in) :: flagMulticolor
    integer, intent(in) :: orderRestriction, orderInterpolation,  &
      preSweeps, postSweeps, numGrids, pivots(:)
    type(MultigridStructure), intent(inout) :: MGStruct(numGrids)
    real(dp), intent(in) :: coarseGridOp(:,:)
    real(dp), intent(in) :: initialGuess(0:,0:)
    real(dp), intent(in) :: RHS(0:,0:)
    real(dp), intent(out) :: approximateSolution(0:,0:)

    MGStruct(1)%approximateSolution(0:,0:) = initialGuess(0:,0:)
    MGStruct(1)%RHS(0:,0:) = RHS(0:,0:)

    call MultigridDescentPortion(MGStruct, flagMulticolor,  &
      orderRestriction, 1, preSweeps, numGrids, pivots, coarseGridOp)

    call MultigridAscentPortion(MGStruct, flagMulticolor, &
      orderInterpolation, 1, postSweeps, numGrids)

    approximateSolution(0:,0:) = MGStruct(1)%approximateSolution(0:,0:)

  end subroutine MultigridCycleV

  subroutine subMultigridCycleV(MGStruct, flagMulticolor, &
    orderRestriction, orderInterpolation, &
    preSweeps, postSweeps, topLevel, numGrids, pivots, coarseGridOp, &
    initialGuessTopLevel, RHSTopLevel, approximateSolutionTopLevel )
    logical, intent(in) :: flagMulticolor
    integer, intent(in) :: orderRestriction, orderInterpolation,  &
      preSweeps, postSweeps, topLevel, numGrids, pivots(:)
    type(MultigridStructure), intent(inout) :: MGStruct(numGrids)
    real(dp), intent(in) :: coarseGridOp(:,:)
    real(dp), intent(in) :: initialGuessTopLevel(0:,0:)
    real(dp), intent(in) :: RHSTopLevel(0:,0:)
    real(dp), intent(out) :: approximateSolutionTopLevel(0:,0:)

    MGStruct(topLevel)%approximateSolution(0:,0:) = initialGuessTopLevel(0:,0:)
    MGStruct(topLevel)%RHS(0:,0:) = RHSTopLevel(0:,0:)

    call MultigridDescentPortion(MGStruct, flagMulticolor,  &
      orderRestriction, topLevel, preSweeps, numGrids, pivots, coarseGridOp)

    call MultigridAscentPortion(MGStruct, flagMulticolor, &
      orderInterpolation, topLevel, postSweeps, numGrids)

    approximateSolutionTopLevel(0:,0:) = MGStruct(topLevel)%approximateSolution(0:,0:)

  end subroutine subMultigridCycleV

  ! Preform the descending portion of a given multigrid cycle. Here topLevel = 1
  ! on the finest grid and topLevel = numGrids on the coarsest grid
  subroutine MultigridDescentPortion(MGStruct, flagMulticolor,  &
    orderRestriction, topLevel, numSweeps, numGrids, pivots, coarseGridOp)
    logical, intent(in) :: flagMulticolor
    integer, intent(in) :: orderRestriction, topLevel, numGrids, &
      numSweeps, pivots(:)
    integer :: i
    type(MultigridStructure), intent(inout) :: MGStruct(numGrids)
    real(dp), intent(in) :: coarseGridOp(:,:)
    real(kind=dp), allocatable :: res(:,:)

    do i = topLevel, numGrids - 1

      ! allocate space for the residaul on the fine grid
      allocate( res(0:(MGStruct(i)%gridPointsPerAxisX + 1), &
        0:(MGStruct(i)%gridPointsPerAxisY + 1)) )

      ! apply the smoother to fine grids system
      call SmootherGaussSeidel(flagMulticolor, &
        MGStruct(i)%gridPointsPerAxisX, &
        MGStruct(i)%gridPointsPerAxisY, &
        numSweeps, &
        MGStruct(i)%stepSizeX, &
        MGStruct(i)%stepSizeY, &
        MGStruct(i)%RHS(0:,0:), &
        MGStruct(i)%approximateSolution(0:,0:))

      !compute the residual on the fine grid
      call getResidualFOCS(MGStruct(i)%gridPointsPerAxisX, &
        MGStruct(i)%gridPointsPerAxisY, &
        MGStruct(i)%stepSizeX, &
        MGStruct(i)%stepSizeY, &
        MGStruct(i)%RHS(0:,0:), &
        MGStruct(i)%approximateSolution(0:,0:), &
        res(0:,0:))

      ! use restriction operator to form the RHS vector of the coarse grid
      ! with the residual of the fine grid
      call RestrictionOperator(orderRestriction, &
        MGStruct(i + 1)%gridPointsPerAxisX, &
        MGStruct(i + 1)%gridPointsPerAxisY, &
        res(0:,0:), &
        MGStruct(i + 1)%RHS(0:,0:))

      deallocate( res )

    end do

    ! Solve the coarsest grid problem  directly
    call solveCoarseSystem(MGStruct(numgrids)%gridPointsPerAxisX, &
      MGStruct(numGrids)%gridPointsPerAxisY, pivots, coarseGridOp, &
      MGStruct(numGrids)%RHS(0:,0:), &
      MGStruct(numGrids)%approximateSolution(0:,0:))

  end subroutine MultigridDescentPortion

  subroutine MultigridAscentPortion(MGStruct, flagMulticolor, &
    orderInterpolation, topLevel, numSweeps, numGrids)
    logical, intent(in) :: flagMulticolor
    integer, intent(in) :: orderInterpolation, topLevel, numSweeps, numGrids
    integer :: i
    type(MultigridStructure), intent(inout) :: MGStruct(numGrids)
    real(kind=dp), allocatable :: interpolationFine(:,:)

    do i = numGrids - 1, topLevel, - 1

      ! allocate space for the correction term on the fine grid
      allocate( interpolationFine(0:(MGStruct(i)%gridPointsPerAxisX + 1), &
        0:(MGStruct(i)%gridPointsPerAxisY + 1)) )

      interpolationFine = 0.0_dp

      ! Use the interpolation operator the transfer from the coarse grid to the
      ! fine grid.
      call interpolationOperator(orderInterpolation, &
        MGStruct(i + 1)%gridPointsPerAxisX, &
        MGStruct(i + 1)%gridPointsPerAxisY, &
        MGStruct(i + 1)%approximateSolution, &
        interpolationFine)

      ! Take the grid function prolongated from the coarse grid and
      ! correct the approximate solution on the fine grid.

      MGStruct(i)%approximateSolution(0:,0:) = &
        MGStruct(i)%approximateSolution(0:,0:) + interpolationFine(0:,0:)

      ! Apply the smoother to further reduce the error
      if ( i .EQ. numGrids - 1) then

        call SmootherGaussSeidel(flagMulticolor, &
          MGStruct(i)%gridPointsPerAxisX, &
          MGStruct(i)%gridPointsPerAxisY, &
          numSweeps + ADDED_SMOOTHING_STEPS_COARSE, MGStruct(i)%stepSizeX, MGStruct(i)%stepSizeY, &
          MGStruct(i)%RHS(0:,0:), MGStruct(i)%approximateSolution(0:,0:))

      else

        call SmootherGaussSeidel(flagMulticolor, &
          MGStruct(i)%gridPointsPerAxisX, &
          MGStruct(i)%gridPointsPerAxisY, &
          numSweeps, MGStruct(i)%stepSizeX, MGStruct(i)%stepSizeY, &
          MGStruct(i)%RHS(0:,0:), MGStruct(i)%approximateSolution(0:,0:))

      end if

      deallocate( interpolationFine )

      ! make the grid function zero on the coarse grid
      MGStruct(i+1)%approximateSolution = 0.0_dp

    end do

  end subroutine MultigridAscentPortion

end module FOCS2DMultigridCycle
