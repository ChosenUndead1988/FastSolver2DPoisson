module FOCS2DFastPoissonSolver

  use FOCS2DSpecifyPrecision
  use FOCS2DConstants
  use FOCS2DMiscellaneous
  use FOCS2DSettings
  use FOCS2DDataStructures
  use FOCS2DCoarseGridSolver
  use FOCS2DMatvecs
  use FOCS2DRestrictionOperator
  use FOCS2DSmoother
  use FOCS2DInterpolationOperator
  use FOCS2DMultigridCycle

  implicit none

  private
  public :: FastPoissonSolverHomogeneousDirichletFOCS

  contains

  subroutine FastPoissonSolverHomogeneousDirichletFOCS(analyticSolution, &
    evaluateRHS, leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY)
    interface
      subroutine evaluateRHS(x, y, RHS)
        use FOCS2DSpecifyPrecision
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: RHS
      end subroutine evaluateRHS

      subroutine analyticSolution(x, y, RHS)
        use FOCS2DSpecifyPrecision
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: RHS
      end subroutine analyticSolution

    end interface
    integer, allocatable :: pivots(:)
    integer :: Nc, numMGCYC
    real(dp), intent(in) :: leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY
    real(dp), allocatable :: RHS_FOCS(:,:), solution(:,:), coarseGridOp(:,:)
    real(dp) :: error
    type(MultigridStructure) :: MGStruct(NUMBER_AUXILIARY_GRIDS)

    call initiateMGStructure(MGStruct, LEADING_FACTOR, &
      LEADING_FACTOR, NUMBER_AUXILIARY_GRIDS, &
      leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY)

    Nc = MGStruct(NUMBER_AUXILIARY_GRIDS)%gridPointsPerAxisX*MGStruct(NUMBER_AUXILIARY_GRIDS)%gridPointsPerAxisY

    allocate( coarseGridOp(Nc,Nc))
    allocate( pivots(Nc) )
    allocate( RHS_FOCS(0:(MGStruct(1)%gridPointsPerAxisX + 1), &
      0:(MGStruct(1)%gridPointsPerAxisY + 1)) )
    allocate( solution(0:(MGStruct(1)%gridPointsPerAxisX + 1), &
      0:(MGStruct(1)%gridPointsPerAxisY + 1)) )

    RHS_FOCS = 0.0_dp
    solution = 0.0_dp
    coarseGridOp = 0.0_dp

    call getCoarseGridOperator(Nc, Nc, &
      MGStruct(NUMBER_AUXILIARY_GRIDS)%stepSizeX, &
      MGStruct(NUMBER_AUXILIARY_GRIDS)%stepSizeY, coarseGridOp)

    call getLUFactorizationCoarseGridOperator(Nc, pivots, coarseGridOp )

    call NestedIterationHomogeneousDirichlet(evaluateRHS, &
      MGStruct, MULTICOLOR_SMOOTHER, &
      leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY, &
      ORDER_RESTRICTION_OPERATOR, ORDER_INTERPOLATION_OPERATOR, &
      NUM_SWEEPS_DESCENT_MULTIGRID, &
      NUM_SWEEPS_ASCENT_MULTIGRID, &
      NUMBER_AUXILIARY_GRIDS, NUM_CYCLES_NESTED_ITERATION, &
      pivots, coarseGridOp, solution)

    call getRHSDirichletBoundaryCondition(evaluateRHS, &
      MGStruct(1)%gridPointsPerAxisX, &
      MGStruct(1)%gridPointsPerAxisY, &
      leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY, &
      RHS_FOCS)

    call MGCYC(MGStruct, MULTICOLOR_SMOOTHER, CYCLE_TYPE_MULTIGRID, &
      ORDER_RESTRICTION_OPERATOR, ORDER_INTERPOLATION_OPERATOR, &
      NUM_SWEEPS_DESCENT_MULTIGRID, &
      NUM_SWEEPS_ASCENT_MULTIGRID, &
      NUMBER_AUXILIARY_GRIDS, &
      MAXIMUM_CYCLES_MULTIGRID, &
      numMGCYC, pivots, RELATIVE_TOLERANCE_MULTIGRID, &
      ABSOLUTE_TOLERANCE_MULTIGRID, coarseGridOp, RHS_FOCS, solution)

    call getInfinityNormError(analyticSolution, &
      MGStruct(1)%gridPointsPerAxisX, MGStruct(1)%gridPointsPerAxisY, &
      leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY, &
      Solution, error)

    print *, 'Infinity Norm:', error

    call deconstructMGStructure(MGStruct, NUMBER_AUXILIARY_GRIDS)

  end subroutine FastPoissonSolverHomogeneousDirichletFOCS

  subroutine initiateMGStructure(MGStruct, leadingFactorX, leadingFactorY, &
    numGrids, leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY)
    integer, intent(in) :: leadingFactorX, leadingFactorY, numGrids
    integer :: loop, count
    real(dp), intent(in) :: leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY
    type(MultigridStructure), intent(inout) :: MGStruct(numGrids)

    do loop = numGrids, 1, -1

      count = numGrids + 1 - loop

      call getGridPointPerAxis(leadingFactorX, loop, MGStruct(count)%gridPointsPerAxisX)
      call getGridPointPerAxis(leadingFactorY, loop, MGStruct(count)%gridPointsPerAxisY)

      call getStepSize(MGStruct(count)%gridPointsPerAxisX, &
        MGStruct(count)%gridPointsPerAxisY, &
        leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY, &
        MGStruct(count)%stepSizeX, MGStruct(count)%stepSizeY)

      allocate( MGStruct(count)%RHS(0:(MGStruct(count)%gridPointsPerAxisX+1), &
        0:(MGStruct(count)%gridPointsPerAxisY+1)) )
      allocate( MGStruct(count)%approximateSolution(0:(MGStruct(count)%gridPointsPerAxisX+1), &
        0:(MGStruct(count)%gridPointsPerAxisY+1)) )

      MGStruct(count)%RHS = 0.0_dp
      MGStruct(count)%approximateSolution = 0.0_dp

    end do

  end subroutine initiateMGStructure

  subroutine deconstructMGStructure(MGStruct, numGrids)
    integer, intent(in) ::  numGrids
    integer :: loop
    type(MultigridStructure), intent(inout) :: MGStruct(numGrids)

    do loop = numGrids, 1, -1

      deallocate( MGStruct(loop)%RHS )
      deallocate( MGStruct(loop)%approximateSolution )

    end do

  end subroutine deconstructMGStructure

  subroutine NestedIterationHomogeneousDirichlet(evaluateRHS, &
    MGStruct, flagMulticolor, &
    leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY, &
    orderRestriction, orderInterpolation, preSweeps, postSweeps, &
    numGrids, numCyclesNI, pivots, coarseGridOp, initialGuess)
    interface
      subroutine evaluateRHS(x, y, RHS)
        use FOCS2DSpecifyPrecision
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: RHS
      end subroutine evaluateRHS

    end interface

    logical, intent(in) :: flagMulticolor
    integer, intent(in) :: orderRestriction, orderInterpolation,  &
      preSweeps, postSweeps, numGrids, numCyclesNI, pivots(:)
    integer :: topLevel, i, j, Nx, Ny
    real(dp), intent(in) :: leftEndpointX, rightEndpointX, leftEndpointY, &
      rightEndpointY, coarseGridOp(:,:)
    real(dp), intent(out) :: initialGuess(0:,0:)
    real(dp), allocatable :: approx(:,:)
    real(dp) :: x, y
    type(MultigridStructure), intent(inout) :: MGStruct(numGrids)

    do topLevel = numGrids, 2, -1
      if (topLevel .eq. numGrids ) then

        call getRHSDirichletBoundaryCondition(evaluateRHS, &
          MGStruct(numGrids)%gridPointsPerAxisX, &
          MGStruct(numGrids)%gridPointsPerAxisY, &
          leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY, &
          MGStruct(numGrids)%RHS(0:,0:))

        call solveCoarseSystem(MGStruct(numGrids)%gridPointsPerAxisX, &
          MGStruct(numGrids)%gridPointsPerAxisY, pivots, coarseGridOp, &
          MGStruct(numGrids)%RHS(0:,0:), &
          MGStruct(numGrids)%approximateSolution(0:,0:))

        call interpolationOperator(ORDER_INTERPOLATION_NESTED_ITERATION, &
          MGStruct(numGrids)%gridPointsPerAxisX, &
          MGStruct(numGrids)%gridPointsPerAxisY, &
          MGStruct(numGrids)%approximateSolution(0:,0:), &
          MGStruct(numGrids - 1)%approximateSolution(0:,0:))

        MGStruct(numGrids)%RHS(0:,0:) = 0.0_dp
        MGStruct(numGrids)%approximateSolution(0:,0:) = 0.0_dp

      else

        call getRHSDirichletBoundaryCondition(evaluateRHS, &
          MGStruct(topLevel)%gridPointsPerAxisX, &
          MGStruct(topLevel)%gridPointsPerAxisY, &
          leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY, &
          MGStruct(topLevel)%RHS(0:,0:))

        allocate( approx(0:(MGStruct(topLevel)%gridPointsPerAxisX + 1), &
          0:(MGStruct(topLevel)%gridPointsPerAxisY + 1))  )

        approx = 0.0_dp

        do i = 1, numCyclesNI

          call subMultigridCycleV(MGStruct, flagMulticolor, &
            orderRestriction, orderInterpolation, &
            preSweeps, postSweeps, topLevel, numGrids, pivots, coarseGridOp, &
            MGStruct(topLevel)%approximateSolution(0:,0:), &
            MGStruct(topLevel)%RHS(0:,0:), approx )

          MGStruct(topLevel)%approximateSolution(0:,0:) = approx(0:,0:)

        end do

        deallocate( approx )

        call interpolationOperator(ORDER_INTERPOLATION_NESTED_ITERATION, &
          MGStruct(topLevel)%gridPointsPerAxisX, &
          MGStruct(topLevel)%gridPointsPerAxisY, &
          MGStruct(topLevel)%approximateSolution(0:,0:), &
          MGStruct(topLevel - 1)%approximateSolution(0:,0:))

        MGStruct(topLevel)%approximateSolution(0:,0:) = 0.0_dp
        MGStruct(topLevel)%RHS(0:,0:) = 0.0_dp

        if (topLevel .eq. 2) then
          initialGuess(0:,0:) = MGStruct(1)%approximateSolution(0:,0:)
        end if

      end if
    end do

  end subroutine NestedIterationHomogeneousDirichlet

  subroutine getRHSDirichletBoundaryCondition(evaluateRHS, &
    gridPointsPerAxisX, gridPointsPerAxisY, &
    leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY, &
    RHS_FOCS)
    interface

      subroutine evaluateRHS(x, y, RHS)
        use FOCS2DSpecifyPrecision
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: RHS
      end subroutine evaluateRHS

    end interface

    integer, intent(in) :: gridPointsPerAxisX, gridPointsPerAxisY
    integer :: i, j
    real(dp), intent(in) :: leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY
    real(dp), intent(out) :: RHS_FOCS(0:,0:)
    real(dp), allocatable :: RHS(:,:)
    real(dp) :: x, y, stepSizeX, stepSizeY

    allocate( RHS(0:(gridPointsPerAxisX + 1), &
      0:(gridPointsPerAxisY + 1)))

    do j = 0, gridPointsPerAxisY + 1
      do i = 0, gridPointsPerAxisX + 1

        call getCoordinates(i, j, gridPointsPerAxisX, gridPointsPerAxisY, &
          leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY, x, y)

        call evaluateRHS(x, y, RHS(i,j))

      end do
    end do

    call getStepSize(gridPointsPerAxisX, gridPointsPerAxisY, &
      leftEndpointX, rightEndpointX, leftEndpointY, rightEndpointY, &
      stepSizeX, stepSizeY)

    call applyRHSOperatorFOCS(gridPointsPerAxisX, gridPointsPerAxisY, &
      stepSizeX, RHS, RHS_FOCS)

    deallocate( RHS )

  end subroutine getRHSDirichletBoundaryCondition

end module FOCS2DFastPoissonSolver
