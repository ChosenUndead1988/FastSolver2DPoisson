program MainProgramTest1FOCS

  use FOCS2DSpecifyPrecision
  use FOCS2DConstants
  use TestProblem1
  use FOCS2DDataStructures
  use FOCS2DSettings
  use FOCS2DCoarseGridSolver
  use FOCS2DMatvecs
  use FOCS2DRestrictionOperator
  use FOCS2DSmoother
  use FOCS2DInterpolationOperator
  use FOCS2DMultigridCycle
  use FOCS2DFastPoissonSolver

  call FastPoissonSolverHomogeneousDirichletFOCS(&
    getAnalyticSolutionTestProblem1, &
    getRHSFunctionTestProblem1, &
    LEFT_ENDPOINT_X_TEST1, &
    RIGHT_ENDPOINT_X_TEST1, &
    LEFT_ENDPOINT_Y_TEST1, &
    RIGHT_ENDPOINT_Y_TEST1)

end program MainProgramTest1FOCS
