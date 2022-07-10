module TestProblem1

  use FOCS2DSpecifyPrecision
  use FOCS2DConstants

  implicit none

  real(dp), parameter :: LEFT_ENDPOINT_X_TEST1 = 0.0_dp
  real(dp), parameter :: RIGHT_ENDPOINT_X_TEST1 = 2.0_dp
  real(dp), parameter :: LEFT_ENDPOINT_Y_TEST1 = 0.0_dp
  real(dp), parameter :: RIGHT_ENDPOINT_Y_TEST1 = 1.0_dp
  real(dp), parameter :: parameter_X = 10.0_dp*pi
  real(dp), parameter :: parameter_Y = 2.0_dp*pi

  contains

  subroutine getAnalyticSolutionTestProblem1(x, y, output)
    real(dp), intent(in) :: x, y
    real(dp), intent(out) :: output
    real(dp) :: a, b

    a = parameter_X
    b = parameter_Y

    output = sin(a*x)*sin(b*y)

  end subroutine getAnalyticSolutionTestProblem1

  subroutine getRHSFunctionTestProblem1(x, y, output)
    real(dp), intent(in) :: x, y
    real(dp), intent(out) :: output
    real(dp) :: a, b

    a = parameter_X
    b = parameter_Y

    output = - (a*a + b*b)*sin(a*x)*sin(b*y)

  end subroutine getRHSFunctionTestProblem1

  subroutine getDirichletBoundaryConditionNorthTestProblem1(x, output)
    real(dp), intent(in) :: x
    real(dp), intent(out) :: output
    real(dp) :: a, b

    a = parameter_X
    b = parameter_Y

    output = sin(a*x)*sin(b*RIGHT_ENDPOINT_Y_TEST1)

  end subroutine getDirichletBoundaryConditionNorthTestProblem1

  subroutine getDirichletBoundaryConditionEastTestProblem1(y, output)
    real(dp), intent(in) :: y
    real(dp), intent(out) :: output
    real(dp) :: a, b

    a = parameter_X
    b = parameter_Y

    output = sin(a*RIGHT_ENDPOINT_X_TEST1)*sin(b*y)

  end subroutine getDirichletBoundaryConditionEastTestProblem1

  subroutine getDirichletBoundaryConditionSouthTestProblem1(x, output)
    real(dp), intent(in) :: x
    real(dp), intent(out) :: output
    real(dp) :: a, b

    a = parameter_X
    b = parameter_Y

    output = sin(a*x)*sin(b*LEFT_ENDPOINT_Y_TEST1)

  end subroutine getDirichletBoundaryConditionSouthTestProblem1

  subroutine getDirichletBoundaryConditionWestTestProblem1(y, output)
    real(dp), intent(in) :: y
    real(dp), intent(out) :: output
    real(dp) :: a, b

    a = parameter_X
    b = parameter_Y

    output = sin(a*LEFT_ENDPOINT_X_TEST1)*sin(b*y)

  end subroutine getDirichletBoundaryConditionWestTestProblem1

end module TestProblem1
