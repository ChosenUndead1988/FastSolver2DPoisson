module FOCS2DInterpolationOperator

  use FOCS2DSpecifyPrecision
  use FOCS2DConstants

  ! precomputed coefficients for the real valued linear interpolation
  real(kind=dp), parameter, dimension(2) :: LINEAR_LAGRANGE = (/ 0.5_dp, 0.5_dp /)
  ! precomputed coefficients for the real valued quadratic interpolation
  real(kind=dp), parameter, dimension(3,2) :: QUADRATIC_LAGRANGE = &
  reshape( (/ 3.0_dp/8.0_dp, & ! (1,1)
          3.0_dp/4.0_dp, & ! (2,1)
        - 1.0_dp/8.0_dp, & ! (3,1)
        - 1.0_dp/8.0_dp, & ! (1,2)
          3.0_dp/4.0_dp, & ! (2,2)
          3.0_dp/8.0_dp /),&
          (/3, 2/) )
  ! precomputed coefficients for the real valued cubic interpolation
  real(kind=dp), parameter, dimension(4,3) :: CUBIC_LAGRANGE = &
    reshape( &
    (/  5.0_dp/16.0_dp, &    ! (1,1)
      15.0_dp/16.0_dp,  &    ! (2,1)
      - 5.0_dp/16.0_dp,   &    ! (3,1)
      1.0_dp/16.0_dp,   &    ! (4,1)
      - 1.0_dp/16.0_dp,   &    ! (1,2)
      9.0_dp/16.0_dp,   &    ! (2,2)
      9.0_dp/16.0_dp,   &    ! (3,2)
      - 1.0_dp/16.0_dp,   &    ! (4,2)
      1.0_dp/16.0_dp,   &    ! (1,3)
      - 5.0_dp/16.0_dp,   &    ! (2,3)
      15.0_dp/16.0_dp,  &    ! (3,3)
      5.0_dp/16.0_dp /),   & ! (4,3)
      (/4, 3 /) )
  ! precomputed coefficients for the real valued quartic interpolation
  real(kind=dp), parameter, dimension(5,4) :: QUARTIC_LAGRANGE = &
    reshape( &
        (/  35.0_dp/128.0_dp, &    ! (1,1)
          140.0_dp/128.0_dp,  &    ! (2,1)
          - 70.0_dp/128.0_dp,   &    ! (3,1)
          28.0_dp/128.0_dp,   &    ! (4,1)
          - 5.0_dp/128.0_dp,   &    ! (5,1)
          - 5.0_dp/128.0_dp,   &    ! (1,2)
          60.0_dp/128.0_dp,   &    ! (2,2)
          90.0_dp/128.0_dp,   &    ! (3,2)
          - 20.0_dp/128.0_dp,   &    ! (4,2)
           3.0_dp/128.0_dp,   &    ! (5,2)
          - 5.0_dp/128.0_dp,   &    ! (1,3)
          28.0_dp/128.0_dp,   &    ! (2,3)
          - 70.0_dp/128.0_dp,  &    ! (3,3)
          140.0_dp/128.0_dp, & ! (4,3)
          + 35.0_dp/128.0_dp, & ! (5, 3)
          3.0_dp/128.0_dp,   &    ! (1,4)
          - 20.0_dp/128.0_dp,   &    ! (2,4)
          90.0_dp/128.0_dp,  &    ! (3,4)
          60.0_dp/128.0_dp, & ! (4,4)
          - 5.0_dp/128.0_dp /), & ! (5,4)
          (/ 5, 4 /)  )

  private
  public :: interpolationOperator

  contains

  subroutine interpolationOperator(orderInterpolation, NCoarseX, NCoarseY, &
    gridFunctionCoarse, gridFunctionFine)
    integer, intent(in) :: orderInterpolation, NCoarseX, NCoarseY
    integer :: highestOrder
    real(kind=dp), intent(in) :: gridFunctionCoarse(0:,0:)
    real(kind=dp), intent(out) :: gridFunctionFine(0:,0:)

    ! if the grid contains a sufficient number of points call the function.
    ! Otherwise use the highest possible order interpolation operator
    if (NCoarseX + 2 .ge. orderInterpolation .and. &
      NCoarseY + 2 .ge. orderInterpolation) then

      call prolongationOperator(orderInterpolation, NCoarseX, NCoarseY, &
        gridFunctionCoarse, gridFunctionFine)

    else

      call checkHighestOrderInterpolation(orderInterpolation, NCoarseX, &
        NCoarseY, highestOrder)

      call prolongationOperator(highestOrder, NCoarseX, NCoarseY, &
        gridFunctionCoarse, gridFunctionFine)

    end if

  end subroutine interpolationOperator

  subroutine checkHighestOrderInterpolation(orderInterpolation, &
    NCoarseX, NCoarseY, highestOrder)
    integer, intent(in) :: orderInterpolation, NCoarseX, NCoarseY
    integer, intent(out) :: highestOrder
    logical :: flag2, flag3, flag4, flag5

    flag2 = .false.
    flag3 = .false.
    flag4 = .false.
    flag5 = .false.

    if (orderInterpolation .eq. 5 .and. NCoarseX .ge. 3 .and. NCoarseY .ge. 3) then
      flag5 = .true.
    elseif (orderInterpolation .eq. 4 .and. NCoarseX .ge. 2 .and. NCoarseY .ge. 2) then
      flag4 = .true.
    elseif (orderInterpolation .eq. 3 .and. NCoarseX .ge. 1 .and. NCoarseY .ge. 1) then
      flag3 = .true.
    else
      flag2 = .true.
    end if

    if (flag5 .eqv. .true.) then
      highestOrder = 5
    elseif (flag4 .eqv. .true.) then
      highestOrder = 4
    elseif (flag3 .eqv. .true.) then
      highestOrder = 3
    else
      highestOrder = 2
    end if

  end subroutine checkHighestOrderInterpolation

  subroutine prolongationOperator(orderInterpolation, NCoarseX, NCoarseY, &
    gridFunctionCoarse, gridFunctionFine)
    integer, intent(in) :: orderInterpolation, NCoarseX, NCoarseY
    integer :: i, j, l, m, nodesX(orderInterpolation), nodesY(orderInterpolation)
    real(kind=dp), intent(in) :: gridFunctionCoarse(0:,0:)
    real(kind=dp), intent(out) :: gridFunctionFine(0:,0:)
    real(kind=dp) :: sum, weightsX(orderInterpolation), &
      weightsY(orderInterpolation), &
      weightsXY(orderInterpolation, orderInterpolation)

    ! on grid points where the fine and coarse grid overlap, simply
    ! copy the result
    do j = 0, NCoarseY + 1
      do i = 0, NCoarseX + 1

         gridFunctionFine(2*i,2*j) = gridFunctionCoarse(i,j)

      end do
    end do

    ! interpolation along the x direction where the coarse and fine grid
    ! overlap along the y-axis of the coarse grid
    do j = 1, NCoarseY
      do i = 0, NCoarseX

        call getNodesUnivariateInterpolation(i, NCoarseX, &
          orderInterpolation, nodesX)

        call getWeightsUnivariateInterpolation(i, NCoarseX, &
          orderInterpolation, weightsX)

        sum = 0.0_dp

        do m = 1, orderInterpolation
          sum = sum + weightsX(m)*gridFunctionCoarse(nodesX(m), j)
        end do

        gridFunctionFine(2*i + 1,2*j) = sum

      end do
    end do

    ! interpolation along the y direction where the coarse and fine grid
    ! overlap along the x-axis of the coarse grid
    do i = 1, NCoarseY
      do j = 0, NCoarseX

        call getNodesUnivariateInterpolation(j, NCoarseY, &
          orderInterpolation, nodesY)

        call getWeightsUnivariateInterpolation(j, NCoarseY, &
          orderInterpolation, weightsY)

        sum = 0.0_dp

        do l = 1, orderInterpolation
          sum = sum + weightsY(l)*gridFunctionCoarse(i, nodesY(l))
        end do

        gridFunctionFine(2*i, 2*j + 1) = sum

      end do
    end do

    ! interpolation along the portion of the fine grid which
    ! doesn't belong to the coarse grid

    do j = 0, NCoarseY
      do i = 0, NCoarseX

        call getNodesUnivariateInterpolation(i, NCoarseX, &
          orderInterpolation, nodesX)

        call getNodesUnivariateInterpolation(j, NCoarseY, &
          orderInterpolation, nodesY)

        call getWeightsBivariateInterpolation(i, j, NCoarseX, NCoarseY, &
          orderInterpolation, weightsXY)

        sum = 0.0_dp

        do l = 1, orderinterpolation
          do m = 1, orderInterpolation
            sum = sum + weightsXY(l, m)*gridFunctionCoarse(nodesX(l), nodesY(m))
          end do
        end do

        gridFunctionFine(2*i + 1, 2*j + 1) = sum

      end do
    end do

  end subroutine prolongationOperator

  subroutine getNodesBivariateInterpolation(i, j, NcoarseX, NcoarseY, &
    orderInterpolation, nodesInterpolationX, nodesInterpolationY)
    integer, intent(in) :: i, j, NcoarseX, NcoarseY, orderInterpolation
    integer, intent(out) ::  nodesInterpolationX(orderInterpolation), &
      nodesInterpolationY(orderInterpolation)
    integer :: loopX, loopY

    call getNodesUnivariateInterpolation(i, Ncoarse, &
      orderInterpolation, nodesInterpolationX)
    call getNodesUnivariateInterpolation(j, Ncoarse, &
      orderInterpolation, nodesInterpolationY)

  end subroutine getNodesBivariateInterpolation

  subroutine getWeightsBivariateInterpolation(i, j, NCoarseX, &
    NCoarseY, orderInterpolation, weightsInterpolation)
    integer, intent(in) :: i, j, NCoarseX, NCoarseY, orderInterpolation
    integer :: l, m
    real(kind=dp), intent(out) :: weightsInterpolation(:,:)
    real(kind=dp) :: weightsInterpolationX(orderInterpolation), &
      weightsInterpolationY(orderInterpolation)

    call getWeightsUnivariateInterpolation(i, NCoarseX, &
      orderInterpolation, weightsInterpolationX)

    call getWeightsUnivariateInterpolation(j, NCoarseY, &
      orderInterpolation, weightsInterpolationY)

    do m = 1, orderInterpolation
      do l = 1, orderInterpolation

        weightsInterpolation(l,m) = &
          weightsInterpolationX(l)*weightsInterpolationY(m)

      end do
    end do

  end subroutine getWeightsBivariateInterpolation

  subroutine getNodesUnivariateInterpolation(i, Ncoarse, &
    orderInterpolation, nodesInterpolation)
    integer, intent(in) :: i, Ncoarse, orderInterpolation
    integer, intent(out) :: nodesInterpolation(orderInterpolation)
    integer :: N

    N = Ncoarse

    select case (orderInterpolation)

      case (2)

        nodesInterpolation = [i, i + 1]

      case (3)

        if (i .eq. 0) then
          nodesInterpolation = [0, 1, 2]
        elseif ((i .eq. N - 1) .or. (i .eq. N)) then
          nodesInterpolation = [N - 1, N, N + 1]
        else
          nodesInterpolation = [i, i + 1, i + 2]
        end if


      case (4)

        if ((i .eq. 0) .or. (i .eq. 1)) then
          nodesInterpolation = [0, 1, 2, 3]
        elseif ((i .eq. N) .or. (i .eq. N - 1)) then
          nodesInterpolation = [N - 2, N - 1, N, N + 1]
        else
          nodesInterpolation = [i - 1, i, i + 1, i + 2]
        end if

      case (5)

        if ((i .eq. 0) .or. (i .eq. 1)) then
          nodesInterpolation = [0, 1, 2, 3, 4]
        elseif ((i .eq. N) .or. (i .eq. N - 1)) then
          nodesInterpolation = [N - 3, N - 2, N - 1, N, N + 1]
        else
          nodesInterpolation = [i - 2, i - 1, i, i + 1, i + 2]
        end if

    end select

  end subroutine getNodesUnivariateInterpolation

  subroutine getWeightsUnivariateInterpolation(i, Ncoarse, &
    orderInterpolation, weightsInterpolation)
    integer, intent(in) :: i, Ncoarse, orderInterpolation
    real(kind=dp), intent(out) :: weightsInterpolation(:)
    integer :: N

    N = Ncoarse

    select case (orderInterpolation)

      case (2)

        weightsInterpolation = LINEAR_LAGRANGE

      case (3)

        if (i .eq. N) then
          weightsInterpolation(:) = QUADRATIC_LAGRANGE(:,2)
        else
          weightsInterpolation(:) = QUADRATIC_LAGRANGE(:,1)
        end if

      case (4)

        if (i .eq. 0) then
          weightsInterpolation(:) = CUBIC_LAGRANGE(:,1)
        elseif (i .eq. N) then
          weightsInterpolation(:) = CUBIC_LAGRANGE(:,3)
        else
          weightsInterpolation(:) = CUBIC_LAGRANGE(:,2)
        end if

      case (5)

        if (i .eq. 0) then
          weightsInterpolation(:) = QUARTIC_LAGRANGE(:,1)
        elseif (i .eq. N) then
          weightsInterpolation(:) = QUARTIC_LAGRANGE(:,4)
        elseif (i .eq. N - 1) then
          weightsInterpolation(:) = QUARTIC_LAGRANGE(:,3)
        else
          weightsInterpolation(:) = QUARTIC_LAGRANGE(:,2)
        end if

    end select

  end subroutine getWeightsUnivariateInterpolation

end module FOCS2DInterpolationOperator
