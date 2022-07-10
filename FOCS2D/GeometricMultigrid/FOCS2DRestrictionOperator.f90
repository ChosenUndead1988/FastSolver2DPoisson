module FOCS2DRestrictionOperator

  use FOCS2DSpecifyPrecision
  use FOCS2DConstants

  real(dp), parameter :: SCALING_RESTRICTION = 4.0_dp

  private
  public :: RestrictionOperator

  contains

  subroutine RestrictionOperator(orderRestriction, NCoarseX, NCoarseY, &
    gridFunctionFine, restrictedGridFunction)
    integer, intent(in) :: orderRestriction, NCoarseX, NCoarseY
    real(kind=dp), intent(in) :: gridFunctionFine(0:,0:)
    real(kind=dp), intent(out) :: restrictedGridFunction(0:,0:)

    if (orderRestriction .eq. 0) then

      call InjectionOperator(NCoarseX, NCoarseY, &
        gridFunctionFine, restrictedGridFunction)

    elseif (orderRestriction .eq. 2) then

      call FullWeightingOperator(NCoarseX, NCoarseY, gridFunctionFine, &
        restrictedGridFunction)

    else

      call FullWeightingOperator(NCoarseX, NCoarseY, gridFunctionFine, &
        restrictedGridFunction)

    end if

  end subroutine RestrictionOperator

  subroutine FullWeightingOperator(NCoarseX, NCoarseY, gridFunctionFine, &
    restrictedGridFunction)
    integer, intent(in) :: NCoarseX, NCoarseY
    integer :: i, j, l, m
    real(dp), intent(in) :: gridFunctionFine(0:,0:)
    real(dp), intent(out) :: restrictedGridFunction(0:,0:)

    restrictedGridFunction = 0.0_dp

    do m = 1, NCoarseY
      do l = 1, NCoarseX

        i = 2*l
        j = 2*m

        restrictedGridFunction(l,m) = SCALING_RESTRICTION*(1.0_dp/16.0_dp)*( &
          4.0_dp*gridFunctionFine(i,j) &
          + 2.0_dp*(gridFunctionFine(i,j-1) + gridFunctionFine(i-1,j) &
          + gridFunctionFine(i+1,j) + gridFunctionFine(i,j+1)) &
          + gridFunctionFine(i-1,j-1) + gridFunctionFine(i+1,j-1) &
          + gridFunctionFine(i-1,j+1) + gridFunctionFine(i+1,j+1))

      end do
    end do

  end subroutine FullWeightingOperator

  subroutine InjectionOperator(NCoarseX, NCoarseY, &
    gridFunctionFine, restrictedGridFunction)
    integer, intent(in) :: NCoarseX, NCoarseY
    integer :: i, j
    real(dp), intent(in) :: gridFunctionFine(0:,0:)
    real(dp), intent(out) :: restrictedGridFunction(0:,0:)

    restrictedGridFunction = 0.0_dp

    do j = 1, NCoarseY
      do i = 1, NCoarseX

        restrictedGridFunction(i, j) = SCALING_RESTRICTION*(&
          gridFunctionFine(2*i,2*j))

      end do
    end do

  end subroutine InjectionOperator

end module FOCS2DRestrictionOperator
