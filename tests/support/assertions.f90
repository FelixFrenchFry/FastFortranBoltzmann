module assertions
    ! imports
    use iso_fortran_env, only: error_unit
    use ieee_arithmetic, only: ieee_is_finite
    use settings, only: FP
    implicit none

    private

    public :: check
    public :: check_close

contains

    subroutine check( &
        condition, name &
        )
        ! inputs
        logical, intent(in) :: condition
        character(len=*), intent(in) :: name

        if (.not. condition) then
            write(error_unit, '(A)') "failed: " // name
            error stop 1
        end if
    end subroutine check


    subroutine check_close( &
        actual, expected, name &
        )
        ! inputs
        real(FP), intent(in) :: actual
        real(FP), intent(in) :: expected
        character(len=*), intent(in) :: name

        ! temp
        real(FP) :: tolerance

        tolerance = 128.0_FP * epsilon(1.0_FP) * max(1.0_FP, abs(expected))

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected) .or. &
            abs(actual - expected) > tolerance) then
            write(error_unit, '(A)') "failed: " // name
            write(error_unit, *) "actual, expected, tolerance: ", actual, expected, tolerance
            error stop 1
        end if
    end subroutine check_close


end module assertions
