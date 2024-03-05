!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the implementation of the [[maria_assert_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Implements the [[maria_assert_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
submodule (maria_assert_mod) maria_assert_sub
implicit none (type, external)

contains
    !------------------------------------------------------------------------------------------------------------------------

    module subroutine assert_linefile &
    (flag, filename, line_number)
    use, intrinsic :: iso_fortran_env, only: &
        STDERR => error_unit
    !-- Input/output arguments -------------------------------------------------
        logical,      intent(in) :: flag
        character(*), intent(in) :: filename
        integer,      intent(in) :: line_number

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: OUT_FORMAT = "('Assertion failed in [', a, '] on line ', i0)"
    
    !-- Executable section -----------------------------------------------------
        if (.not. flag) then
            write(STDERR, OUT_FORMAT) filename, line_number
            stop 1
        end if
    end subroutine assert_linefile

    !------------------------------------------------------------------------------------------------------------------------
end submodule maria_assert_sub
