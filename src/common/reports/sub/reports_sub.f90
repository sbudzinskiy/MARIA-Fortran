!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the implementation of the [[maria_reports_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Implements the [[maria_reports_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
submodule (maria_reports_mod) maria_reports_sub
implicit none (type, external)

contains
    !------------------------------------------------------------------------------------------------------------------------

    module subroutine report_bad_arg &
    (srname, pos)
    use, intrinsic :: iso_fortran_env, only: &
        STDERR => error_unit
    !-- Input/output arguments -------------------------------------------------
        character(*), intent(in) :: srname
        integer,      intent(in) :: pos

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: OUT_FORMAT = "(a, ': parameter ', i0, ' had an illegal value')"

    !-- Executable section -----------------------------------------------------
        write(STDERR, OUT_FORMAT) srname, pos
    end subroutine report_bad_arg

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine report_runtime_err &
    (srname, err_code)
    use, intrinsic :: iso_fortran_env, only: &
        STDERR => error_unit
    !-- Input/output arguments -------------------------------------------------
        character(*), intent(in) :: srname
        integer,      intent(in) :: err_code

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: OUT_FORMAT = "(a, ': ', a)"

    !-- Executable section -----------------------------------------------------
        select case (err_code)
            case (SINGULAR_MATRIX_ERR_CODE)
                write(STDERR, OUT_FORMAT) srname, "trying to invert a singular matrix"
            case (MAX_SEED_EXCEEDED_ERR_CODE)
                write(STDERR, OUT_FORMAT) srname, "the seed size required by the system exceeds MAX_SEED_SIZE"
            case (IND_OUT_OF_RANGE_ERR_CODE)
                write(STDERR, OUT_FORMAT) srname, "index of an array is out of range"
            case (SVD_FAILED_ERR_CODE)
                write(STDERR, OUT_FORMAT) srname, "the computational routine for SVD did not converge"
            case default   
                write(STDERR, OUT_FORMAT) srname, "unknown runtime error occured"
        end select 
    end subroutine report_runtime_err

    !------------------------------------------------------------------------------------------------------------------------
end submodule maria_reports_sub
