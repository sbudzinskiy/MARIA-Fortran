!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_reports_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Provides access to subroutines for reporting errors related to input arguments and runtime errors.
!----------------------------------------------------------------------------------------------------------------------------
module maria_reports_mod
implicit none (type, external)

! Procedures
public :: report_bad_arg,    &
          report_runtime_err
! Varibales
public :: SINGULAR_MATRIX_ERR_CODE,   &
          MAX_SEED_EXCEEDED_ERR_CODE, & 
          IND_OUT_OF_RANGE_ERR_CODE,  &
          SVD_FAILED_ERR_CODE
private

!> Signifies a runtime error due to inverting a singular matrix.
integer, parameter :: SINGULAR_MATRIX_ERR_CODE = 1
integer, parameter :: MAX_SEED_EXCEEDED_ERR_CODE = 2
integer, parameter :: IND_OUT_OF_RANGE_ERR_CODE = 3
integer, parameter :: SVD_FAILED_ERR_CODE = 4

interface
    !------------------------------------------------------------------------------------------------------------------------

    !> Prints an error message, saying that an argument at position `pos`
    !> of a subroutine called `srname` had an illegal value.
    module subroutine report_bad_arg &
    (srname, pos)
    implicit none
        !> Name of the subroutine
        character(*), intent(in) :: srname
        !> Position of the argument
        integer,      intent(in) :: pos
    end subroutine report_bad_arg

    !------------------------------------------------------------------------------------------------------------------------

    !> Prints an error message, saying that a specific runtime error was encountered in a subroutine called `srname`.
    module subroutine report_runtime_err &
    (srname, err_code)
    implicit none
        !> Name of the subroutine
        character(*), intent(in) :: srname
        !> Error code:
        !>
        !> - `err_code = [[maria_reports_mod(module):SINGULAR_MATRIX_ERR_CODE(variable)]]`
        !> if attempted to invert a singular matrix
        integer,      intent(in) :: err_code
    end subroutine report_runtime_err

    !------------------------------------------------------------------------------------------------------------------------
end interface

end module maria_reports_mod
