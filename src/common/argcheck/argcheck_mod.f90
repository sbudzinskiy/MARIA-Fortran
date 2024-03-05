!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_argcheck_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Provides access to subroutines for checking the sanity of input arguments.
!----------------------------------------------------------------------------------------------------------------------------
module maria_argcheck_mod
implicit none (type, external)

! Parameters
public :: BAD_IF_LESS, &
          BAD_IF_MORE, &
          BAD_IF_SAME
! Interfaces
public :: arg_is_bad
! Procedures
public :: iarg_is_bad, &
          sarg_is_bad, &
          darg_is_bad
private

!> Signifies that an argument cannot be smaller than a reference value.
character(1), parameter :: BAD_IF_LESS = 'L'
!> Signifies that an argument cannot be greater than a reference value.
character(1), parameter :: BAD_IF_MORE = 'M'
!> Signifies that an argument cannot be equal to a reference value.
character(1), parameter :: BAD_IF_SAME = 'S'

!> Compares a numerical value with a reference value.
interface arg_is_bad
    procedure iarg_is_bad
    procedure sarg_is_bad
    procedure darg_is_bad
end interface arg_is_bad

interface
    !------------------------------------------------------------------------------------------------------------------------

    !> Compares an input integer value with a reference value.
    !>
    !> Returns `compare(val,ref)`, where
    !>
    !> - `compare(val,ref)` is `val < ref`, or
    !> - `compare(val,ref)` is `val > ref`, or
    !> - `compare(val,ref)` is `val = ref`.
    module function iarg_is_bad &
    (bad_if, val, ref)
    implicit none
        !> Specifies how to compare the two values:
        !> 
        !> - `bad_if = [[maria_argcheck_mod(module):BAD_IF_LESS(variable)]]` then `compare` is `<`
        !> - `bad_if = [[maria_argcheck_mod(module):BAD_IF_MORE(variable)]]` then `compare` is `>`
        !> - `bad_if = [[maria_argcheck_mod(module):BAD_IF_SAME(variable)]]` then `compare` is `=`
        character(1), intent(in) :: bad_if
        !> Input value
        integer,      intent(in) :: val
        !> Reference value
        integer,      intent(in) :: ref
        logical                  :: iarg_is_bad
    end function iarg_is_bad

    !------------------------------------------------------------------------------------------------------------------------

    !> Compares an input single-precision value with a reference value.
    !>
    !> Returns `compare(val,ref)`, where
    !>
    !> - `compare(val,ref)` is `val < ref`, or
    !> - `compare(val,ref)` is `val > ref`, or
    !> - `compare(val,ref)` is `val = ref`.
    module function sarg_is_bad &
    (bad_if, val, ref)
    use maria_kinds_mod, only: & 
        WP => SP
    implicit none
        !> Specifies how to compare the two values:
        !> 
        !> - `bad_if = [[maria_argcheck_mod(module):BAD_IF_LESS(variable)]]` then `compare` is `<`
        !> - `bad_if = [[maria_argcheck_mod(module):BAD_IF_MORE(variable)]]` then `compare` is `>`
        !> - `bad_if = [[maria_argcheck_mod(module):BAD_IF_SAME(variable)]]` then `compare` is `=`
        character(1), intent(in) :: bad_if
        !> Input value
        real(WP),     intent(in) :: val
        !> Reference value
        real(WP),     intent(in) :: ref
        logical                  :: sarg_is_bad
    end function sarg_is_bad

    !------------------------------------------------------------------------------------------------------------------------

    !> Compares an input double-precision value with a reference value.
    !>
    !> Returns `compare(val,ref)`, where
    !>
    !> - `compare(val,ref)` is `val < ref`, or
    !> - `compare(val,ref)` is `val > ref`, or
    !> - `compare(val,ref)` is `val = ref`.
    module function darg_is_bad &
    (bad_if, val, ref)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Specifies how to compare the two values:
        !> 
        !> - `bad_if = [[maria_argcheck_mod(module):BAD_IF_LESS(variable)]]` then `compare` is `<`
        !> - `bad_if = [[maria_argcheck_mod(module):BAD_IF_MORE(variable)]]` then `compare` is `>`
        !> - `bad_if = [[maria_argcheck_mod(module):BAD_IF_SAME(variable)]]` then `compare` is `=`
        character(1), intent(in) :: bad_if
        !> Input value
        real(WP),     intent(in) :: val
        !> Reference value
        real(WP),     intent(in) :: ref
        logical                  :: darg_is_bad
    end function darg_is_bad

    !------------------------------------------------------------------------------------------------------------------------
end interface

end module maria_argcheck_mod
