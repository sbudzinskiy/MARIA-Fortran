!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the implementation of the [[maria_argcheck_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Implements the [[maria_argcheck_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
submodule (maria_argcheck_mod) maria_argcheck_sub
implicit none (type, external)

contains
    !------------------------------------------------------------------------------------------------------------------------

    module function iarg_is_bad &
    (bad_if, val, ref)
    use maria_reports_mod, only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        character(1), intent(in) :: bad_if
        integer,      intent(in) :: val
        integer,      intent(in) :: ref
        logical                  :: iarg_is_bad

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'IARG_IS_BAD'

    !-- Executable section -----------------------------------------------------
        select case (bad_if)
            case (BAD_IF_LESS)
                iarg_is_bad = (val < ref)
            case (BAD_IF_MORE)
                iarg_is_bad = (val > ref)
            case (BAD_IF_SAME)
                iarg_is_bad = (val == ref)
            case default
                call report_bad_arg(SRNAME, 1)
                iarg_is_bad = .true.
        end select
    end function iarg_is_bad

    !------------------------------------------------------------------------------------------------------------------------

    module function sarg_is_bad &
    (bad_if, val, ref)
    use maria_kinds_mod,      only: &
        WP => SP
    use maria_comparison_mod, only: &
        safe_less,                  &
        safe_eq
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        character(1), intent(in) :: bad_if
        real(WP),     intent(in) :: val
        real(WP),     intent(in) :: ref
        logical                  :: sarg_is_bad

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SARG_IS_BAD'

    !-- Executable section -----------------------------------------------------
        select case (bad_if)
            case (BAD_IF_LESS)
                sarg_is_bad = safe_less(val, ref)
            case (BAD_IF_MORE)
                sarg_is_bad = safe_less(ref, val)
            case (BAD_IF_SAME)
                sarg_is_bad = safe_eq(val, ref)
            case default
                call report_bad_arg(SRNAME, 1)
                sarg_is_bad = .true.
        end select
    end function sarg_is_bad

    !------------------------------------------------------------------------------------------------------------------------

    module function darg_is_bad &
    (bad_if, val, ref)
    use maria_kinds_mod,      only: &
        WP => DP
    use maria_comparison_mod, only: &
        safe_less,                  &
        safe_eq
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        character(1), intent(in) :: bad_if
        real(WP),     intent(in) :: val
        real(WP),     intent(in) :: ref
        logical                  :: darg_is_bad

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DARG_IS_BAD'

    !-- Executable section -----------------------------------------------------
        select case (bad_if)
            case (BAD_IF_LESS)
                darg_is_bad = safe_less(val, ref)
            case (BAD_IF_MORE)
                darg_is_bad = safe_less(ref, val)
            case (BAD_IF_SAME)
                darg_is_bad = safe_eq(val, ref)
            case default
                call report_bad_arg(SRNAME, 1)
                darg_is_bad = .true.
        end select
    end function darg_is_bad

    !------------------------------------------------------------------------------------------------------------------------
end submodule maria_argcheck_sub
