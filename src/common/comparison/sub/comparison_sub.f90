!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the implementation of the [[maria_comparison_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Implements the [[maria_comparison_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
submodule (maria_comparison_mod) maria_comparison_sub
implicit none (type, external)

contains
    !------------------------------------------------------------------------------------------------------------------------

    module function dsafe_less &
    (a, b)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        EPS => D_MACHTOL
    !-- Input/output arguments -------------------------------------------------
        real(WP), intent(in) :: a
        real(WP), intent(in) :: b
        logical              :: dsafe_less
        
    !-- Executable section -----------------------------------------------------
        dsafe_less = (a - b <= -EPS)
    end function dsafe_less

    !------------------------------------------------------------------------------------------------------------------------

    module function ssafe_less &
    (a, b)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        EPS => S_MACHTOL
    !-- Input/output arguments -------------------------------------------------
        real(WP), intent(in) :: a
        real(WP), intent(in) :: b
        logical              :: ssafe_less

    !-- Executable section -----------------------------------------------------
        ssafe_less = (a - b <= -EPS)
    end function ssafe_less

    !------------------------------------------------------------------------------------------------------------------------

    module function dsafe_leq &
    (a, b)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        EPS => D_MACHTOL
    !-- Input/output arguments -------------------------------------------------
        real(WP), intent(in) :: a
        real(WP), intent(in) :: b
        logical              :: dsafe_leq

    !-- Executable section -----------------------------------------------------
        dsafe_leq = (a - b < EPS)
    end function dsafe_leq

    !------------------------------------------------------------------------------------------------------------------------

    module function ssafe_leq &
    (a, b)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        EPS => S_MACHTOL
    !-- Input/output arguments -------------------------------------------------
        real(WP), intent(in) :: a
        real(WP), intent(in) :: b
        logical              :: ssafe_leq

    !-- Executable section -----------------------------------------------------
        ssafe_leq = (a - b < EPS)
    end function ssafe_leq

    !------------------------------------------------------------------------------------------------------------------------

    module function dsafe_eq &
    (a, b)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        EPS => D_MACHTOL
    !-- Input/output arguments -------------------------------------------------
        real(WP), intent(in) :: a
        real(WP), intent(in) :: b
        logical              :: dsafe_eq

    !-- Executable section -----------------------------------------------------
        dsafe_eq = (abs(a - b) < EPS)
    end function dsafe_eq

    !------------------------------------------------------------------------------------------------------------------------

    module function ssafe_eq &
    (a, b)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        EPS => S_MACHTOL
    !-- Input/output arguments -------------------------------------------------
        real(WP), intent(in) :: a
        real(WP), intent(in) :: b
        logical              :: ssafe_eq

    !-- Executable section -----------------------------------------------------
        ssafe_eq = (abs(a - b) < EPS)
    end function ssafe_eq

    !------------------------------------------------------------------------------------------------------------------------
end submodule maria_comparison_sub
