!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_comparison_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Provides access to safe comparison subroutines for floating-point numbers.
!----------------------------------------------------------------------------------------------------------------------------
module maria_comparison_mod
implicit none (type, external)

! Interfaces
public :: safe_less, &
          safe_leq,  &
          safe_eq
! Procedures
public :: dsafe_less, &
          ssafe_less, &
          dsafe_leq,  &
          ssafe_leq,  &
          dsafe_eq,   &
          ssafe_eq
private

!> Checks if one floating-point number is strictly less than the other one, using unit roundoff.
interface safe_less
    procedure dsafe_less
    procedure ssafe_less
end interface safe_less

!> Checks if one floating-point number is less than or equal to the other one, using unit roundoff.
interface safe_leq
    procedure dsafe_leq
    procedure ssafe_leq
end interface safe_leq

!> Checks if one floating-point number is equal to the other one, using unit roundoff.
interface safe_eq
    procedure dsafe_eq
    procedure ssafe_eq
end interface safe_eq

interface
    !------------------------------------------------------------------------------------------------------------------------

    !> Checks if one double-precision number is strictly less than the other one, using unit roundoff.
    !>
    !> Returns `.true.` if \( a - b \leq -\textbf{u} \) and `.false.` otherwise.
    module function dsafe_less &
    (a, b)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> First number
        real(WP), intent(in) :: a
        !> Second number
        real(WP), intent(in) :: b
        logical              :: dsafe_less
    end function dsafe_less

    !------------------------------------------------------------------------------------------------------------------------

    !> Checks if one single-precision number is strictly less than the other one, using unit roundoff.
    !>
    !> Returns `.true.` if \( a - b \leq -\textbf{u} \) and `.false.` otherwise.
    module function ssafe_less &
    (a, b)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> First number
        real(WP), intent(in) :: a
        !> Second number
        real(WP), intent(in) :: b
        logical              :: ssafe_less
    end function ssafe_less

    !------------------------------------------------------------------------------------------------------------------------

    !> Checks if one double-precision number is less than or equal to the other one, using unit roundoff.
    !>
    !> Returns `.true.` if \( a - b < \textbf{u} \) and `.false.` otherwise.
    module function dsafe_leq &
    (a, b)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> First number
        real(WP), intent(in) :: a
        !> Second number
        real(WP), intent(in) :: b
        logical              :: dsafe_leq
    end function dsafe_leq

    !------------------------------------------------------------------------------------------------------------------------

    !> Checks if one single-precision number is less than or equal to the other one, using unit roundoff.
    !>
    !> Returns `.true.` if \( a - b < \textbf{u} \) and `.false.` otherwise.
    module function ssafe_leq &
    (a, b)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> First number
        real(WP), intent(in) :: a
        !> Second number
        real(WP), intent(in) :: b
        logical              :: ssafe_leq
    end function ssafe_leq

    !------------------------------------------------------------------------------------------------------------------------

    !> Checks if one double-precision number is equal to the other one, using unit roundoff.
    !>
    !> Returns `.true.` if \( |a - b| < \textbf{u} \) and `.false.` otherwise.
    module function dsafe_eq &
    (a, b)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> First number
        real(WP), intent(in) :: a
        !> Second number
        real(WP), intent(in) :: b
        logical              :: dsafe_eq
    end function dsafe_eq

    !------------------------------------------------------------------------------------------------------------------------

    !> Checks if one single-precision number is equal to the other one, using unit roundoff.
    !>
    !> Returns `.true.` if \( |a - b| < \textbf{u} \) and `.false.` otherwise.
    module function ssafe_eq &
    (a, b)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> First number
        real(WP), intent(in) :: a
        !> Second number
        real(WP), intent(in) :: b
        logical              :: ssafe_eq
    end function ssafe_eq

    !------------------------------------------------------------------------------------------------------------------------
end interface

end module maria_comparison_mod
