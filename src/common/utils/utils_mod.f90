!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_utils_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Provides access to general utility subroutines.
!----------------------------------------------------------------------------------------------------------------------------
module maria_utils_mod
use maria_kinds_mod, only: &
    SP,                    &
    DP
implicit none (type, external)

! Interfaces
public :: optional_val,          &
          are_close,             &
          swap_pair,             &
          linspace,              &
          arange,                &
          median_of_3,           &
          sort,                  &
          select_no_replacement, &
          permute
! Procedures
public :: loptional_val,          &
          ioptional_val,          &
          soptional_val,          &
          doptional_val,          &
          sare_close,             &
          dare_close,             &
          iswap_pair,             &
          sswap_pair,             &
          dswap_pair,             &
          slinspace,              &
          dlinspace,              &
          iarange,                &
          sarange,                &
          darange,                &
          imedian_of_3,           &
          isort,                  &
          iselect_no_replacement, &
          ipermute,               &
          srealloc,               &
          drealloc
private

!> Processes the value of an optional argument:
!>
!> - if `val` is present, returns `val`,
!> - if `val` is missing, returns `default_val`.
interface optional_val
    procedure loptional_val
    procedure ioptional_val
    procedure soptional_val
    procedure doptional_val
end interface optional_val

!> Checks if two floating-point numbers are approximately equal 
!> with given absolute and/or relative tolerance.
interface are_close
    procedure sare_close
    procedure dare_close
end interface are_close

!> Swaps the values of two variables.
interface swap_pair
    procedure iswap_pair
    procedure sswap_pair
    procedure dswap_pair
end interface swap_pair

!> Generates an array of equidistant points in a given interval.
interface linspace
    procedure slinspace
    procedure dlinspace
end interface linspace

!> Generates an array of equidistant points with a given step size.
interface arange
    procedure iarange
    procedure sarange
    procedure darange
end interface arange

!> Computes the median of three numbers. 
interface median_of_3
    procedure imedian_of_3
end interface median_of_3

!> Sorts an array of numbers in increasing or decreasing order.
!>
!> Uses quicksort with median-of-3 pivoting and switches to insertion sort for smaller arrays.
interface sort
    procedure isort
end interface sort

!> Selects multiple random items without replacement from a given array.
interface select_no_replacement
    procedure iselect_no_replacement
end interface select_no_replacement

!> Randomly permutes the elements of an array,
interface permute
    procedure ipermute
end interface permute

interface
    !------------------------------------------------------------------------------------------------------------------------

    !> Processes the value of an optional logical argument:
    !>
    !> - if `val` is present, returns `val`,
    !> - if `val` is missing, returns `default_val`.
    module function loptional_val &
    (default_val, val)
    implicit none
        !> Default value of the optional argument
        logical, intent(in)           :: default_val
        !> Optional argument
        logical, intent(in), optional :: val
        logical                       :: loptional_val
    end function loptional_val

    !------------------------------------------------------------------------------------------------------------------------

    !> Processes the value of an optional integer argument:
    !>
    !> - if `val` is present, returns `val`,
    !> - if `val` is missing, returns `default_val`.
    module function ioptional_val &
    (default_val, val)
    implicit none
        !> Default value of the optional argument
        integer, intent(in)           :: default_val
        !> Optional argument
        integer, intent(in), optional :: val
        integer                       :: ioptional_val
    end function ioptional_val

    !------------------------------------------------------------------------------------------------------------------------

    !> Processes the value of an optional single-precision argument:
    !>
    !> - if `val` is present, returns `val`,
    !> - if `val` is missing, returns `default_val`.
    module function soptional_val &
    (default_val, val)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Default value of the optional argument
        real(WP), intent(in)           :: default_val
        !> Optional argument
        real(WP), intent(in), optional :: val
        real(WP)                       :: soptional_val
    end function soptional_val

    !------------------------------------------------------------------------------------------------------------------------

    !> Processes the value of an optional double-precision argument:
    !>
    !> - if `val` is present, returns `val`,
    !> - if `val` is missing, returns `default_val`.
    module function doptional_val &
    (default_val, val)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Default value of the optional argument
        real(WP), intent(in)           :: default_val
        !> Optional argument
        real(WP), intent(in), optional :: val
        real(WP)                       :: doptional_val
    end function doptional_val

    !------------------------------------------------------------------------------------------------------------------------

    !> Checks if two single-precision numbers are approximately equal 
    !> with given absolute and/or relative tolerance:
    !>
    !> \[ |a - b| < \max\big( \text{atol}, \text{rtol} \cdot \max(|a|, |b|)  \big) \].
    module function sare_close &
    (a, b, info, atol, rtol)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> First number
        real(WP), intent(in)           :: a
        !> Second number
        real(WP), intent(in)           :: b
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)          :: info
        !> Absolute tolerance
        !>
        !> Possible values: \( \text{atol} \geq 0 \)
        !>
        !> Default value: [[maria_constants_mod(module):S_MACHTOL(variable)]]
        real(WP), intent(in), optional :: atol
        !> Relative tolerance
        !>
        !> Possible values: \( \text{rtol} \geq 0 \)
        !>
        !> Default value: \( 0 \)
        real(WP), intent(in), optional :: rtol
        logical                        :: sare_close
    end function sare_close

    !------------------------------------------------------------------------------------------------------------------------

    !> Checks if two double-precision numbers are approximately equal 
    !> with given absolute and/or relative tolerance:
    !>
    !> \[ |a - b| < \max\big( \text{atol}, \text{rtol} \cdot \max(|a|, |b|)  \big) \].
    module function dare_close &
    (a, b, info, atol, rtol)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> First number
        real(WP), intent(in)           :: a
        !> Second number
        real(WP), intent(in)           :: b
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)          :: info
        !> Absolute tolerance
        !>
        !> Possible values: \( \text{atol} \geq 0 \)
        !>
        !> Default value: [[maria_constants_mod(module):D_MACHTOL(variable)]]
        real(WP), intent(in), optional :: atol
        !> Relative tolerance
        !>
        !> Possible values: \( \text{rtol} \geq 0 \)
        !>
        !> Default value: \( 0 \)
        real(WP), intent(in), optional :: rtol
        logical                        :: dare_close
    end function dare_close

    !------------------------------------------------------------------------------------------------------------------------

    !> Swaps the values of two integer variables.
    module subroutine iswap_pair &
    (a, b)
    implicit none
        !> On exit, is overwritten with the value of `b`
        integer, intent(inout) :: a
        !> On exit, is overwritten with the value of `a`
        integer, intent(inout) :: b
    end subroutine iswap_pair

    !------------------------------------------------------------------------------------------------------------------------

    !> Swaps the values of two single-precision variables.
    module subroutine sswap_pair &
    (a, b)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> On exit, is overwritten with the value of `b`
        real(WP), intent(inout) :: a
        !> On exit, is overwritten with the value of `a`
        real(WP), intent(inout) :: b
    end subroutine sswap_pair

    !------------------------------------------------------------------------------------------------------------------------

    !> Swaps the values of two double-precision variables.
    module subroutine dswap_pair &
    (a, b)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> On exit, is overwritten with the value of `b`
        real(WP), intent(inout) :: a
        !> On exit, is overwritten with the value of `a`
        real(WP), intent(inout) :: b
    end subroutine dswap_pair

    !------------------------------------------------------------------------------------------------------------------------

    !> Generates an array of \( n \) equidistant points in an interval \( [a, b] \).
    !>
    !> The end point \( b \) can optionally be excluded based on the value of `include_b`.
    module subroutine slinspace &
    (n, x, a, b, info, include_b, step)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Number of points to generate
        !>
        !> Possible values: \( n \geq 0 \)
        integer,  intent(in)              :: n
        !> Array of size \( n \)
        !>
        !> On exit, contains the equidistant points such that
        !> \[
        !>      x(j) = a + (j-1)h, \quad j = 1, \ldots, n,
        !> \]
        !> where 
        !> \[
        !>      h = 
        !>      \begin{cases}
        !>          \frac{b - a}{n - 1}, & `include_b = .true.` \\
        !>          \frac{b - a}{n}, & `include_b = .false.`
        !>      \end{cases}
        !> \]
        real(WP), intent(out), contiguous :: x(:)
        !> Left boundary 
        real(WP), intent(in)              :: a
        !> Right boundary
        !>
        !> Possible values: \( b > a \)
        real(WP), intent(in)              :: b
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)             :: info
        !> Specifies if \( b \) is included or excluded.
        !>
        !> Default value: `.true.` 
        logical,  intent(in),  optional   :: include_b
        !> Distance between adjacent points
        real(WP), intent(out), optional   :: step
    end subroutine slinspace

    !------------------------------------------------------------------------------------------------------------------------

    !> Generates an array of \( n \) equidistant points in an interval \( [a, b] \).
    module subroutine dlinspace &
    (n, x, a, b, info, include_b, step)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Number of points to generate
        !>
        !> Possible values: \( n \geq 0 \)
        integer,  intent(in)              :: n
        !> Array of size \( n \)
        !>
        !> On exit, contains the equidistant points such that
        !> \[
        !>      x(j) = a + (j-1)h, \quad j = 1, \ldots, n,
        !> \]
        !> where 
        !> \[
        !>      h = 
        !>      \begin{cases}
        !>          \frac{b - a}{n - 1}, & `include_b = .true.` \\
        !>          \frac{b - a}{n}, & `include_b = .false.`
        !>      \end{cases}
        !> \]
        real(WP), intent(out), contiguous :: x(:)
        !> Left boundary 
        real(WP), intent(in)              :: a
        !> Right boundary
        !>
        !> Possible values: \( b > a \)
        real(WP), intent(in)              :: b
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)             :: info
        !> Specifies if \( b \) is included or excluded.
        !>
        !> Default value: `.true.` 
        logical,  intent(in),  optional   :: include_b
        !> Distance between adjacent points
        real(WP), intent(out), optional   :: step
    end subroutine dlinspace

    !------------------------------------------------------------------------------------------------------------------------

    !> Generates an array of \( n \) equidistant points with a given step size.
    module subroutine iarange &
    (n, x, a, step, info)
    implicit none
        !> Number of points to generate
        !>
        !> Possible values: \( n \geq 0 \)
        integer, intent(in)              :: n
        !> Array of size \( n \)
        !>
        !> On exit, contains the equidistant points such that
        !> \[
        !>      x(j) = a + (j-1) \text{step}, \quad j = 1, \ldots, n,
        !> \]
        integer, intent(out), contiguous :: x(:)
        !> Starting point
        integer, intent(in)              :: a
        !> Step
        integer, intent(in)              :: step
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)            :: info
    end subroutine iarange

    !------------------------------------------------------------------------------------------------------------------------

    !> Generates an array of \( n \) equidistant points with a given step size.
    module subroutine sarange &
    (n, x, a, step, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Number of points to generate
        !>
        !> Possible values: \( n \geq 0 \)
        integer,  intent(in)              :: n
        !> Array of size \( n \)
        !>
        !> On exit, contains the equidistant points such that
        !> \[
        !>      x(j) = a + (j-1) \text{step}, \quad j = 1, \ldots, n,
        !> \]
        real(WP), intent(out), contiguous :: x(:)
        !> Starting point
        real(WP), intent(in)              :: a
        !> Step
        real(WP), intent(in)              :: step
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)             :: info
    end subroutine sarange

    !------------------------------------------------------------------------------------------------------------------------

    !> Generates an array of \( n \) equidistant points with a given step size.
    module subroutine darange &
    (n, x, a, step, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Number of points to generate
        !>
        !> Possible values: \( n \geq 0 \)
        integer,  intent(in)              :: n
        !> Array of size \( n \)
        !>
        !> On exit, contains the equidistant points such that
        !> \[
        !>      x(j) = a + (j-1) \text{step}, \quad j = 1, \ldots, n,
        !> \]
        real(WP), intent(out), contiguous :: x(:)
        !> Starting point
        real(WP), intent(in)              :: a
        !> Step
        real(WP), intent(in)              :: step
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)             :: info
    end subroutine darange

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the median of three integer numbers. 
    module integer function imedian_of_3 &
    (a, b, c)
    implicit none
        !> First number
        integer, intent(in) :: a
        !> Second number
        integer, intent(in) :: b
        !> Third number
        integer, intent(in) :: c 
    end function imedian_of_3

    !------------------------------------------------------------------------------------------------------------------------

    !> Sorts an array of integer numbers in increasing or decreasing order.
    !>
    !> Uses quicksort with median-of-3 pivoting and switches to insertion sort for smaller arrays.
    module subroutine isort &
    (id, n, x, info)
    implicit none
        !> Specifies the sorting order:
        !>
        !> - `id = "I" or "i"`: sort in increasing order
        !> - `id`= "D" or "d"`: sort in decreasing order
        character(1), intent(in)                :: id
        !> Size of the array.
        !>
        !> Possible values: \( n \geq 0 \)
        integer,      intent(in)                :: n
        !> Array to be sorted.
        !>
        !> On exit, the array sorted in increasing or decreasing order
        integer,      intent(inout), contiguous :: x(:)
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,      intent(out)               :: info
    end subroutine isort

    !------------------------------------------------------------------------------------------------------------------------

    !> Selects \( k \) random items without replacement from a given integer array.
    module subroutine iselect_no_replacement &
    (rng, n, x, k, info)
    use maria_prng_mod, only: &
        prng
    implicit none
        !> Pseudorandom number generator
        class(prng), intent(in)                :: rng
        !> Size of the array
        !>
        !> Possible values: \( n \geq 0 \)
        integer,     intent(in)                :: n
        !> Array of size \( n \)
        !>
        !> On entry, contains the values that are selected from.
        !>
        !> On exit, a permutation of the initial array such that the first `k` items
        !> have been selected at random without replacement.
        integer,     intent(inout), contiguous :: x(:)
        !> Number of items to select
        !>
        !> Possible values: \( 0 \leq k \leq \max(n-1, 0) \)
        integer,     intent(in)                :: k
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,     intent(out)               :: info
    end subroutine iselect_no_replacement

    !------------------------------------------------------------------------------------------------------------------------

    !> Randomly permutes the entries of a given integer array.
    module subroutine ipermute &
    (rng, n, x, info)
    use maria_prng_mod, only: &
        prng
    implicit none
        !> Pseudorandom number generator
        class(prng), intent(in)                :: rng
        !> Size of the array
        !>
        !> Possible values: \( n \geq 0 \)
        integer,     intent(in)                :: n
        !> Array of size \( n \)
        !>
        !> On entry, contains the values that are permuted.
        !>
        !> On exit, a random permutation of the initial array.
        integer,     intent(inout), contiguous :: x(:)
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,     intent(out)               :: info
    end subroutine ipermute

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine srealloc &
    (copy, x, n, ierr, factor)
    use maria_kinds_mod,    only:   &
        WP => SP
    implicit none
        logical,    intent(in   )               :: copy
        real(WP),   intent(inout),  allocatable :: x(:)
        integer,    intent(in   )               :: n
        integer,    intent(  out)               :: ierr
        integer,    intent(in   ),  optional    :: factor
    end subroutine srealloc

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine drealloc &
    (copy, x, n, ierr, factor)
    use maria_kinds_mod,    only:   &
        WP => DP
    implicit none
        logical,    intent(in   )               :: copy
        real(DP),   intent(inout),  allocatable :: x(:)
        integer,    intent(in   )               :: n
        integer,    intent(  out)               :: ierr
        integer,    intent(in   ),  optional    :: factor
    end subroutine drealloc

    !------------------------------------------------------------------------------------------------------------------------
end interface

end module maria_utils_mod
