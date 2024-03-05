program suniform_test
#include "maria_assert_mod.h"
#ifdef MARIA_MKL
use maria_prng_mkl_mod,     only: &
    prng => prng_mkl
#else
use maria_prng_builtin_mod, only: &
    prng => prng_builtin
#endif
use maria_kinds_mod,        only: &
    WP => SP
use maria_constants_mod,    only: &
    ZERO => S_ZERO,               &
    HALF => S_HALF,               &
    ONE => S_ONE,                 &
    EPS => S_MACHTOL
use maria_comparison_mod,   only: &
    safe_eq,                      &
    safe_less
implicit none (type, external)

    character(100) :: seed_str
    integer        :: seed

    call get_command_argument(1, seed_str)
    read(seed_str, *) seed

    call test_sanity(seed)
    call test_quick(seed)
    call test(seed)

contains
    subroutine test_sanity(seed)
        integer, intent(in) :: seed

        integer    :: info
        type(prng) :: rng
        real(WP)   :: foo(1)

        call rng%init(seed, info)
        ASSERT(info == 0)

        call rng%suniform(-1, foo, ZERO, -EPS, info)
        ASSERT(info == -2)
 
        call rng%suniform(10, foo, ZERO, -EPS, info)
        ASSERT(info == -5)

        call rng%suniform(10, foo, ZERO, ZERO, info)
        ASSERT(info == -5)

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test_sanity

    subroutine test_quick(seed)
        integer, intent(in) :: seed

        integer    :: info
        real(WP)   :: foo(1)
        type(prng) :: rng

        call rng%init(seed, info)
        ASSERT(info == 0)

        call rng%suniform(0, foo, ZERO, ONE, info)
        ASSERT(info == 0)

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer               :: info, n
        real(WP)              :: a, b, emax, emin, prob, tol
        real(WP), allocatable :: x(:)
        type(prng)            :: rng

        call rng%init(seed, info)
        ASSERT(info == 0)

        n = 10000
        a = -HALF
        b = 3*HALF

        allocate(x(n + 1))
        x = b + ONE

        call rng%suniform(n, x, a, b, info)
        ASSERT(info == 0)

    !-- Check range ------------------------------------------------------------
        ASSERT(safe_eq(x(n + 1), b + ONE))

    !-- Estimate max -----------------------------------------------------------
        emax = maxval(x(1:n))
        ! Pr( emax <= b - tol ) = (b - a - tol)**n / (b - a)**n
        prob = 1e-2_WP
        tol = (b - a) * (ONE - prob**(ONE/n))
        ASSERT(safe_less(b - tol, emax))

    !-- Estimate min -----------------------------------------------------------
        emin = minval(x(1:n))
        ! Pr( emin >= a + tol ) = (b - a - tol)**n / (b - a)**n
        ASSERT(safe_less(emin, a + tol))

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test   
end program suniform_test
