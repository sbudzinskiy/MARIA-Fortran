program snormal_test
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
    safe_eq
use maria_utils_mod,        only: &
    are_close
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

        call rng%snormal(-1, foo, ZERO, -EPS, info)
        ASSERT(info == -2)

        call rng%snormal(10, foo, ZERO, -EPS, info)
        ASSERT(info == -5)

        call rng%snormal(10, foo, ZERO, ZERO, info)
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

        call rng%snormal(0, foo, ZERO, ONE, info)
        ASSERT(info == 0)

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer               :: info, n
        real(WP)              :: mean, std, emean, estd, prob, tol
        real(WP), allocatable :: x(:)
        type(prng)            :: rng

        call rng%init(seed, info)
        ASSERT(info == 0)

        n = 1000000
        mean = ONE
        std = HALF

        allocate(x(n + 1))
        x = ZERO

        call rng%snormal(n, x, mean, std, info)
        ASSERT(info == 0)

    !-- Check range ------------------------------------------------------------
        ASSERT(safe_eq(x(n+1), ZERO))

    !-- Estimate mean ----------------------------------------------------------
        emean = sum(x(1:n)) / n
        ! From Chebyshev inequality,
        ! Pr( |emean - mean| >= tol ) <= prob, tol = std / sqrt(prob * n)
        prob = 1e-2_WP
        tol = std / sqrt(prob * n)
        ASSERT(are_close(emean, mean, info, atol=tol))

    !-- Estimate std -----------------------------------------------------------
        estd = sqrt(sum((x(1:n) - emean)**2 / (n-1)))
        ASSERT(are_close(estd, std, info, atol=tol))
    end subroutine test   
end program snormal_test
