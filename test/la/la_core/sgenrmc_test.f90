program sgenrmc_test
#include "maria_assert_mod.h"
#ifdef MARIA_MKL
use maria_prng_mkl_mod,     only: &
    prng => prng_mkl
#else
use maria_prng_builtin_mod, only: &
    prng => prng_builtin
#endif
use maria_la_core_mod,      only: &
    sgenrmc
use maria_kinds_mod,        only: &
    WP => SP
use maria_constants_mod,    only: &
    ZERO => S_ZERO,               &
    ONE => S_ONE
use maria_comparison_mod,   only: &
    safe_leq,                     &
    safe_eq
use maria_utils_mod,        only: &
    are_close
implicit none (type, external)

    character(100) :: seed_str
    integer        :: seed

    call get_command_argument(1, seed_str)
    read(seed_str, *) seed

    call test_sanity()
    call test_quick()
    call test(seed)

contains
    subroutine test_sanity()
        integer  :: info
        real(WP) :: foo(1), nrm
        
        nrm = sgenrmc(-1, -1, foo, 0, info)
        ASSERT(info == -1)

        nrm = sgenrmc(1, -1, foo, 0, info)
        ASSERT(info == -2)

        nrm = sgenrmc(1, 1, foo, 0, info)
        ASSERT(info == -4)
    end subroutine test_sanity

    subroutine test_quick()
        integer  :: info, pos(2)
        real(WP) :: foo(1), nrm
        
        nrm = sgenrmc(0, 1, foo, 1, info, pos=pos)
        ASSERT(info == 0)
        ASSERT(safe_eq(nrm, ZERO))
        ASSERT(pos(1) == 1)
        ASSERT(pos(2) == 1)

        nrm = sgenrmc(1, 0, foo, 1, info, pos=pos)
        ASSERT(info == 0)
        ASSERT(safe_eq(nrm, ZERO))
        ASSERT(pos(1) == 1)
        ASSERT(pos(2) == 1)
    end subroutine test_quick

    subroutine test(seed)
    use, intrinsic :: ieee_arithmetic, only: &
        ieee_value, IEEE_QUIET_NAN
        integer, intent(in) :: seed

        integer               :: info, m, n, ldA, i, j, pos(2)
        real(WP)              :: nrm, nan
        type(prng)            :: rng
        real(WP), allocatable :: A(:)

        call rng%init(seed, info)
        ASSERT(info == 0)

    !---------------------------------------------------------------------------
        m = 10
        n = 20
        lda = m + 1
        allocate(A(lda*n))
        call rng%suniform(lda*n, A, -ONE, ONE, info)
        nrm = sgenrmc(m, n, A, lda, info)
        ASSERT(info == 0)
        ASSERT(safe_leq(nrm, ONE))
        i = 3
        j = 7
        A(i + (j-1)*ldA) = 2 * ONE
        A(i + j*ldA) = 2 * ONE
        nrm = sgenrmc(m, n, A, lda, info, pos=pos)
        ASSERT(safe_eq(nrm, abs(A(i + (j-1)*ldA))))
        ASSERT(pos(1) == i)
        ASSERT(pos(2) == j)
        deallocate(A)

    !-- NaN --------------------------------------------------------------------
        m = 10
        n = 20
        lda = m + 1
        allocate(A(lda*n))
        call rng%suniform(lda*n, A, -ONE, ONE, info)
        i = 4
        j = 5
        nan = ieee_value(nan, IEEE_QUIET_NAN)
        A(i + (j-1)*ldA) = nan
        A(i + j*ldA) = nan
        nrm = sgenrmc(m, n, A, lda, info, pos=pos)
        ASSERT(info == 0)
        ASSERT(isnan(nrm))
        ASSERT(pos(1) == i)
        ASSERT(pos(2) == j)

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test
end program sgenrmc_test
