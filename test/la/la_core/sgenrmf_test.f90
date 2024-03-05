program sgenrmf_test
#include "maria_assert_mod.h"
#ifdef MARIA_MKL
use maria_prng_mkl_mod,     only: &
    prng => prng_mkl
#else
use maria_prng_builtin_mod, only: &
    prng => prng_builtin
#endif
use maria_la_core_mod,      only: &
    sgenrmf,                      &
    sgedotf,                      &
    slaset
use maria_kinds_mod,        only: &
    WP => SP
use maria_constants_mod,    only: &
    ZERO => S_ZERO,               &
    ONE => S_ONE
use maria_comparison_mod,   only: &
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
        
        nrm = sgenrmf(-1, -1, foo, 0, info)
        ASSERT(info == -1)

        nrm = sgenrmf(1, -1, foo, 0, info)
        ASSERT(info == -2)

        nrm = sgenrmf(1, 1, foo, 0, info)
        ASSERT(info == -4)
    end subroutine test_sanity

    subroutine test_quick()
        integer  :: info
        real(WP) :: foo(1), nrm
        
        nrm = sgenrmf(0, 1, foo, 1, info)
        ASSERT(info == 0)
        ASSERT(safe_eq(nrm, ZERO))

        nrm = sgenrmf(1, 0, foo, 1, info)
        ASSERT(info == 0)
        ASSERT(safe_eq(nrm, ZERO))
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer               :: info, m, n, ldA
        real(WP)              :: nrm1, nrm2
        type(prng)            :: rng
        real(WP), allocatable :: A(:)

        call rng%init(seed, info)
        ASSERT(info == 0)

    !-- ldA > m ----------------------------------------------------------------
        m = 10
        n = 20
        lda = m + 1
        allocate(A(lda*n))
        call rng%snormal(lda*n, A, ZERO, ONE, info)
        nrm1 = sgenrmf(m, n, A, lda, info)
        ASSERT(info == 0)
        nrm2 = sqrt(sgedotf('n', m, n, A, lda, A, lda, info))
        ASSERT(are_close(nrm1, nrm2, info, rtol=1e-5_WP))
        deallocate(A)

    !-- ldA == m ---------------------------------------------------------------
        m = 10
        n = 20
        lda = m + 1
        allocate(A(lda*n))
        call rng%snormal(lda*n, A, ZERO, ONE, info)
        call slaset('a', lda-m, n, ZERO, ZERO, A(m+1:), ldA)
        nrm1 = sgenrmf(m, n, A, lda, info)
        ASSERT(info == 0)
        nrm2 = sgenrmf(lda, n, A, lda, info)
        ASSERT(info == 0)
        ASSERT(are_close(nrm1, nrm2, info, rtol=1e-5_WP))

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test
end program sgenrmf_test
