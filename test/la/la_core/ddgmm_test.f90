program ddgmm_test
#include "maria_assert_mod.h"
#ifdef MARIA_MKL
use maria_prng_mkl_mod,     only: &
    prng => prng_mkl
#else
use maria_prng_builtin_mod, only: &
    prng => prng_builtin
#endif
use maria_la_core_mod,      only: &
    ddgmm,                        &
    dscal,                        &
    dlacpy,                       &
    dgescal
use maria_kinds_mod,        only: &
    WP => DP
use maria_constants_mod,    only: &
    ZERO => D_ZERO,               &
    ONE => D_ONE
use maria_la_utils_mod,     only: &
    all_close,                    &
    geall_close
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
        real(WP) :: foo(1)
        
        call ddgmm('x', -1, -1, foo, 0, foo, 0, info)
        ASSERT(info == -1)

        call ddgmm('l', -1, -1, foo, 0, foo, 0, info)
        ASSERT(info == -2)

        call ddgmm('l', 1, -1, foo, 0, foo, 0, info)
        ASSERT(info == -3)

        call ddgmm('l', 1, 1, foo, 0, foo, 0, info)
        ASSERT(info == -5)
    end subroutine test_sanity

    subroutine test_quick()
        integer  :: info
        real(WP) :: foo(1)
        
        call ddgmm('l', 0, 1, foo, 1, foo, 0, info)
        ASSERT(info == 0)

        call ddgmm('l', 1, 0, foo, 1, foo, 0, info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        character(1)          :: side
        integer               :: info, m, n, ldA, incD, lenD, i
        type(prng)            :: rng
        real(WP), allocatable :: A(:), D(:), tmp(:)

        call rng%init(seed, info)
        ASSERT(info == 0)

    !-- side = 'l', incD > 0 ---------------------------------------------------
        side = 'l'
        m = 10
        n = 20
        lda = m + 1
        incD = 2
        lenD = 1 + (m-1) * abs(incD)
        allocate(A(ldA*n), D(lenD), tmp(m*n))
        call rng%dnormal(lda*n, A, ZERO, ONE, info)
        call dlacpy('a', m, n, A, lda, tmp, m)
        call rng%dnormal(lend, D, ZERO, ONE, info)
        call ddgmm(side, m, n, A, ldA, D, incD, info)
        ASSERT(info == 0)
        do i = 1, m
            call dscal(n, D(1 + (i-1)*incD), tmp(i:), m)
            ASSERT(all_close(n, tmp(i:), m, A(i:), ldA, info, rtol=1e-5_WP))
        end do
        deallocate(A, D, tmp)

    !-- side = 'l', incD < 0 ---------------------------------------------------
        side = 'l'
        m = 10
        n = 20
        lda = m + 1
        incD = -3
        lenD = 1 + (m-1) * abs(incD)
        allocate(A(ldA*n), D(lenD), tmp(m*n))
        call rng%dnormal(lda*n, A, ZERO, ONE, info)
        call dlacpy('a', m, n, A, lda, tmp, m)
        call rng%dnormal(lend, D, ZERO, ONE, info)
        call ddgmm(side, m, n, A, ldA, D, incD, info)
        ASSERT(info == 0)
        do i = 1, m
            call dscal(n, D(lend + (i-1)*incD), tmp(i:), m)
            ASSERT(all_close(n, tmp(i:), m, A(i:), ldA, info, rtol=1e-5_WP))
        end do
        deallocate(A, D, tmp)

    !-- side = 'r', incD > 0 ---------------------------------------------------
        side = 'r'
        m = 10
        n = 20
        lda = m + 1
        incD = 2
        lenD = 1 + (n-1) * abs(incD)
        allocate(A(ldA*n), D(lenD), tmp(m*n))
        call rng%dnormal(lda*n, A, ZERO, ONE, info)
        call dlacpy('a', m, n, A, lda, tmp, m)
        call rng%dnormal(lend, D, ZERO, ONE, info)
        call ddgmm(side, m, n, A, ldA, D, incD, info)
        ASSERT(info == 0)
        do i = 1, n
            call dscal(m, D(1 + (i-1)*incD), tmp(1 + (i-1)*m:), 1)
            ASSERT(all_close(m, tmp(1 + (i-1)*m:), 1, A(1 + (i-1)*ldA:), 1, info, rtol=1e-5_WP))
        end do
        deallocate(A, D, tmp)

    !-- side = 'r', incD < 0 ---------------------------------------------------
        side = 'r'
        m = 10
        n = 20
        lda = m + 1
        incD = -3
        lenD = 1 + (n-1) * abs(incD)
        allocate(A(ldA*n), D(lenD), tmp(m*n))
        call rng%dnormal(lda*n, A, ZERO, ONE, info)
        call dlacpy('a', m, n, A, lda, tmp, m)
        call rng%dnormal(lend, D, ZERO, ONE, info)
        call ddgmm(side, m, n, A, ldA, D, incD, info)
        ASSERT(info == 0)
        do i = 1, n
            call dscal(m, D(lenD + (i-1)*incD), tmp(1 + (i-1)*m:), 1)
            ASSERT(all_close(m, tmp(1 + (i-1)*m:), 1, A(1 + (i-1)*ldA:), 1, info, rtol=1e-5_WP))
        end do
        deallocate(A, D, tmp)

    !-- side = 'l', incD = 0 ---------------------------------------------------
        side = 'l'
        m = 10
        n = 20
        lda = m + 1
        incD = 0
        lenD = 1
        allocate(A(ldA*n), D(lenD), tmp(m*n))
        call rng%dnormal(lda*n, A, ZERO, ONE, info)
        call dlacpy('a', m, n, A, lda, tmp, m)
        call rng%dnormal(lend, D, ZERO, ONE, info)
        call ddgmm(side, m, n, A, ldA, D, incD, info)
        ASSERT(info == 0)
        call dgescal(m, n, D(1), tmp, m, info)
        ASSERT(geall_close(m, n, A, ldA, tmp, m, info, rtol=1e-5_WP))

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test
end program ddgmm_test
