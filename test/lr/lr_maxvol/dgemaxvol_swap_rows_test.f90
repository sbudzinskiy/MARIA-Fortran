program dgemaxvol_swap_rows_test
#include "maria_assert_mod.h"
#ifdef MARIA_MKL
use maria_prng_mkl_mod,     only: &
    prng => prng_mkl
#else
use maria_prng_builtin_mod, only: &
    prng => prng_builtin
#endif
use maria_kinds_mod,        only: &
    WP => DP
use maria_constants_mod,    only: &
    ONE => D_ONE,                 &
    ZERO => D_ZERO,               &
    EPS => D_MACHTOL
use maria_comparison_mod,   only: &
    safe_eq,                      &
    safe_less
use maria_utils_mod,        only: &
    are_close,                    &
    arange
use maria_la_core_mod,      only: &
    dlasrt,                       &
    dlacpy,                       &
    dgepiv,                       &
    dlaset
use maria_lr_maxvol_mod,    only: &
    dgevolume,                    &
    dgemaxvol_swap_rows
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
        integer                :: info, ifoo(1)
        real(WP)               :: foo(1)
        
        call dgemaxvol_swap_rows('a', -1, -1, foo, 0, ifoo, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -1)

        call dgemaxvol_swap_rows('y', -1, -1, foo, 0, ifoo, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -2)

        call dgemaxvol_swap_rows('y', 1, -1, foo, 0, ifoo, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -3)

        call dgemaxvol_swap_rows('y', 1, 2, foo, 0, ifoo, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -3)

        call dgemaxvol_swap_rows('y', 1, 1, foo, 0, ifoo, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -5)

        call dgemaxvol_swap_rows('y', 1, 1, foo, 1, ifoo, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -7)
    end subroutine test_sanity

    subroutine test_quick()
        integer                :: info, ifoo(1)
        real(WP)               :: foo(1)

        call dgemaxvol_swap_rows('y', 1, 0, foo, 1, ifoo, ONE+EPS, foo, 1, ifoo, 1, info)
        ASSERT(info == 0)

        call dgemaxvol_swap_rows('y', 1, 1, foo, 1, ifoo, ONE+EPS, foo, 1, ifoo, 1, info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer                :: m, n, lda, lwork, liwork, info, ifoo(1)
        real(WP)               :: vol0, vol1, foo(1), thresh, growth
        type(prng)             :: rng
        integer,  allocatable  :: iwork(:), irow(:)
        real(WP), allocatable  :: A(:), work(:), tmp(:)
        
        call rng%init(seed, info)
        ASSERT(info == 0)

    !-- check
        m = 50
        n = 10
        ldA = m + 1
        thresh = 1.02_WP

        call dgemaxvol_swap_rows('n', m, n, foo, lda, ifoo, thresh, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)
        call dgemaxvol_swap_rows('y', m, n, foo, lda, ifoo, thresh, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        vol0 = dgevolume(n, n, n, foo, lda, foo, -1, ifoo, -1, info)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
        
        allocate(A(lda*n), irow(n), work(lwork), iwork(liwork), tmp(m*n))

        call rng%dnormal(lda*n, A, ZERO, ONE, info)
        call arange(n, irow, 1, 1, info)
        call dlacpy('a', m, n, A, lda, tmp, m)
        vol0 = dgevolume(n, n, n, A, lda, work, lwork, iwork, liwork, info)

        call dlacpy('a', m, n, tmp, m, A, lda)
        call dgemaxvol_swap_rows('n', m, n, A, lda, irow, thresh, work, lwork, iwork, liwork, info)
        growth = work(1)
        call dlacpy('a', m, n, tmp, m, A, lda)
        call dgepiv('r', 'f', m, n, A, lda, 1, n, irow, info)
        vol1 = dgevolume(n, n, n, A, lda, work, lwork, iwork, liwork, info)

        ASSERT(safe_less(ONE, growth))
        ASSERT(are_close(vol1, vol0 * growth, info, rtol=1e-12_WP))

        call dlacpy('a', m, n, tmp, m, A, lda)
        call arange(n, irow, 1, 1, info)
        call dgemaxvol_swap_rows('y', m, n, A, lda, irow, thresh, work, lwork, iwork, liwork, info)
        call dlacpy('a', m, n, tmp, m, A, lda)
        call dgepiv('r', 'f', m, n, A, lda, 1, n, irow, info)
        vol0 = dgevolume(n, n, n, A, lda, work, lwork, iwork, liwork, info)
        ASSERT(are_close(vol1, vol0, info, rtol=1e-12_WP))

        deallocate(A, irow, work, iwork, tmp)

    !-- runtime error
        m = 50
        n = 10
        ldA = m + 1
        thresh = 1.02_WP

        call dgemaxvol_swap_rows('n', m, n, foo, lda, ifoo, thresh, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)
        call dgemaxvol_swap_rows('y', m, n, foo, lda, ifoo, thresh, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        allocate(A(lda*n), irow(n), work(lwork), iwork(liwork), tmp(m*n))

        call rng%dnormal(lda*n, A, ZERO, ONE, info)
        call dlaset('a', n, n, ZERO, ZERO, A, lda)
        call arange(n, irow, 1, 1, info)

        call dlacpy('a', m, n, A, lda, tmp, m)
        call dgemaxvol_swap_rows('n', m, n, A, lda, irow, thresh, work, lwork, iwork, liwork, info)
        ASSERT(info > 0)

        call dlacpy('a', m, n, A, lda, tmp, m)
        call dgemaxvol_swap_rows('y', m, n, A, lda, irow, thresh, work, lwork, iwork, liwork, info)
        ASSERT(info > 0)

        deallocate(A, irow, work, iwork, tmp)
    end subroutine test
end program dgemaxvol_swap_rows_test
