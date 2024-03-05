program sgemaxvol_swap_rows_test
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
    ONE => S_ONE,                 &
    ZERO => S_ZERO,               &
    EPS => S_MACHTOL
use maria_comparison_mod,   only: &
    safe_eq,                      &
    safe_less
use maria_utils_mod,        only: &
    are_close,                    &
    arange
use maria_la_core_mod,      only: &
    slasrt,                       &
    slacpy,                       &
    sgepiv,                       &
    slaset
use maria_lr_maxvol_mod,    only: &
    sgevolume,                    &
    sgemaxvol_swap_rows
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
        
        call sgemaxvol_swap_rows('a', -1, -1, foo, 0, ifoo, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -1)

        call sgemaxvol_swap_rows('y', -1, -1, foo, 0, ifoo, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -2)

        call sgemaxvol_swap_rows('y', 1, -1, foo, 0, ifoo, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -3)

        call sgemaxvol_swap_rows('y', 1, 2, foo, 0, ifoo, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -3)

        call sgemaxvol_swap_rows('y', 1, 1, foo, 0, ifoo, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -5)

        call sgemaxvol_swap_rows('y', 1, 1, foo, 1, ifoo, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -7)
    end subroutine test_sanity

    subroutine test_quick()
        integer                :: info, ifoo(1)
        real(WP)               :: foo(1)

        call sgemaxvol_swap_rows('y', 1, 0, foo, 1, ifoo, ONE+EPS, foo, 1, ifoo, 1, info)
        ASSERT(info == 0)

        call sgemaxvol_swap_rows('y', 1, 1, foo, 1, ifoo, ONE+EPS, foo, 1, ifoo, 1, info)
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

        call sgemaxvol_swap_rows('n', m, n, foo, lda, ifoo, thresh, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)
        call sgemaxvol_swap_rows('y', m, n, foo, lda, ifoo, thresh, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        vol0 = sgevolume(n, n, n, foo, lda, foo, -1, ifoo, -1, info)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
        
        allocate(A(lda*n), irow(n), work(lwork), iwork(liwork), tmp(m*n))

        call rng%snormal(lda*n, A, ZERO, ONE, info)
        call arange(n, irow, 1, 1, info)
        call slacpy('a', m, n, A, lda, tmp, m)
        vol0 = sgevolume(n, n, n, A, lda, work, lwork, iwork, liwork, info)

        call slacpy('a', m, n, tmp, m, A, lda)
        call sgemaxvol_swap_rows('n', m, n, A, lda, irow, thresh, work, lwork, iwork, liwork, info)
        growth = work(1)
        call slacpy('a', m, n, tmp, m, A, lda)
        call sgepiv('r', 'f', m, n, A, lda, 1, n, irow, info)
        vol1 = sgevolume(n, n, n, A, lda, work, lwork, iwork, liwork, info)

        ASSERT(safe_less(ONE, growth))
        ASSERT(are_close(vol1, vol0 * growth, info, rtol=1e-5_WP))

        call slacpy('a', m, n, tmp, m, A, lda)
        call arange(n, irow, 1, 1, info)
        call sgemaxvol_swap_rows('y', m, n, A, lda, irow, thresh, work, lwork, iwork, liwork, info)
        call slacpy('a', m, n, tmp, m, A, lda)
        call sgepiv('r', 'f', m, n, A, lda, 1, n, irow, info)
        vol0 = sgevolume(n, n, n, A, lda, work, lwork, iwork, liwork, info)
        ASSERT(are_close(vol1, vol0, info, rtol=1e-5_WP))

        deallocate(A, irow, work, iwork, tmp)

    !-- runtime error
        m = 50
        n = 10
        ldA = m + 1
        thresh = 1.02_WP

        call sgemaxvol_swap_rows('n', m, n, foo, lda, ifoo, thresh, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)
        call sgemaxvol_swap_rows('y', m, n, foo, lda, ifoo, thresh, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        allocate(A(lda*n), irow(n), work(lwork), iwork(liwork), tmp(m*n))

        call rng%snormal(lda*n, A, ZERO, ONE, info)
        call slaset('a', n, n, ZERO, ZERO, A, lda)
        call arange(n, irow, 1, 1, info)

        call slacpy('a', m, n, A, lda, tmp, m)
        call sgemaxvol_swap_rows('n', m, n, A, lda, irow, thresh, work, lwork, iwork, liwork, info)
        ASSERT(info > 0)

        call slacpy('a', m, n, A, lda, tmp, m)
        call sgemaxvol_swap_rows('y', m, n, A, lda, irow, thresh, work, lwork, iwork, liwork, info)
        ASSERT(info > 0)

        deallocate(A, irow, work, iwork, tmp)
    end subroutine test
end program sgemaxvol_swap_rows_test