program sgemaxvol_rect_add_rows_test
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
    sgemaxvol_rect_add_rows
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
        integer     :: ierr, ifoo(1)
        real(WP)    :: foo(1)
        
        call sgemaxvol_rect_add_rows('a', -1, -1, foo, 0, -1, ifoo, foo, 0, ifoo, 0, ierr)
        ASSERT(ierr == -1)

        call sgemaxvol_rect_add_rows('y', -1, -1, foo, 0, -1, ifoo, foo, 0, ifoo, 0, ierr)
        ASSERT(ierr == -2)

        call sgemaxvol_rect_add_rows('y', 1, -1, foo, 0, -1, ifoo, foo, 0, ifoo, 0, ierr)
        ASSERT(ierr == -3)

        call sgemaxvol_rect_add_rows('y', 1, 2, foo, 0, -1, ifoo, foo, 0, ifoo, 0, ierr)
        ASSERT(ierr == -3)

        call sgemaxvol_rect_add_rows('y', 1, 1, foo, 0, -1, ifoo, foo, 0, ifoo, 0, ierr)
        ASSERT(ierr == -5)

        call sgemaxvol_rect_add_rows('y', 3, 2, foo, 3, 1, ifoo, foo, 0, ifoo, 0, ierr)
        ASSERT(ierr == -6)

        call sgemaxvol_rect_add_rows('y', 3, 2, foo, 3, 4, ifoo, foo, 0, ifoo, 0, ierr)
        ASSERT(ierr == -6)
    end subroutine test_sanity

    subroutine test_quick()
        integer     :: ierr, ifoo(1), irow(3), i
        real(WP)    :: foo(1)

        call sgemaxvol_rect_add_rows('y', 0, 0, foo, 1, 0, ifoo, foo, 0, ifoo, 0, ierr)
        ASSERT(ierr == 0)

        call sgemaxvol_rect_add_rows('y', 1, 0, foo, 1, 0, ifoo, foo, 0, ifoo, 0, ierr)
        ASSERT(ierr == 0)

        call sgemaxvol_rect_add_rows('y', 3, 2, foo, 3, 3, irow, foo, 0, ifoo, 0, ierr)
        ASSERT(ierr == 0)
        do i = 1, 3
            ASSERT(irow(i) == i)
        end do

        irow(1) = 3
        irow(2) = 2
        call sgemaxvol_rect_add_rows('y', 3, 2, foo, 3, 2, irow, foo, 0, ifoo, 0, ierr)
        ASSERT(ierr == 0)
        ASSERT(irow(1) == 2)
        ASSERT(irow(2) == 3)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer     :: m, n, lda, k, lwork, liwork, ierr, ifoo(1)
        real(WP)    :: vol0, vol1, foo(1), growth
        type(prng)  :: rng

        integer,  allocatable :: iwork(:), irow(:)
        real(WP), allocatable :: A(:), work(:), tmp(:)
        
        call rng%init(seed, ierr)
        ASSERT(ierr == 0)

    !-- check
        m = 50
        n = 10
        ldA = m + 1
        k = 15

        call sgemaxvol_rect_add_rows('n', m, n, foo, lda, k, ifoo, foo, -1, ifoo, -1, ierr)
        ASSERT(ierr == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)
        call sgemaxvol_rect_add_rows('y', m, n, foo, lda, k, ifoo, foo, -1, ifoo, -1, ierr)
        ASSERT(ierr == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
        vol0 = sgevolume(n, n, n, foo, lda, foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
        vol0 = sgevolume(k, n, n, foo, lda, foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
       
        allocate(A(lda*n), irow(k), work(lwork), iwork(liwork), tmp(m*n))

        call rng%snormal(lda*n, A, ZERO, ONE, ierr)
        call arange(n, irow, 1, 1, ierr)
        call slacpy('a', m, n, A, lda, tmp, m)
        vol0 = sgevolume(n, n, n, A, lda, work, lwork, iwork, liwork, ierr)

        call slacpy('a', m, n, tmp, m, A, lda)
        call sgemaxvol_rect_add_rows('n', m, n, A, lda, k, irow, work, lwork, iwork, liwork, ierr)
        growth = work(1)
        call slacpy('a', m, n, tmp, m, A, lda)
        call sgepiv('r', 'f', m, n, A, lda, 1, k, irow, ierr)
        vol1 = sgevolume(k, n, n, A, lda, work, lwork, iwork, liwork, ierr)

        ASSERT(safe_less(ONE, growth))
        ASSERT(are_close(vol1, vol0 * growth, ierr, rtol=1e-5_WP))

        call slacpy('a', m, n, tmp, m, A, lda)
        call arange(n, irow, 1, 1, ierr)
        call sgemaxvol_rect_add_rows('y', m, n, A, lda, k, irow, work, lwork, iwork, liwork, ierr)
        call slacpy('a', m, n, tmp, m, A, lda)
        call sgepiv('r', 'f', m, n, A, lda, 1, k, irow, ierr)
        vol0 = sgevolume(k, n, n, A, lda, work, lwork, iwork, liwork, ierr)
        ASSERT(are_close(vol1, vol0, ierr, rtol=1e-5_WP))

        deallocate(A, irow, work, iwork, tmp)

    !-- runtime error
        m = 50
        n = 10
        ldA = m + 1
        k = 15

        call sgemaxvol_rect_add_rows('n', m, n, foo, lda, k, ifoo, foo, -1, ifoo, -1, ierr)
        ASSERT(ierr == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)
        call sgemaxvol_rect_add_rows('y', m, n, foo, lda, k, ifoo, foo, -1, ifoo, -1, ierr)
        ASSERT(ierr == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        allocate(A(lda*n), irow(k), work(lwork), iwork(liwork), tmp(m*n))

        call rng%snormal(lda*n, A, ZERO, ONE, ierr)
        call slaset('a', n, n, ZERO, ZERO, A, lda)
        call arange(n, irow, 1, 1, ierr)

        call slacpy('a', m, n, A, lda, tmp, m)
        call sgemaxvol_rect_add_rows('n', m, n, A, lda, k, irow, work, lwork, iwork, liwork, ierr)
        ASSERT(ierr > 0)

        call slacpy('a', m, n, A, lda, tmp, m)
        call sgemaxvol_rect_add_rows('y', m, n, A, lda, k, irow, work, lwork, iwork, liwork, ierr)
        ASSERT(ierr > 0)

        deallocate(A, irow, work, iwork, tmp)
    end subroutine test
end program sgemaxvol_rect_add_rows_test
