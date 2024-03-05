program sgemaxvol_rect_swap_cols_test
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
    scopy,                        &
    sgepiv,                       &
    slaset
use maria_lr_maxvol_mod,    only: &
    sgevolume,                    &
    sgemaxvol_rect_swap_cols, sgemaxvol_swap_rows
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
        
        call sgemaxvol_rect_swap_cols(-1, -1, foo, 0, -1, ifoo, ONE-EPS, foo, 0, ifoo, 0, ierr)
        ASSERT(ierr == -1)

        call sgemaxvol_rect_swap_cols(1, -1, foo, 0, -1, ifoo, ONE-EPS, foo, 0, ifoo, 0, ierr)
        ASSERT(ierr == -2)

        call sgemaxvol_rect_swap_cols(1, 1, foo, 0, -1, ifoo, ONE-EPS, foo, 0, ifoo, 0, ierr)
        ASSERT(ierr == -4)

        call sgemaxvol_rect_swap_cols(1, 1, foo, 1, -1, ifoo, ONE-EPS, foo, 0, ifoo, 0, ierr)
        ASSERT(ierr == -5)

        call sgemaxvol_rect_swap_cols(1, 1, foo, 1, 2, ifoo, ONE-EPS, foo, 0, ifoo, 0, ierr)
        ASSERT(ierr == -5)

        call sgemaxvol_rect_swap_cols(1, 1, foo, 1, 1, ifoo, ONE-EPS, foo, 0, ifoo, 0, ierr)
        ASSERT(ierr == -7)
    end subroutine test_sanity

    subroutine test_quick()
        integer     :: ierr, ifoo(1), icol(2), i
        real(WP)    :: foo(1)

        call sgemaxvol_rect_swap_cols(0, 1, foo, 1, 0, ifoo, ONE, foo, 0, ifoo, 0, ierr)
        ASSERT(ierr == 0)

        call sgemaxvol_rect_swap_cols(1, 0, foo, 1, 0, ifoo, ONE, foo, 0, ifoo, 0, ierr)
        ASSERT(ierr == 0)

        call sgemaxvol_rect_swap_cols(1, 1, foo, 1, 0, ifoo, ONE, foo, 0, ifoo, 0, ierr)
        ASSERT(ierr == 0)

        call sgemaxvol_rect_swap_cols(3, 2, foo, 3, 2, icol, ONE, foo, 0, ifoo, 0, ierr)
        ASSERT(ierr == 0)
        do i = 1, 2
            ASSERT(icol(i) == i)
        end do
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer     :: m, n, lda, k, lwork, liwork, ierr, ifoo(1), i
        real(WP)    :: vol0, vol1, foo(1), growth, thresh
        type(prng)  :: rng

        integer,  allocatable :: iwork(:), icol(:)
        real(WP), allocatable :: A(:), work(:), tmp(:)
        
        call rng%init(seed, ierr)
        ASSERT(ierr == 0)

    !-- check
        m = 50
        n = 20
        ldA = m + 1
        k = 10
        thresh = 1.02_WP

        call sgemaxvol_rect_swap_cols(m, n, foo, lda, k, ifoo, thresh, foo, -1, ifoo, -1, ierr)
        ASSERT(ierr == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)
        vol0 = sgevolume(m, k, k, foo, lda, foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
       
        allocate(A(lda*n), icol(k), work(lwork), iwork(liwork), tmp(m*n))

        call rng%snormal(lda*n, A, ZERO, ONE, ierr)
        call arange(k, icol, 1, 1, ierr)
        call slacpy('a', m, n, A, lda, tmp, m)
        vol0 = sgevolume(m, k, k, A, lda, work, lwork, iwork, liwork, ierr)

        call slacpy('a', m, n, tmp, m, A, lda)
        call sgemaxvol_rect_swap_cols(m, n, A, lda, k, icol, thresh, work, lwork, iwork, liwork, ierr)
        growth = work(1)
        call slacpy('a', m, n, tmp, m, A, lda)
        call sgepiv('c', 'f', m, n, A, lda, 1, k, icol, ierr)
        vol1 = sgevolume(m, k, k, A, lda, work, lwork, iwork, liwork, ierr)
        ASSERT(safe_less(ONE, growth))
        ASSERT(are_close(vol1, vol0 * growth, ierr, rtol=1e-5_WP))
        deallocate(A, icol, work, iwork, tmp)
 
    !-- check
        m = 50
        n = 100
        ldA = m + 1
        k = 10
        thresh = 1.02_WP

        call sgemaxvol_rect_swap_cols(m, n, foo, lda, k, ifoo, thresh, foo, -1, ifoo, -1, ierr)
        ASSERT(ierr == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)
        vol0 = sgevolume(m, k, k, foo, lda, foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
       
        allocate(A(lda*n), icol(k), work(lwork), iwork(liwork), tmp(m*n))

        call rng%snormal(lda*n, A, ZERO, ONE, ierr)
        call arange(k, icol, 1, 1, ierr)
        call slacpy('a', m, n, A, lda, tmp, m)
        vol0 = sgevolume(m, k, k, A, lda, work, lwork, iwork, liwork, ierr)

        call slacpy('a', m, n, tmp, m, A, lda)
        call sgemaxvol_rect_swap_cols(m, n, A, lda, k, icol, thresh, work, lwork, iwork, liwork, ierr)
        growth = work(1)
        call slacpy('a', m, n, tmp, m, A, lda)
        call sgepiv('c', 'f', m, n, A, lda, 1, k, icol, ierr)
        vol1 = sgevolume(m, k, k, A, lda, work, lwork, iwork, liwork, ierr)
        ASSERT(safe_less(ONE, growth))
        ASSERT(are_close(vol1, vol0 * growth, ierr, rtol=1e-5_WP))
        deallocate(A, icol, work, iwork, tmp)
 
    !-- check
        m = 50
        n = 100
        ldA = m + 1
        k = 1
        thresh = 1.02_WP

        call sgemaxvol_rect_swap_cols(m, n, foo, lda, k, ifoo, thresh, foo, -1, ifoo, -1, ierr)
        ASSERT(ierr == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)
        vol0 = sgevolume(m, k, k, foo, lda, foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
       
        allocate(A(lda*n), icol(k), work(lwork), iwork(liwork), tmp(m*n))

        call rng%snormal(lda*n, A, ZERO, ONE, ierr)
        call arange(k, icol, 1, 1, ierr)
        call slacpy('a', m, n, A, lda, tmp, m)
        vol0 = sgevolume(m, k, k, A, lda, work, lwork, iwork, liwork, ierr)

        call slacpy('a', m, n, tmp, m, A, lda)
        call sgemaxvol_rect_swap_cols(m, n, A, lda, k, icol, thresh, work, lwork, iwork, liwork, ierr)
        growth = work(1)
        call slacpy('a', m, n, tmp, m, A, lda)
        call sgepiv('c', 'f', m, n, A, lda, 1, k, icol, ierr)
        vol1 = sgevolume(m, k, k, A, lda, work, lwork, iwork, liwork, ierr)
        ASSERT(safe_less(ONE, growth))
        ASSERT(are_close(vol1, vol0 * growth, ierr, rtol=1e-5_WP))
        deallocate(A, icol, work, iwork, tmp)
 
    !-- check
        m = 10
        n = 100
        ldA = m + 1
        k = 10
        thresh = 1.02_WP

        call sgemaxvol_rect_swap_cols(m, n, foo, lda, k, ifoo, thresh, foo, -1, ifoo, -1, ierr)
        ASSERT(ierr == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)
        call sgemaxvol_swap_rows('n', n, m, foo, n, ifoo, thresh, foo, -1, ifoo, -1, ierr)
        ASSERT(ierr == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
        vol0 = sgevolume(m, k, k, foo, lda, foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
       
        allocate(A(lda*n), icol(k), work(lwork), iwork(liwork), tmp(m*n))

        call rng%snormal(lda*n, A, ZERO, ONE, ierr)
        call arange(k, icol, 1, 1, ierr)
        call slacpy('a', m, n, A, lda, tmp, m)

        call slacpy('a', m, n, tmp, m, A, lda)
        call sgemaxvol_rect_swap_cols(m, n, A, lda, k, icol, thresh, work, lwork, iwork, liwork, ierr)
        call slacpy('a', m, n, tmp, m, A, lda)
        call sgepiv('c', 'f', m, n, A, lda, 1, k, icol, ierr)
        vol0 = sgevolume(m, k, k, A, lda, work, lwork, iwork, liwork, ierr)

        call arange(k, icol, 1, 1, ierr)
        do i = 1, n
            call scopy(m, tmp(1 + (i-1)*m:), 1, A(i:), n)
        end do
        call sgemaxvol_swap_rows('n', n, m, A, n, icol, thresh, work, lwork, iwork, liwork, ierr)
        call slacpy('a', m, n, tmp, m, A, lda)
        call sgepiv('c', 'f', m, n, A, lda, 1, k, icol, ierr)
        vol1 = sgevolume(m, k, k, A, lda, work, lwork, iwork, liwork, ierr)
      
        ASSERT(are_close(vol0, vol1, ierr, rtol=1e-5_WP))
        deallocate(A, icol, work, iwork, tmp)
    
    !-- runtime error
        m = 50
        n = 20
        ldA = m + 1
        k = 10
        thresh = 1.02_WP

        call sgemaxvol_rect_swap_cols(m, n, foo, lda, k, ifoo, thresh, foo, -1, ifoo, -1, ierr)
        ASSERT(ierr == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)

        allocate(A(lda*n), icol(k), work(lwork), iwork(liwork), tmp(m*n))

        call rng%snormal(lda*n, A, ZERO, ONE, ierr)
        call slaset('a', m, k, ZERO, ZERO, A, lda)
        call arange(k, icol, 1, 1, ierr)

        call sgemaxvol_rect_swap_cols(m, n, A, lda, k, icol, thresh, work, lwork, iwork, liwork, ierr)
        ASSERT(ierr > 0)

        deallocate(A, icol, work, iwork, tmp)
    end subroutine test
end program sgemaxvol_rect_swap_cols_test
