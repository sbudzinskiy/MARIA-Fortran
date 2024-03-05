program sgemaxvol_test
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
use maria_access_matrix_mod, only: &
    MS => smatslc
use maria_comparison_mod,   only: &
    safe_leq,                     &
    safe_less
use maria_utils_mod,        only: &
    arange,                       &
    permute
use maria_la_utils_mod,     only: &
    all_close
use maria_la_core_mod,      only: &
    slacpy,                       &
    slaset,                       &
    sgemm,                        &
    sgenrmc
use maria_lr_cross_mod,     only: &
    smatcross
use maria_lr_maxvol_mod,   only: &
    sgemaxvol
implicit none (type, external)

    character(100) :: seed_str
    integer        :: seed
    real(WP), allocatable :: A(:)

    call get_command_argument(1, seed_str)
    read(seed_str, *) seed

    call test_sanity()
    call test_quick()
    call test(seed)

contains
    subroutine test_sanity()
        integer                :: info, ifoo(1)
        real(WP)               :: foo(1)
        procedure(MS), pointer :: fun

        call sgemaxvol(-1, -1, fun, -1, ifoo, ifoo, -1, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -1)

        call sgemaxvol(1, -1, fun, -1, ifoo, ifoo, -1, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -2)

        call sgemaxvol(1, 1, fun, -1, ifoo, ifoo, -1, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -4)

        call sgemaxvol(1, 1, fun, 1, ifoo, ifoo, -1, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -7)

        call sgemaxvol(1, 1, fun, 1, ifoo, ifoo, 1, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -8)
    end subroutine test_sanity

    subroutine test_quick()
        integer                :: info, ifoo(1)
        real(WP)               :: foo(1)
        procedure(MS), pointer :: fun

        call sgemaxvol(0, 1, fun, 0, ifoo, ifoo, 1, ONE, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)

        call sgemaxvol(1, 0, fun, 0, ifoo, ifoo, 1, ONE, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer                :: m, n, r, niter, ldu, ldvt, lwork, liwork, info, ifoo(1)
        real(WP)               :: thresh, foo(1), err0, err1
        type(prng)             :: rng
        procedure(MS), pointer :: fun
        integer,  allocatable  :: iwork(:), icol(:), irow(:)
        real(WP), allocatable  :: tmp(:), U(:), VT(:), work(:)
        
        call rng%init(seed, info)
        ASSERT(info == 0)

        fun => matslc

    !-- check
        m = 15
        n = 20
        r = 5
        niter = 10
        thresh = 1.02_WP
        ldu = m + 1
        ldvt = r + 1

        call sgemaxvol(m, n, fun, r, ifoo, ifoo, niter, thresh, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)
        call smatcross(m, n, fun, r, ifoo, r, ifoo, r, foo, ldu, foo, ldvt, foo, -1, ifoo, -1, info)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        allocate(A(m*n), tmp(m*n), icol(n), irow(m), U(ldu*r), VT(ldvt*n), work(lwork), iwork(liwork))
        
        call rng%snormal(m*n, A, ZERO, ONE, info)
        call arange(n, icol, 1, 1, info)
        call arange(m, irow, 1, 1, info)
        call permute(rng, n, icol, info)
        call permute(rng, m, irow, info)
         
        call smatcross(m, n, fun, r, icol, r, irow, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
        call slacpy('a', m, n, A, m, tmp, m)
        call sgemm('n', 'n', m, n, r, -ONE, U, ldu, VT, ldvt, ONE, tmp, m)
        err0 = sgenrmc(m, n, tmp, m, info)
        
        call sgemaxvol(m, n, fun, r, irow, icol, niter, thresh, work, lwork, iwork, liwork, info)
        ASSERT(info == 0)
        call smatcross(m, n, fun, r, icol, r, irow, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
        call slacpy('a', m, n, A, m, tmp, m)
        call sgemm('n', 'n', m, n, r, -ONE, U, ldu, VT, ldvt, ONE, tmp, m)
        err1 = sgenrmc(m, n, tmp, m, info)
        ASSERT(safe_less(err1, err0))

        deallocate(A, tmp, icol, irow, U, VT, work, iwork)
    end subroutine test

    real(WP) function matval &
    (m, n, i, j, info)
        integer, intent(in) :: m
        integer, intent(in) :: n
        integer, intent(in) :: i
        integer, intent(in) :: j
        integer, intent(out) :: info

        matval = A(i + (j-1)*m)
        info = 0 * n
    end function matval

    subroutine matslc &
    (m, n, mode, ind, x, incx, info)
    use maria_access_matrix_mod, only: &
        smatval, smatval2slc
        integer, intent(in) :: m
        integer, intent(in) :: n
        integer, intent(in) :: mode
        integer, intent(in) :: ind
        real(WP), intent(out), contiguous :: x(:)
        integer, intent(in) :: incx
        integer, intent(out) :: info
        
        procedure(smatval), pointer :: fun

        fun => matval
        
        call smatval2slc(fun, m, n, mode, ind, x, incx, info)
    end subroutine
end program sgemaxvol_test
