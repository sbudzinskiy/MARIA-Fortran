program sgemaxvol_proj_test
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
    MS => smatslc,                 &
    smatval_ballistic
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
    sgenrmc, sgenrmf
use maria_lr_cross_mod,     only: &
    smatcross
use maria_lr_maxvol_mod,   only: &
    sgemaxvol, sgemaxvol_proj
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

        call sgemaxvol_proj(-1, -1, fun, -1, -1, ifoo, ifoo, -1, ifoo, ifoo, -1, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -1)

        call sgemaxvol_proj(1, -1, fun, -1, -1, ifoo, ifoo, -1, ifoo, ifoo, -1, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -2)

        call sgemaxvol_proj(1, 2, fun, -1, -1, ifoo, ifoo, -1, ifoo, ifoo, -1, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -4)

        call sgemaxvol_proj(1, 2, fun, 2, -1, ifoo, ifoo, -1, ifoo, ifoo, -1, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -4)

        call sgemaxvol_proj(1, 2, fun, 1, 0, ifoo, ifoo, -1, ifoo, ifoo, -1, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -5)

        call sgemaxvol_proj(1, 2, fun, 1, 2, ifoo, ifoo, -1, ifoo, ifoo, -1, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -5)

        call sgemaxvol_proj(1, 2, fun, 1, 1, ifoo, ifoo, 0, ifoo, ifoo, -1, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -8)

        call sgemaxvol_proj(1, 2, fun, 1, 1, ifoo, ifoo, 3, ifoo, ifoo, -1, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -8)

        call sgemaxvol_proj(1, 2, fun, 1, 1, ifoo, ifoo, 2, ifoo, ifoo, -1, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -11)

        call sgemaxvol_proj(1, 2, fun, 1, 1, ifoo, ifoo, 2, ifoo, ifoo, 1, ONE-EPS, foo, 0, ifoo, 0, info)
        ASSERT(info == -12)
    end subroutine test_sanity

    subroutine test_quick()
        integer                :: info, ifoo(1)
        real(WP)               :: foo(1)
        procedure(MS), pointer :: fun

        call sgemaxvol_proj(0, 1, fun, 0, 0, ifoo, ifoo, 1, ifoo, ifoo, 1, ONE, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)

        call sgemaxvol_proj(1, 0, fun, 0, 1, ifoo, ifoo, 0, ifoo, ifoo, 1, ONE, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer                :: m, n, r, kr, kc, niter, ldu, ldvt, lwork, liwork, info, ifoo(1), i, j
        real(WP)               :: thresh, foo(1), err0, err1
        type(prng)             :: rng
        procedure(MS), pointer :: fun
        integer,  allocatable  :: iwork(:), icol(:), irow(:), icol_short(:), irow_short(:)
        real(WP), allocatable  :: tmp(:), U(:), VT(:), work(:)
        
        call rng%init(seed, info)
        ASSERT(info == 0)

        fun => matslc

    !-- improves with iterations
        m = 20
        n = 20
        r = 2
        kr = 2*r
        kc = 2*r
        niter = 20
        thresh = 1.02_WP
        ldu = m + 1
        ldvt = r + 1

        call sgemaxvol(m, n, fun, r, ifoo, ifoo, niter, thresh, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)
        call sgemaxvol_proj(m, n, fun, r, kr, ifoo, ifoo, kc, ifoo, ifoo, niter, thresh, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
        call smatcross(m, n, fun, kc, ifoo, kr, ifoo, r, foo, ldu, foo, ldvt, foo, -1, ifoo, -1, info)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        allocate(A(m*n), tmp(m*n), icol(n), irow(m), icol_short(n), irow_short(m), &
            U(ldu*r), VT(ldvt*n), work(lwork), iwork(liwork))
        
        do j = 1, n
            do i = 1, m
                A(i + (j-1)*m) = smatval_ballistic(m, n, i, j, info)
            end do
        end do

        call arange(n, icol, 1, 1, info)
        call arange(m, irow, 1, 1, info) 
        call permute(rng, n, icol, info)
        call permute(rng, m, irow, info)
        call smatcross(m, n, fun, kc, icol, kr, irow, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
        call slacpy('a', m, n, A, m, tmp, m)
        call sgemm('n', 'n', m, n, r, -ONE, U, ldu, VT, ldvt, ONE, tmp, m)
        err0 = sgenrmc(m, n, tmp, m, info)

        call arange(n, icol_short, 1, 1, info)
        call arange(m, irow_short, 1, 1, info)
        call permute(rng, m, irow_short, info)
        call permute(rng, n, icol_short, info)
        call sgemaxvol_proj(m, n, fun, r, kr, irow, icol_short, kc, icol, irow_short, niter, thresh, work, lwork, iwork, liwork, info)
        ASSERT(info == 0)
        call smatcross(m, n, fun, kc, icol, kr, irow, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
        call slacpy('a', m, n, A, m, tmp, m)
        call sgemm('n', 'n', m, n, r, -ONE, U, ldu, VT, ldvt, ONE, tmp, m)
        err1 = sgenrmc(m, n, tmp, m, info)
        ASSERT(safe_less(err1, err0))

        deallocate(A, tmp, icol, irow, icol_short, irow_short, U, VT, work, iwork)
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
end program sgemaxvol_proj_test
