program smatcross_aca_test
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
    MV => smatval,                 &
    smatval2slc,                   &
    smatval2full,                  &
    smatval_ballistic
use maria_comparison_mod,   only: &
    safe_leq,                     &
    safe_less
use maria_utils_mod,        only: &
    arange,                       &
    permute
use maria_la_utils_mod,     only: &
    sgeall_const
use maria_la_core_mod,      only: &
    sgemm,                        &
    sgenrmf
use maria_lr_la_mod,        only: &
    slrnrmf
use maria_lr_cross_mod,     only: &
    smatcross_aca,                &
    smatcross
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
        integer                :: info(2), ifoo(1)
        real(WP)               :: foo(1)
        procedure(MV), pointer :: matval

        call smatcross_aca(-1, -1, matval, ifoo, -1, -1, ifoo(1), foo, 0, foo, 0, ifoo, ifoo, foo, 0, info)
        ASSERT(info(1) == -1)

        call smatcross_aca(1, -1, matval, ifoo, -1, -1, ifoo(1), foo, 0, foo, 0, ifoo, ifoo, foo, 0, info)
        ASSERT(info(1) == -2)

        call smatcross_aca(1, 1, matval, ifoo, -1, -1, ifoo(1), foo, 0, foo, 0, ifoo, ifoo, foo, 0, info)
        ASSERT(info(1) == -5)

        call smatcross_aca(1, 1, matval, ifoo, 2, -1, ifoo(1), foo, 0, foo, 0, ifoo, ifoo, foo, 0, info)
        ASSERT(info(1) == -5)

        call smatcross_aca(1, 1, matval, ifoo, 1, -1, ifoo(1), foo, 0, foo, 0, ifoo, ifoo, foo, 0, info)
        ASSERT(info(1) == -6)

        call smatcross_aca(1, 1, matval, ifoo, 1, 1, ifoo(1), foo, 0, foo, 0, ifoo, ifoo, foo, 0, info)
        ASSERT(info(1) == -9)

        call smatcross_aca(1, 1, matval, ifoo, 1, 1, ifoo(1), foo, 1, foo, 0, ifoo, ifoo, foo, 0, info)
        ASSERT(info(1) == -11)

        call smatcross_aca(1, 1, matval, ifoo, 1, 1, ifoo(1), foo, 1, foo, 1, ifoo, ifoo, foo, 0, info, rtol=-EPS)
        ASSERT(info(1) == -17)
    end subroutine test_sanity

    subroutine test_quick()
        integer                :: info(2), ifoo(1)
        real(WP)               :: foo(1)
        procedure(MV), pointer :: matval

        call smatcross_aca(0, 1, matval, ifoo, 0, 1, ifoo(1), foo, 1, foo, 1, ifoo, ifoo, foo, 0, info)
        ASSERT(info(1) == 0)
        ASSERT(ifoo(1) == 0)

        call smatcross_aca(1, 0, matval, ifoo, 0, 1, ifoo(1), foo, 1, foo, 1, ifoo, ifoo, foo, 0, info)
        ASSERT(info(1) == 0)
        ASSERT(ifoo(1) == 0)

        call smatcross_aca(1, 1, matval, ifoo, 0, 1, ifoo(1), foo, 1, foo, 1, ifoo, ifoo, foo, 0, info)
        ASSERT(info(1) == 0)
        ASSERT(ifoo(1) == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer                :: m, n, maxr, niter, r, ldu, ldvt, lwork, info(2), ifoo(1), liwork
        real(WP)               :: foo(1), rtol, err1, err2, nrm, nrmdiff
        type(prng)             :: rng
        procedure(MV), pointer :: matval
        procedure(MS), pointer :: matslc
        integer,  allocatable  :: col_order(:), icol(:), irow(:), iwork(:)
        real(WP), allocatable  :: A(:), U(:), VT(:), work(:), UU(:), VVT(:)
        
        call rng%init(seed, info(1))
        ASSERT(info(1) == 0)

        matval => smatval_ballistic
        matslc => matslc_ballistic

    !-- Effect of niter
        m = 15
        n = 20
        maxr = 8
        niter = 0
        ldu = m + 1
        ldvt = maxr + 1

        call smatcross_aca(m, n, matval, ifoo, maxr, niter, r, foo, ldu, foo, ldvt, ifoo, ifoo, foo, -1, info)
        ASSERT(info(1) == 0)
        lwork = int(foo(1))

        allocate(A(m*n), col_order(n), U(ldu*maxr), VT(ldvt*n), irow(maxr), icol(maxr), work(lwork))
        
        call arange(n, col_order, 1, 1, info(1))
        call permute(rng, n, col_order, info(1))

        call smatcross_aca(m, n, matval, col_order, maxr, niter, r, U, ldu, VT, ldVT, irow, icol, work, lwork, info)
        ASSERT(info(1) == 0)
        ASSERT(r == maxr)
        call smatval2full(matval, 'n', m, n, A, m, info(1))
        call sgemm('n', 'n', m, n, r, -ONE, U, ldu, VT, ldvt, ONE, A, m)
        err1 = sgenrmf(m, n, A, m, info(1))

        niter = 5        
        call smatcross_aca(m, n, matval, col_order, maxr, niter, r, U, ldu, VT, ldVT, irow, icol, work, lwork, info)
        ASSERT(info(1) == 0)
        ASSERT(r == maxr)
        call smatval2full(matval, 'n', m, n, A, m, info(1))
        call sgemm('n', 'n', m, n, r, -ONE, U, ldu, VT, ldvt, ONE, A, m)
        err2 = sgenrmf(m, n, A, m, info(1))

        ASSERT(safe_leq(err2, err1))
        deallocate(A, col_order, U, VT, irow, icol, work)
    
    !-- Check indices
        m = 15
        n = 20
        maxr = 8
        niter = 5
        ldu = m + 1
        ldvt = maxr + 1

        call smatcross_aca(m, n, matval, ifoo, maxr, niter, r, foo, ldu, foo, ldvt, ifoo, ifoo, foo, -1, info)
        ASSERT(info(1) == 0)
        lwork = int(foo(1))
        call smatcross(m, n, matslc, maxr, ifoo, maxr, ifoo, maxr, foo, m, foo, maxr, foo, -1, ifoo, -1, info(1))
        lwork = max(lwork, int(foo(1)))
        liwork = ifoo(1)

        allocate(A(m*n), col_order(n), U(ldu*maxr), VT(ldvt*n), irow(maxr), icol(maxr), work(lwork), iwork(liwork), &
            UU(m*maxr), VVT(maxr*n))
        
        call arange(n, col_order, 1, 1, info(1))
        call permute(rng, n, col_order, info(1))

        call smatcross_aca(m, n, matval, col_order, maxr, niter, r, U, ldu, VT, ldVT, irow, icol, work, lwork, info)
        ASSERT(info(1) == 0)
        ASSERT(r == maxr)

        call smatcross(m, n, matslc, r, icol, r, irow, r, UU, m, VVT, maxr, work, lwork, iwork, liwork, info(1))

        call sgemm('n', 'n', m, n, r, ONE, U, ldu, VT, ldvt, ZERO, A, m)
        call sgemm('n', 'n', m, n, r, -ONE, UU, m, VVT, maxr, ONE, A, m)
        err1 = sgenrmf(m, n, A, m, info(1))

        call smatval2full(matval, 'n', m, n, A, m, info(1))
        nrm = sgenrmf(m, n, A, m, info(1))
        ASSERT(safe_less(err1, 1e-5_WP * nrm))
        
        deallocate(A, col_order, U, VT, irow, icol, work, iwork, UU, VVT)

    !-- rtol
        m = 15
        n = 20
        maxr = 15
        niter = 5
        ldu = m + 1
        ldvt = maxr + 1
        rtol = 1e-5_WP

        call smatcross_aca(m, n, matval, ifoo, maxr, niter, r, foo, ldu, foo, ldvt, ifoo, ifoo, foo, -1, info, rtol)
        ASSERT(info(1) == 0)
        lwork = int(foo(1))
        nrm = slrnrmf(m, n, maxr, foo, ldu, foo, ldvt, foo, -1, info(1))
        lwork = max(lwork, int(foo(1)))

        allocate(A(m*n), col_order(n), U(ldu*maxr), VT(ldvt*n), irow(maxr), icol(maxr), work(lwork))
        
        call arange(n, col_order, 1, 1, info(1))
        call permute(rng, n, col_order, info(1))

        call smatcross_aca(m, n, matval, col_order, maxr, niter, r, U, ldu, VT, ldVT, irow, icol, work, lwork, info, rtol)
        ASSERT(info(1) == 0)
        ASSERT(r < maxr)

        call sgemm('n', 'n', m, n, r, -ONE, U, ldu, VT, ldvt, ONE, A, m)
        nrm = slrnrmf(m, n, r-1, U, ldu, VT, ldvt, work, lwork, info(1))
        nrmdiff = slrnrmf(m, n, 1, U(1 + (r-1)*ldu:), ldu, VT(r:), ldvt, work, lwork, info(1))
        ASSERT(safe_leq(nrmdiff, rtol * nrm))
        deallocate(A, col_order, U, VT, irow, icol, work)
    end subroutine test

    subroutine matslc_ballistic &
    (m, n, mode, ind, x, incx, info)
        integer, intent(in) :: m
        integer, intent(in) :: n
        integer, intent(in) :: mode
        integer, intent(in) :: ind
        real(WP), intent(out), contiguous :: x(:)
        integer, intent(in) :: incx
        integer, intent(out) :: info
        
        procedure(MV), pointer :: fun

        fun => smatval_ballistic
        
        call smatval2slc(fun, m, n, mode, ind, x, incx, info)
    end subroutine
end program smatcross_aca_test
