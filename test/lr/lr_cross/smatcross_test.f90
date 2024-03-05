program smatcross_test
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
    ZERO => S_ZERO
use maria_access_matrix_mod, only: &
    MS => smatslc
use maria_comparison_mod,   only: &
    safe_leq,                     &
    safe_less
use maria_utils_mod,        only: &
    arange,                       &
    permute
use maria_la_utils_mod,     only: &
    sgeall_const
use maria_la_core_mod,      only: &
    slacpy,                       &
    slaset,                       &
    sgemm,                        &
    sgenrmc
use maria_lr_cross_mod,     only: &
    smatcross
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

        call smatcross(-1, -1, fun, -1, ifoo, -1, ifoo, -1, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -1)

        call smatcross(1, -1, fun, -1, ifoo, -1, ifoo, -1, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -2)

        call smatcross(1, 1, fun, -1, ifoo, -1, ifoo, -1, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -4)

        call smatcross(1, 1, fun, 2, ifoo, -1, ifoo, -1, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -4)

        call smatcross(1, 1, fun, 1, ifoo, -1, ifoo, -1, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -6)

        call smatcross(1, 1, fun, 1, ifoo, 2, ifoo, -1, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -6)

        call smatcross(1, 1, fun, 1, ifoo, 1, ifoo, -1, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -8)

        call smatcross(1, 1, fun, 1, ifoo, 1, ifoo, 2, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -8)

        call smatcross(1, 1, fun, 1, ifoo, 1, ifoo, 1, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -10)

        call smatcross(1, 1, fun, 1, ifoo, 1, ifoo, 1, foo, 1, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -12)
    end subroutine test_sanity

    subroutine test_quick()
        integer                :: info, ifoo(1)
        real(WP)               :: foo(1)
        procedure(MS), pointer :: fun

        call smatcross(0, 1, fun, 1, ifoo, 0, ifoo, 0, foo, 1, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)

        call smatcross(1, 0, fun, 0, ifoo, 1, ifoo, 0, foo, 1, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)

        call smatcross(1, 1, fun, 0, ifoo, 1, ifoo, 0, foo, 1, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)

        call smatcross(1, 1, fun, 1, ifoo, 0, ifoo, 0, foo, 1, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)

        call smatcross(1, 1, fun, 1, ifoo, 1, ifoo, 0, foo, 1, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer                :: m, n, nc, nr, r, ldu, ldvt, lwork, liwork, info, ifoo(1)
        real(WP)               :: nrmc, foo(1)
        type(prng)             :: rng
        procedure(MS), pointer :: fun
        integer,  allocatable  :: iwork(:), icol(:), irow(:)
        real(WP), allocatable  :: U(:), VT(:), work(:)
        
        call rng%init(seed, info)
        ASSERT(info == 0)

        fun => matslc

    !-- nc = nr = r
        m = 15
        n = 20
        nc = 4
        nr = nc
        r = nc
        ldu = m + 1
        ldvt = r + 1

        call smatcross(m, n, fun, nc, ifoo, nr, ifoo, r, foo, ldu, foo, ldvt, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)

        allocate(A(m*n), icol(n), irow(m), U(ldu*r), VT(ldvt*n), work(lwork), iwork(liwork))
        
        call rng%snormal(ldu*r, U, ZERO, ONE, info)
        call rng%snormal(ldvt*n, VT, ZERO, ONE, info)
        call arange(n, icol, 1, 1, info)
        call arange(m, irow, 1, 1, info)
        call permute(rng, n, icol, info)
        call permute(rng, m, irow, info)

        call sgemm('n', 'n', m, n, r, ONE, U, ldu, VT, ldvt, ZERO, A, m)
        nrmc = sgenrmc(m, n, A, m, info)
        
        call smatcross(m, n, fun, nc, icol, nr, irow, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
        ASSERT(info == 0)
        call sgemm('n', 'n', m, n, r, -ONE, U, ldu, VT, ldvt, ONE, A, m)
        ASSERT(sgeall_const(m, n, A, m, ZERO, ZERO, info, atol=1e-5_WP * nrmc))

        deallocate(A, icol, irow, U, VT, work, iwork)

    !-- r < nc < nr
        m = 15
        n = 20
        nc = 4
        nr = nc + 1
        r = nc - 1
        ldu = m + 1
        ldvt = r + 1

        call smatcross(m, n, fun, nc, ifoo, nr, ifoo, r, foo, ldu, foo, ldvt, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)

        allocate(A(m*n), icol(n), irow(m), U(ldu*r), VT(ldvt*n), work(lwork), iwork(liwork))
        
        call rng%snormal(ldu*r, U, ZERO, ONE, info)
        call rng%snormal(ldvt*n, VT, ZERO, ONE, info)
        call arange(n, icol, 1, 1, info)
        call arange(m, irow, 1, 1, info)
        call permute(rng, n, icol, info)
        call permute(rng, m, irow, info)

        call sgemm('n', 'n', m, n, r, ONE, U, ldu, VT, ldvt, ZERO, A, m)
        nrmc = sgenrmc(m, n, A, m, info)
        
        call smatcross(m, n, fun, nc, icol, nr, irow, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
        ASSERT(info == 0)
        call sgemm('n', 'n', m, n, r, -ONE, U, ldu, VT, ldvt, ONE, A, m)
        ASSERT(sgeall_const(m, n, A, m, ZERO, ZERO, info, atol=1e-5_WP * nrmc))

        deallocate(A, icol, irow, U, VT, work, iwork)

    !-- nc = nr = r, runtimer error
        m = 15
        n = 20
        nc = 4
        nr = nc
        r = nc
        ldu = m + 1
        ldvt = r + 1

        call smatcross(m, n, fun, nc, ifoo, nr, ifoo, r, foo, ldu, foo, ldvt, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)

        allocate(A(m*n), icol(n), irow(m), U(ldu*r), VT(ldvt*n), work(lwork), iwork(liwork))
 
        call rng%snormal(ldu*r, U, ZERO, ONE, info)
        call rng%snormal(ldvt*n, VT, ZERO, ONE, info)       
        call slaset('a', r, r, ZERO, ZERO, U, ldu)
        call sgemm('n', 'n', m, n, r, ONE, U, ldu, VT, ldvt, ZERO, A, m)

        call arange(n, icol, 1, 1, info)
        call arange(m, irow, 1, 1, info)

        call smatcross(m, n, fun, nc, icol, nr, irow, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
        ASSERT(info > 0)

        deallocate(A, icol, irow, U, VT, work, iwork)

    !-- r < nc < nr, runtimer error
        m = 15
        n = 20
        nc = 4
        nr = nc + 1
        r = nc - 1
        ldu = m + 1
        ldvt = r + 1

        call smatcross(m, n, fun, nc, ifoo, nr, ifoo, r, foo, ldu, foo, ldvt, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)

        allocate(A(m*n), icol(n), irow(m), U(ldu*r), VT(ldvt*n), work(lwork), iwork(liwork))
 
        call rng%snormal(ldu*r, U, ZERO, ONE, info)
        call rng%snormal(ldvt*n, VT, ZERO, ONE, info)       
        call slaset('a', r, r, ZERO, ZERO, U, ldu)
        call sgemm('n', 'n', m, n, r, ONE, U, ldu, VT, ldvt, ZERO, A, m)

        call arange(n, icol, 1, 1, info)
        call arange(m, irow, 1, 1, info)

        call smatcross(m, n, fun, nc, icol, nr, irow, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
        ASSERT(info > 0)

        deallocate(A, icol, irow, U, VT, work, iwork)
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
end program smatcross_test
