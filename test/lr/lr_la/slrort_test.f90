program slrort_test
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
use maria_la_utils_mod,     only: &
    geall_const
use maria_la_core_mod,      only: &
    sgemm,                        &
    sorgqr,                       &
    sorglq
use maria_lr_la_mod,        only: &
    slrort,                       &
    lrort_rank
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
        integer    :: info
        real(WP)   :: foo(1)

        call slrort('x', -1, -1, -1, foo, 0, foo, 0, foo, foo, 0, info)
        ASSERT(info == -1)

        call slrort('l', -1, -1, -1, foo, 0, foo, 0, foo, foo, 0, info)
        ASSERT(info == -2)

        call slrort('l', 1, -1, -1, foo, 0, foo, 0, foo, foo, 0, info)
        ASSERT(info == -3)

        call slrort('l', 1, 1, -1, foo, 0, foo, 0, foo, foo, 0, info)
        ASSERT(info == -4)

        call slrort('l', 1, 1, 1, foo, 0, foo, 0, foo, foo, 0, info)
        ASSERT(info == -6)

        call slrort('l', 1, 1, 1, foo, 1, foo, 0, foo, foo, 0, info)
        ASSERT(info == -8)

        call slrort('l', 1, 1, 1, foo, 1, foo, 1, foo, foo, 0, info)
        ASSERT(info == -11)
    end subroutine test_sanity

    subroutine test_quick()
        integer    :: info
        real(WP)   :: foo(1)

        call slrort('l', 0, 1, 1, foo, 1, foo, 1, foo, foo, 0, info)
        ASSERT(info == 0)

        call slrort('l', 1, 0, 1, foo, 1, foo, 1, foo, foo, 0, info)
        ASSERT(info == 0)

        call slrort('l', 1, 1, 0, foo, 1, foo, 1, foo, foo, 0, info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        character(1)          :: side
        integer               :: info, m, n, r, newr, ldU, ldVT, lwork
        real(WP)              :: foo(1)
        type(prng)            :: rng       
        real(WP), allocatable :: U(:), VT(:), work(:), A(:), tau(:)
 
        call rng%init(seed, info)
        ASSERT(info == 0)

    !-- side = 'l', m > r ------------------------------------------------------
        side = 'l'
        m = 20
        n = 10
        r = 15
        ldU = m + 1
        ldVT = r + 1
        call lrort_rank(side, m, n, r, newr, info)
        call slrort(side, m, n, r, foo, ldU, foo, ldVT, foo, foo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        call sorgqr(m, newr, newr, foo, ldU, foo, foo, -1, info)
        lwork = max(lwork, int(foo(1)))
        allocate(U(ldU*r), VT(ldVT*n), A(m*n), work(lwork), tau(newr))
        call rng%snormal(ldU*r, U, ZERO, ONE, info)
        call rng%snormal(ldVT*n, VT, ZERO, ONE, info)
        call sgemm('n', 'n', m, n, r, ONE, U, ldU, VT, ldVT, ZERO, A, m)
        call slrort(side, m, n, r, U, ldU, VT, ldVT, tau, work, lwork, info)
        call sorgqr(m, newr, newr, U, ldU, tau, work, lwork, info)
        call sgemm('n', 'n', m, n, newr, ONE, U, ldU, VT, ldVT, -ONE, A, m)
        ASSERT(geall_const(m, n, A, m, ZERO, ZERO, info, atol=1e-5_WP))
        deallocate(U, VT, A, work, tau)

    !-- side = 'l', m < r ------------------------------------------------------
        side = 'l'
        m = 20
        n = 10
        r = m + 1
        ldU = m + 1
        ldVT = r + 1
        call lrort_rank(side, m, n, r, newr, info)
        call slrort(side, m, n, r, foo, ldU, foo, ldVT, foo, foo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        call sorgqr(m, newr, newr, foo, ldU, foo, foo, -1, info)
        lwork = max(lwork, int(foo(1)))
        allocate(U(ldU*r), VT(ldVT*n), A(m*n), work(lwork), tau(newr))
        call rng%snormal(ldU*r, U, ZERO, ONE, info)
        call rng%snormal(ldVT*n, VT, ZERO, ONE, info)
        call sgemm('n', 'n', m, n, r, ONE, U, ldU, VT, ldVT, ZERO, A, m)
        call slrort(side, m, n, r, U, ldU, VT, ldVT, tau, work, lwork, info)
        call sorgqr(m, newr, newr, U, ldU, tau, work, lwork, info)
        call sgemm('n', 'n', m, n, newr, ONE, U, ldU, VT, ldVT, -ONE, A, m)
        ASSERT(geall_const(m, n, A, m, ZERO, ZERO, info, atol=1e-5_WP))
        deallocate(U, VT, A, work, tau)

    !-- side = 'r', n > r ------------------------------------------------------
        side = 'r'
        m = 20
        n = 10
        r = 5
        ldU = m + 1
        ldVT = r + 1
        call lrort_rank(side, m, n, r, newr, info)
        call slrort(side, m, n, r, foo, ldU, foo, ldVT, foo, foo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        call sorgqr(m, newr, newr, foo, ldU, foo, foo, -1, info)
        lwork = max(lwork, int(foo(1)))
        allocate(U(ldU*r), VT(ldVT*n), A(m*n), work(lwork), tau(newr))
        call rng%snormal(ldU*r, U, ZERO, ONE, info)
        call rng%snormal(ldVT*n, VT, ZERO, ONE, info)
        call sgemm('n', 'n', m, n, r, ONE, U, ldU, VT, ldVT, ZERO, A, m)
        call slrort(side, m, n, r, U, ldU, VT, ldVT, tau, work, lwork, info)
        call sorglq(newr, n, newr, VT, ldVT, tau, work, lwork, info)
        call sgemm('n', 'n', m, n, newr, ONE, U, ldU, VT, ldVT, -ONE, A, m)
        ASSERT(geall_const(m, n, A, m, ZERO, ZERO, info, atol=1e-5_WP))
        deallocate(U, VT, A, work, tau)

    !-- side = 'r', n < r ------------------------------------------------------
        side = 'r'
        m = 20
        n = 10
        r = n + 1
        ldU = m + 1
        ldVT = r + 1
        call lrort_rank(side, m, n, r, newr, info)
        call slrort(side, m, n, r, foo, ldU, foo, ldVT, foo, foo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        call sorgqr(m, newr, newr, foo, ldU, foo, foo, -1, info)
        lwork = max(lwork, int(foo(1)))
        allocate(U(ldU*r), VT(ldVT*n), A(m*n), work(lwork), tau(newr))
        call rng%snormal(ldU*r, U, ZERO, ONE, info)
        call rng%snormal(ldVT*n, VT, ZERO, ONE, info)
        call sgemm('n', 'n', m, n, r, ONE, U, ldU, VT, ldVT, ZERO, A, m)
        call slrort(side, m, n, r, U, ldU, VT, ldVT, tau, work, lwork, info)
        call sorglq(newr, n, newr, VT, ldVT, tau, work, lwork, info)
        call sgemm('n', 'n', m, n, newr, ONE, U, ldU, VT, ldVT, -ONE, A, m)
        ASSERT(geall_const(m, n, A, m, ZERO, ZERO, info, atol=1e-5_WP))

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test
end program slrort_test
