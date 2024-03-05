program srandort_test
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
use maria_la_core_mod,      only: &
    sgemm
use maria_la_utils_mod,     only: &
    srandort,                     &
    geall_const
implicit none (type, external)

    character(100) :: seed_str
    integer        :: seed

    call get_command_argument(1, seed_str)
    read(seed_str, *) seed

    call test_sanity(seed)
    call test_quick(seed)
    call test(seed)

contains
    subroutine test_sanity(seed)
        integer, intent(in) :: seed

        integer    :: info
        real(WP)   :: foo(1)
        type(prng) :: rng

        call rng%init(seed, info)
        ASSERT(info == 0)

        call srandort(rng, -1, -1, foo, 0, foo, -1, info)
        ASSERT(info == -2)

        call srandort(rng, 2, -1, foo, 0, foo, -1, info)
        ASSERT(info == -3)

        call srandort(rng, 2, 1, foo, 2, foo, 0, info)
        ASSERT(info == -7)

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test_sanity

    subroutine test_quick(seed)
        integer, intent(in) :: seed

        integer    :: info
        real(WP)   :: foo(1)
        type(prng) :: rng

        call rng%init(seed, info)
        ASSERT(info == 0)

        call srandort(rng, 1, 0, foo, 1, foo, 0, info)
        ASSERT(info == 0)

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer               :: info, m, n, ldQ, lwork
        real(WP)              :: foo(1)
        type(prng)            :: rng       
        real(WP), allocatable :: Q(:), work(:), QTQ(:)
 
        call rng%init(seed, info)
        ASSERT(info == 0)

    !-- ldQ > m, m > n ----------------------------------------------------------------
        m = 20
        n = 10
        ldQ = m + 1
        call srandort(rng, m, n, foo, ldQ, foo, -1, info)
        lwork = int(foo(1))
        allocate(Q(ldQ*n), work(lwork), QTQ(n*n))
        call srandort(rng, m, n, Q, ldQ, work, lwork, info)
        ASSERT(info == 0)
        call sgemm('t', 'n', n, n, m, ONE, Q, ldQ, Q, ldQ, ZERO, QTQ, n)
        ASSERT(geall_const(n, n, QTQ, n, ZERO, ONE, info, atol=1e-5_WP))
        deallocate(Q, work, QTQ)

    !-- ldQ = m, m > n ----------------------------------------------------------------
        m = 20
        n = 10
        ldQ = m
        call srandort(rng, m, n, foo, ldQ, foo, -1, info)
        lwork = int(foo(1))
        allocate(Q(ldQ*n), work(lwork), QTQ(n*n))
        call srandort(rng, m, n, Q, ldQ, work, lwork, info)
        ASSERT(info == 0)
        call sgemm('t', 'n', n, n, m, ONE, Q, ldQ, Q, ldQ, ZERO, QTQ, n)
        ASSERT(geall_const(n, n, QTQ, n, ZERO, ONE, info, atol=1e-5_WP))
        deallocate(Q, work, QTQ)

    !-- ldQ > m, m < n ----------------------------------------------------------------
        m = 10
        n = 20
        ldQ = m + 1
        call srandort(rng, m, n, foo, ldQ, foo, -1, info)
        lwork = int(foo(1))
        allocate(Q(ldQ*n), work(lwork), QTQ(m*m))
        call srandort(rng, m, n, Q, ldQ, work, lwork, info)
        ASSERT(info == 0)
        call sgemm('n', 't', m, m, n, ONE, Q, ldQ, Q, ldQ, ZERO, QTQ, m)
        ASSERT(geall_const(m, m, QTQ, m, ZERO, ONE, info, atol=1e-5_WP))
        deallocate(Q, work, QTQ)

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test
end program srandort_test
