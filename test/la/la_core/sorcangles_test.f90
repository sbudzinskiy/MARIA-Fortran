program sorcangles_test
#include "maria_assert_mod.h"
#ifdef MARIA_MKL
use maria_prng_mkl_mod,     only: &
    prng => prng_mkl
#else
use maria_prng_builtin_mod, only: &
    prng => prng_builtin
#endif
use maria_la_core_mod,      only: &
    sorcangles,                   &
    sgeqrf,                       &
    sorgqr,                       &
    sgelqf,                       &
    sorglq
use maria_kinds_mod,        only: &
    WP => SP
use maria_constants_mod,    only: &
    ZERO => S_ZERO,               &
    ONE => S_ONE
use maria_la_utils_mod,     only: &
    all_const
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
        integer  :: info, ifoo(1)
        real(WP) :: foo(1)
        
        call sorcangles('a', -1, -1, foo, 1, foo, 1, foo, foo, 0, ifoo, 0, info)
        ASSERT(info == -1)

        call sorcangles('c', -1, -1, foo, 1, foo, 1, foo, foo, 0, ifoo, 0, info)
        ASSERT(info == -2)

        call sorcangles('c', 2, -1, foo, 1, foo, 1, foo, foo, 0, ifoo, 0, info)
        ASSERT(info == -3)

        call sorcangles('c', 2, 3, foo, 1, foo, 1, foo, foo, 0, ifoo, 0, info)
        ASSERT(info == -3)

        call sorcangles('r', 2, 1, foo, 1, foo, 1, foo, foo, 0, ifoo, 0, info)
        ASSERT(info == -3)

        call sorcangles('c', 2, 2, foo, 1, foo, 1, foo, foo, 0, ifoo, 0, info)
        ASSERT(info == -5)

        call sorcangles('c', 2, 2, foo, 2, foo, 1, foo, foo, 0, ifoo, 0, info)
        ASSERT(info == -7)
    end subroutine test_sanity

    subroutine test_quick()
        integer  :: info, ifoo(1)
        real(WP) :: foo(1)
        
        call sorcangles('r', 0, 1, foo, 1, foo, 1, foo, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)

        call sorcangles('c', 1, 0, foo, 1, foo, 1, foo, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        character(1)          :: what
        integer               :: info, m, n, ldA, lwork, liwork, k, ifoo(1)
        real(WP)              :: foo(1)
        type(prng)            :: rng
        integer,  allocatable :: iwork(:)
        real(WP), allocatable :: A(:), work(:), cosines(:), tau(:)

        call rng%init(seed, info)
        ASSERT(info == 0)

    !-- what = 'c' ----------------------------------------------------------------
        what = 'c'
        m = 10
        n = 8
        k = n / 2
        lda = m + 1
        call sgeqrf(m, n, foo, ldA, foo, foo, -1, info)
        lwork = int(foo(1))
        call sorgqr(m, n, n, foo, ldA, foo, foo, -1, info)
        lwork = max(lwork, int(foo(1)))
        call sorcangles(what, m, n, foo, ldA, foo, ldA, foo, foo, -1, ifoo, -1, info)
        lwork = max(lwork, int(foo(1)))
        liwork = ifoo(1)
        allocate(A(lda*n), tau(n), work(lwork), iwork(liwork), cosines(k))
        call rng%snormal(lda*n, A, ZERO, ONE, info)
        call sgeqrf(m, n, A, ldA, tau, work, lwork, info)
        call sorgqr(m, n, n, A, ldA, tau, work, lwork, info)
        associate(U => A, V => A(1 + k*ldA:))
            call sorcangles(what, m, k, U, ldA, V, ldA, cosines, work, lwork, iwork, liwork, info)
        end associate
        ASSERT(info == 0)
        ASSERT(all_const(k, cosines, 1, ZERO, info, atol=1e-5_WP))
        associate(U => A, V => A)
            call sorcangles(what, m, k, U, ldA, V, ldA, cosines, work, lwork, iwork, liwork, info)
        end associate
        ASSERT(info == 0)
        ASSERT(all_const(k, cosines, 1, ONE, info, atol=1e-5_WP))
        deallocate(A, tau, work, iwork, cosines)

    !-- what = 'r' ----------------------------------------------------------------
        what = 'r'
        m = 8
        n = 10
        k = m / 2
        lda = m + 1
        call sgelqf(m, n, foo, ldA, foo, foo, -1, info)
        lwork = int(foo(1))
        call sorglq(m, n, m, foo, ldA, foo, foo, -1, info)
        lwork = max(lwork, int(foo(1)))
        call sorcangles(what, k, n, foo, ldA, foo, ldA, foo, foo, -1, ifoo, -1, info)
        lwork = max(lwork, int(foo(1)))
        liwork = ifoo(1)
        allocate(A(lda*n), tau(m), work(lwork), iwork(liwork), cosines(k))
        call rng%snormal(lda*n, A, ZERO, ONE, info)
        call sgelqf(m, n, A, ldA, tau, work, lwork, info)
        call sorglq(m, n, m, A, ldA, tau, work, lwork, info)
        associate(U => A, V => A(k+1:))
            call sorcangles(what, k, n, U, ldA, V, ldA, cosines, work, lwork, iwork, liwork, info)
        end associate
        ASSERT(info == 0)
        ASSERT(all_const(k, cosines, 1, ZERO, info, atol=1e-5_WP))
        associate(U => A, V => A)
            call sorcangles(what, k, n, U, ldA, V, ldA, cosines, work, lwork, iwork, liwork, info)
        end associate
        ASSERT(info == 0)
        ASSERT(all_const(k, cosines, 1, ONE, info, atol=1e-5_WP))

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test
end program sorcangles_test
