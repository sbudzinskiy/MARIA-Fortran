program srandsvd_test
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
    ONE => S_ONE
use maria_utils_mod,        only: &
    arange
use maria_la_core_mod,      only: &
    sgesdd_q
use maria_la_utils_mod,     only: &
    srandsvd,                     &
    all_close
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

        call srandsvd(rng, -1, -1, foo, 0, foo, foo, -1, info)
        ASSERT(info == -2)

        call srandsvd(rng, 2, -1, foo, 0, foo, foo, -1, info)
        ASSERT(info == -3)

        call srandsvd(rng, 2, 1, foo, 0, foo, foo, 0, info)
        ASSERT(info == -5)

        call srandsvd(rng, 2, 1, foo, 1, foo, foo, 0, info)
        ASSERT(info == -5)

        call srandsvd(rng, 2, 1, foo, 2, foo, foo, 0, info)
        ASSERT(info == -8)

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

        call srandsvd(rng, 1, 0, foo, 1, foo, foo, 0, info)
        ASSERT(info == 0)

        call srandsvd(rng, 0, 1, foo, 1, foo, foo, 0, info)
        ASSERT(info == 0)

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer               :: info, m, n, ldA, lwork, liwork, ifoo(1)
        real(WP)              :: foo(1)
        type(prng)            :: rng       
        integer,  allocatable :: iwork(:)
        real(WP), allocatable :: A(:), S(:), SS(:), work(:)
 
        call rng%init(seed, info)
        ASSERT(info == 0)

        m = 20
        n = 10
        ldA = m + 1
        call srandsvd(rng, m, n, foo, ldA, foo, foo, -1, info)
        lwork = int(foo(1))
        call sgesdd_q('n', m, n, foo, ldA, foo, foo, 1, foo, 1, foo, -1, ifoo, -1, info)
        lwork = max(lwork, int(foo(1)))
        liwork = ifoo(1)
        allocate(A(lda*n), S(n), SS(n), work(lwork), iwork(liwork))
        call arange(n, S, ONE, ONE, info)
        S = ONE / S
        call srandsvd(rng, m, n, A, ldA, S, work, lwork, info)
        ASSERT(info == 0)
        call sgesdd_q('n', m, n, A, ldA, SS, foo, 1, foo, 1, work, lwork, iwork, liwork, info)
        write(*, '(10es12.5)') S
        write(*, '(10es12.5)') SS
        ASSERT(all_close(n, S, 1, SS, 1, info, rtol=1e-5_WP))

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test
end program srandsvd_test
