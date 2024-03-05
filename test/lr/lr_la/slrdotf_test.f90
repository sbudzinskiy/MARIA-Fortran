program slrdotf_test
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
use maria_utils_mod,        only: &
    are_close
use maria_la_core_mod,      only: &
    sgemm,                        &
    sgedotf
use maria_lr_la_mod,        only: &
    slrdotf
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
        real(WP)   :: dot, foo(1)

        dot = slrdotf(-1, -1, -1, foo, 0, foo, 0, -1, foo, 0, foo, 0, foo, 0, info)
        ASSERT(info == -1)

        dot = slrdotf(1, -1, -1, foo, 0, foo, 0, -1, foo, 0, foo, 0, foo, 0, info)
        ASSERT(info == -2)

        dot = slrdotf(1, 1, -1, foo, 0, foo, 0, -1, foo, 0, foo, 0, foo, 0, info)
        ASSERT(info == -3)

        dot = slrdotf(1, 1, 1, foo, 0, foo, 0, -1, foo, 0, foo, 0, foo, 0, info)
        ASSERT(info == -5)

        dot = slrdotf(1, 1, 1, foo, 1, foo, 0, -1, foo, 0, foo, 0, foo, 0, info)
        ASSERT(info == -7)

        dot = slrdotf(1, 1, 1, foo, 1, foo, 1, -1, foo, 0, foo, 0, foo, 0, info)
        ASSERT(info == -8)

        dot = slrdotf(1, 1, 1, foo, 1, foo, 1, 1, foo, 0, foo, 0, foo, 0, info)
        ASSERT(info == -10)

        dot = slrdotf(1, 1, 1, foo, 1, foo, 1, 1, foo, 1, foo, 0, foo, 0, info)
        ASSERT(info == -12)
    end subroutine test_sanity

    subroutine test_quick()
        integer    :: info
        real(WP)   :: dot, foo(1)

        dot = slrdotf(0, 1, 1, foo, 1, foo, 1, 1, foo, 1, foo, 1, foo, 0, info)
        ASSERT(info == 0)

        dot = slrdotf(1, 0, 1, foo, 1, foo, 1, 1, foo, 1, foo, 1, foo, 0, info)
        ASSERT(info == 0)

        dot = slrdotf(1, 1, 0, foo, 1, foo, 1, 1, foo, 1, foo, 1, foo, 0, info)
        ASSERT(info == 0)

        dot = slrdotf(1, 1, 1, foo, 1, foo, 1, 0, foo, 1, foo, 1, foo, 0, info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer               :: info, m, n, r, k, ldU, ldVT, ldA, ldBT, lwork
        real(WP)              :: foo(1), dot1, dot2
        type(prng)            :: rng       
        real(WP), allocatable :: U(:), VT(:), A(:), BT(:), work(:), C(:), D(:)
 
        call rng%init(seed, info)
        ASSERT(info == 0)

        m = 20
        n = 10
        r = 15
        k = 8
        ldU = m + 1
        ldVT = r + 1
        ldA = m + 1
        ldBT = k + 1
        dot1 = slrdotf(m, n, r, foo, ldU, foo, ldVT, k, foo, ldA, foo, ldBT, foo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        allocate(U(ldU*r), VT(ldVT*n), A(ldA*k), BT(ldBT*n), work(lwork), C(m*n), D(m*n))
        call rng%snormal(ldU*r, U, ZERO, ONE, info)
        call rng%snormal(ldVT*n, VT, ZERO, ONE, info)
        call rng%snormal(ldA*k, A, ZERO, ONE, info)
        call rng%snormal(ldBT*n, BT, ZERO, ONE, info)
        dot1 = slrdotf(m, n, r, U, ldU, VT, ldVT, k, A, ldA, BT, ldBT, work, lwork, info)
        ASSERT(info == 0)
        call sgemm('n', 'n', m, n, r, ONE, U, ldU, VT, ldVT, ZERO, C, m)
        call sgemm('n', 'n', m, n, k, ONE, A, ldA, BT, ldBT, ZERO, D, m)
        dot2 = sgedotf('n', m, n, C, m, D, m, info)
        ASSERT(are_close(dot1, dot2, info, rtol=1e-5_WP))

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test
end program slrdotf_test
