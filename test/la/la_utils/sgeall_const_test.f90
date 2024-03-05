program sgeall_const_test
#include "maria_assert_mod.h"
use maria_la_utils_mod,  only: &
    sgeall_const
use maria_kinds_mod,     only: &
    WP => SP
use maria_constants_mod, only: &
    HALF => S_HALF,            &
    EPS => S_MACHTOL,          &
    ZERO => S_ZERO
use maria_la_core_mod,   only: &
    slaset
implicit none (type, external)

    call test_sanity()
    call test_quick()
    call test()

contains
    subroutine test_sanity()
        integer  :: info
        real(WP) :: foo(1)

        ASSERT(.not. sgeall_const(-1, -1, foo, 0, ZERO, HALF, info))
        ASSERT(info == -1)

        ASSERT(.not. sgeall_const(1, -1, foo, 0, ZERO, HALF, info))
        ASSERT(info == -2)

        ASSERT(.not. sgeall_const(1, 1, foo, 0, ZERO, HALF, info))
        ASSERT(info == -4)

        ASSERT(.not. sgeall_const(1, 1, foo, 1, ZERO, HALF, info, atol=-EPS))
        ASSERT(info == -8)

        ASSERT(.not. sgeall_const(1, 1, foo, 1, ZERO, HALF, info, rtol=-EPS))
        ASSERT(info == -9)
    end subroutine test_sanity

    subroutine test_quick()
        integer  :: info
        real(WP) :: foo(1)

        ASSERT(sgeall_const(0, 1, foo, 1, ZERO, HALF, info))
        ASSERT(info == 0)

        ASSERT(sgeall_const(1, 0, foo, 1, ZERO, HALF, info))
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test()
    implicit none
        integer               :: m, n, ldA, info
        real(WP)              :: alpha, beta
        real(WP), allocatable :: A(:)

        m = 3
        n = 2
        ldA = m + 1

        allocate(A(ldA*n))
        A = HALF + EPS
        alpha = HALF
        beta = EPS
        call slaset('a', m, n, alpha, beta, A, ldA)

        ASSERT(sgeall_const(m, n, A, ldA, alpha, beta, info))
        ASSERT(info == 0)

        ASSERT(.not. sgeall_const(ldA, n, A, ldA, alpha, beta, info))
        ASSERT(info == 0)

        ASSERT(sgeall_const(ldA, n, A, ldA, alpha, beta, info, atol=2*EPS))
        ASSERT(info == 0)

        ASSERT(sgeall_const(ldA, n, A, ldA, alpha, beta, info, rtol=4*EPS))
        ASSERT(info == 0)
    end subroutine test
end program sgeall_const_test
