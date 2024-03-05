program sgeall_close_test
#include "maria_assert_mod.h"
use maria_la_utils_mod,  only: &
    sgeall_close
use maria_kinds_mod,     only: &
    WP => SP
use maria_constants_mod, only: &
    HALF => S_HALF,            &
    EPS => S_MACHTOL
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

        ASSERT(.not. sgeall_close(-1, -1, foo, 0, foo, 0, info))
        ASSERT(info == -1)

        ASSERT(.not. sgeall_close(1, -1, foo, 0, foo, 0, info))
        ASSERT(info == -2)

        ASSERT(.not. sgeall_close(1, 1, foo, 0, foo, 0, info))
        ASSERT(info == -4)

        ASSERT(.not. sgeall_close(1, 1, foo, 1, foo, 0, info))
        ASSERT(info == -6)

        ASSERT(.not. sgeall_close(1, 1, foo, 1, foo, 1, info, atol=-EPS))
        ASSERT(info == -8)

        ASSERT(.not. sgeall_close(1, 1, foo, 1, foo, 1, info, rtol=-EPS))
        ASSERT(info == -9)
    end subroutine test_sanity

    subroutine test_quick()
        integer  :: info
        real(WP) :: foo(1)

        ASSERT(sgeall_close(0, 1, foo, 1, foo, 1, info))
        ASSERT(info == 0)

        ASSERT(sgeall_close(1, 0, foo, 1, foo, 1, info))
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test()
        integer               :: m, n, ldA, ldB, info
        real(WP), allocatable :: A(:), B(:)

        m = 3
        n = 2
        ldA = m + 2
        ldB = m + 2

        allocate(A(ldA*n), B(ldB*n))
        A = HALF
        B = HALF
        call slaset('a', ldA-m, n, HALF+EPS, HALF+EPS, A(m+1:), ldA)
        call slaset('a', ldB-m, n, HALF-EPS, HALF-EPS, B(m+1:), ldB)

        ASSERT(sgeall_close(m, n, A, ldA, B, ldB, info))
        ASSERT(info == 0)

        ASSERT(.not. sgeall_close(m+1, n, A, ldA, B, ldB, info))
        ASSERT(info == 0)

        ASSERT(sgeall_close(m+1, n, A, ldA, B, ldB, info, atol=3*EPS))
        ASSERT(info == 0)

        ASSERT(sgeall_close(m+1, n, A, ldA, B, ldB, info, rtol=6*EPS))
        ASSERT(info == 0)

        ASSERT(.not. sgeall_close(ldA, n, A, ldA, B, ldB, info))
        ASSERT(info == 0)

        ASSERT(sgeall_close(ldA, n, A, ldA, B, ldB, info, atol=3*EPS))
        ASSERT(info == 0)

        ASSERT(sgeall_close(ldA, n, A, ldA, B, ldB, info, rtol=6*EPS))
        ASSERT(info == 0)
    end subroutine test
end program sgeall_close_test
