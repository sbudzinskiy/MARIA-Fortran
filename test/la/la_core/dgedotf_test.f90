program dgedotf_test
#include "maria_assert_mod.h"
use maria_la_core_mod,    only: &
    dgedotf,                    &
    dlaset
use maria_kinds_mod,      only: &
    WP => DP
use maria_constants_mod,  only: &
    ZERO => D_ZERO
use maria_comparison_mod, only: &
    safe_eq
use maria_utils_mod,      only: &
    are_close
implicit none (type, external)

    call test_sanity()
    call test_quick()
    call test()

contains
    subroutine test_sanity()
        integer  :: info
        real(WP) :: foo(1), dot
        
        dot = dgedotf('a', -1, -1, foo, 0, foo, 0, info)
        ASSERT(info == -1)

        dot = dgedotf('n', -1, -1, foo, 0, foo, 0, info)
        ASSERT(info == -2)

        dot = dgedotf('n', 1, -1, foo, 0, foo, 0, info)
        ASSERT(info == -3)

        dot = dgedotf('n', 1, 1, foo, 0, foo, 0, info)
        ASSERT(info == -5)

        dot = dgedotf('n', 1, 1, foo, 1, foo, 0, info)
        ASSERT(info == -7)
    end subroutine test_sanity

    subroutine test_quick()
        integer  :: info
        real(WP) :: foo(1), dot
        
        dot = dgedotf('n', 0, 1, foo, 1, foo, 1, info)
        ASSERT(info == 0)
        ASSERT(safe_eq(dot, ZERO))

        dot = dgedotf('n', 1, 0, foo, 1, foo, 1, info)
        ASSERT(info == 0)
        ASSERT(safe_eq(dot, ZERO))
    end subroutine test_quick

    subroutine test()
        character(1)          :: transA
        integer               :: info, m, n, ldA, ldB, i
        real(WP)              :: dot1, dot2
        real(WP), allocatable :: A(:), B(:)

    !-- transA = 'n' -----------------------------------------------------------
        transA = 'n'
        m = 2
        n = 3
        ldA = m + 1
        ldB = m + 2
        allocate(A(ldA*n), B(ldb*n))
        do i = 1, lda*n
            A(i) = real(i * (-1)**i, WP)
        end do
        do i = 1, ldb*n
            B(i) = real(i**2, WP)
        end do
        dot1 = dgedotf(transA, m, n, A, lda, B, ldb, info)
        ASSERT(info == 0)
        dot2 = real(-1 + 2*2**2 + 4*5**2 - 5*6**2 - 7*9**2 + 8*10**2, WP)
        ASSERT(are_close(dot1, dot2, info, rtol=1e-12_WP))
        deallocate(A, B)

    !-- transA = 't' -----------------------------------------------------------
        transA = 't'
        m = 2
        n = 3
        ldA = n + 2
        ldB = m + 2
        allocate(A(ldA*m), B(ldb*n))
        do i = 1, lda*m
            A(i) = real(i * (-1)**i, WP)
        end do
        do i = 1, ldb*n
            B(i) = real(i, WP)
        end do
        dot1 = dgedotf(transA, m, n, A, lda, B, ldb, info)
        ASSERT(info == 0)
        dot2 = real(-1 + 6*2 + 2*5 - 7*6 - 3*9 + 8*10, WP)
        ASSERT(are_close(dot1, dot2, info, rtol=1e-12_WP))
        deallocate(A, B)

    !-- transA = 'n', ldA == ldB == m  -----------------------------------------
        transA = 'n'
        m = 10
        n = 20
        ldA = m + 2
        ldB = ldA
        allocate(A(ldA*n), B(ldb*n))
        do i = 1, lda*n
            A(i) = real(i * (-1)**i, WP)
        end do
        do i = 1, ldb*n
            B(i) = sin(real(i, WP))
        end do
        call dlaset('a', lda-m, n, ZERO, ZERO, A(m+1:), ldA)
        dot1 = dgedotf(transA, m, n, A, ldA, B, ldB, info)
        dot2 = dgedotf(transA, ldA, n, A, ldA, B, ldB, info)
        ASSERT(are_close(dot1, dot2, info, rtol=1e-12_WP))
    end subroutine test
end program dgedotf_test
