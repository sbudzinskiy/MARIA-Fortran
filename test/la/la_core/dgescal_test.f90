program dgescal_test
#include "maria_assert_mod.h"
use maria_la_core_mod,   only: &
    dgescal
use maria_kinds_mod,     only: &
    WP => DP
use maria_constants_mod, only: &
    ONE => D_ONE, &
    TWO => D_TWO
use maria_la_utils_mod,  only: &
    geall_const
implicit none (type, external)

    call test_sanity()
    call test_quick()
    call test()

contains
    subroutine test_sanity()
        integer  :: info
        real(WP) :: foo(1)
        
        call dgescal(-1, -1, foo(1), foo, 0, info)
        ASSERT(info == -1)

        call dgescal(1, -1, foo(1), foo, 0, info)
        ASSERT(info == -2)

        call dgescal(1, 1, foo(1), foo, 0, info)
        ASSERT(info == -5)
    end subroutine test_sanity

    subroutine test_quick()
        integer  :: info
        real(WP) :: foo(1)
        
        call dgescal(1, 1, ONE, foo, 1, info)
        ASSERT(info == 0)

        call dgescal(0, 1, TWO, foo, 1, info)
        ASSERT(info == 0)

        call dgescal(1, 0, TWO, foo, 1, info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test()
        integer               :: info, m, n, ldA
        real(WP)              :: alpha
        real(WP), allocatable :: A(:)

        m = 3
        n = 2
        ldA = m + 1
        allocate(A(ldA*n))
        A = ONE
        alpha = TWO

        call dgescal(m, n, alpha, A, ldA, info)
        ASSERT(info == 0)
        ASSERT(geall_const(m, n, A, ldA, alpha, alpha, info))
        ASSERT(.not. geall_const(lda, n, A, ldA, alpha, alpha, info))

        A = ONE
        call dgescal(ldA, n, alpha, A, ldA, info)
        ASSERT(info == 0)
        ASSERT(geall_const(ldA, n, A, ldA, alpha, alpha, info))
    end subroutine test
end program dgescal_test
