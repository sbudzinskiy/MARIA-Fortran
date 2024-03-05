program dall_const_test
#include "maria_assert_mod.h"
use maria_la_utils_mod,  only: &
    dall_const
use maria_kinds_mod,     only: &
    WP => DP
use maria_constants_mod, only: &
    HALF => D_HALF,            &
    EPS => D_MACHTOL
implicit none (type, external)

    call test_sanity()
    call test_quick()
    call test()

contains
    subroutine test_sanity()
        integer  :: info
        real(WP) :: foo(1)

        ASSERT(.not. dall_const(-1, foo, 0, foo(1), info))
        ASSERT(info == -1)

        ASSERT(.not. dall_const(1, foo, 0, foo(1), info, atol=-EPS))
        ASSERT(info == -6)

        ASSERT(.not. dall_const(1, foo, 0, foo(1), info, rtol=-EPS))
        ASSERT(info == -7)
    end subroutine test_sanity

    subroutine test_quick()
        integer  :: info
        real(WP) :: foo(1)

        ASSERT(dall_const(0, foo, 0, foo(1), info))
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test()
        integer               :: n, incx, info, lenx
        real(WP), allocatable :: x(:), a

        n = 5

    !-- Case 1 -----------------------------------------------------------------
        incx = 2
        lenx = 1 + (n-1) * abs(incx)
        allocate(x(lenx))

        x = HALF
        a = HALF
        ASSERT(dall_const(n, x, incx, a, info))
        ASSERT(info == 0)

        x(1 + incx - 1) = 2*HALF
        ASSERT(dall_const(n, x, incx, a, info))
        ASSERT(info == 0)

        x(1 + incx) = x(1 + incx) + 2*EPS
        ASSERT(.not. dall_const(n, x, incx, a, info))
        ASSERT(info == 0)

        ASSERT(dall_const(n, x, incx, a, info, atol=3*EPS))
        ASSERT(info == 0)

        ASSERT(dall_const(n, x, incx, a, info, rtol=6*EPS))
        ASSERT(info == 0)

        deallocate(x)

    !-- Case 2 -----------------------------------------------------------------
        incx = -4
        lenx = 1 + (n-1) * abs(incx)
        allocate(x(lenx))

        x = HALF
        a = HALF
        ASSERT(dall_const(n, x, incx, a, info))
        ASSERT(info == 0)

        x(lenx + incx - 1) = 2*HALF
        ASSERT(dall_const(n, x, incx, a, info))
        ASSERT(info == 0)

        x(lenx + incx) = x(lenx + incx) + 2*EPS
        ASSERT(.not. dall_const(n, x, incx, a, info))
        ASSERT(info == 0)

        ASSERT(dall_const(n, x, incx, a, info, atol=3*EPS))
        ASSERT(info == 0)

        ASSERT(dall_const(n, x, incx, a, info, rtol=6*EPS))
        ASSERT(info == 0)

        deallocate(x)

    !-- Case 3 -----------------------------------------------------------------
        incx = 0
        allocate(x(2))

        x = HALF
        a = HALF
        ASSERT(dall_const(n, x, incx, a, info))
        ASSERT(info == 0)

        x(2) = 2*HALF
        ASSERT(dall_const(n, x, incx, a, info))
        ASSERT(info == 0)

        x(1) = x(1) + 2*EPS
        ASSERT(.not. dall_const(n, x, incx, a, info))
        ASSERT(info == 0)

        ASSERT(dall_const(n, x, incx, a, info, atol=3*EPS))
        ASSERT(info == 0)

        ASSERT(dall_const(n, x, incx, a, info, rtol=6*EPS))
        ASSERT(info == 0)
    end subroutine test
end program dall_const_test
