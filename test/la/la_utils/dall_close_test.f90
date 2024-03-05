program dall_close_test
#include "maria_assert_mod.h"
use maria_la_utils_mod,  only: &
    dall_close
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

        ASSERT(.not. dall_close(-1, foo, 0, foo, 0, info))
        ASSERT(info == -1)

        ASSERT(.not. dall_close(1, foo, 0, foo, 0, info, atol=-EPS))
        ASSERT(info == -7)

        ASSERT(.not. dall_close(1, foo, 0, foo, 0, info, rtol=-EPS))
        ASSERT(info == -8)
    end subroutine test_sanity

    subroutine test_quick()
        integer  :: info
        real(WP) :: foo(1)

        ASSERT(dall_close(0, foo, 0, foo, 0, info))
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test()
        integer               :: n, incx, incy, info, lenx, leny
        real(WP), allocatable :: x(:), y(:)

        n = 5

    !-- Case 1 -----------------------------------------------------------------
        incx = 2
        incy = -3
        lenx = 1 + (n-1) * abs(incx)
        leny = 1 + (n-1) * abs(incy)
        allocate(x(lenx), y(leny))

        x = HALF
        y = HALF
        ASSERT(dall_close(n, x, incx, y, incy, info))
        ASSERT(info == 0)

        x(1 + incx - 1) = 2*HALF
        ASSERT(dall_close(n, x, incx, y, incy, info))
        ASSERT(info == 0)

        x(1 + incx) = x(1 + incx) + 2*EPS
        ASSERT(.not. dall_close(n, x, incx, y, incy, info))
        ASSERT(info == 0)

        ASSERT(dall_close(n, x, incx, y, incy, info, atol=3*EPS))
        ASSERT(info == 0)

        ASSERT(dall_close(n, x, incx, y, incy, info, rtol=6*EPS))
        ASSERT(info == 0)

        deallocate(x, y)

    !-- Case 2 -----------------------------------------------------------------
        incx = -4
        incy = 5
        lenx = 1 + (n-1) * abs(incx)
        leny = 1 + (n-1) * abs(incy)
        allocate(x(lenx), y(leny))

        x = HALF
        y = HALF
        ASSERT(dall_close(n, x, incx, y, incy, info))
        ASSERT(info == 0)

        y(1 + incy - 1) = 2*HALF
        ASSERT(dall_close(n, x, incx, y, incy, info))
        ASSERT(info == 0)

        y(1 + incy) = y(1 + incy) + 2*EPS
        ASSERT(.not. dall_close(n, x, incx, y, incy, info))
        ASSERT(info == 0)

        ASSERT(dall_close(n, x, incx, y, incy, info, atol=3*EPS))
        ASSERT(info == 0)

        ASSERT(dall_close(n, x, incx, y, incy, info, rtol=6*EPS))
        ASSERT(info == 0)

        deallocate(x, y)

    !-- Case 3 -----------------------------------------------------------------
        incx = 0
        incy = 0
        allocate(x(2), y(2))

        x = HALF
        y = HALF
        ASSERT(dall_close(n, x, incx, y, incy, info))
        ASSERT(info == 0)

        y(2) = 2*HALF
        ASSERT(dall_close(n, x, incx, y, incy, info))
        ASSERT(info == 0)

        y(1) = y(1) + 2*EPS
        ASSERT(.not. dall_close(n, x, incx, y, incy, info))
        ASSERT(info == 0)

        ASSERT(dall_close(n, x, incx, y, incy, info, atol=3*EPS))
        ASSERT(info == 0)

        ASSERT(dall_close(n, x, incx, y, incy, info, rtol=6*EPS))
        ASSERT(info == 0)
    end subroutine test
end program dall_close_test
