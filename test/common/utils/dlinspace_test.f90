program dlinspace_test
#include "maria_assert_mod.h"
use maria_utils_mod,      only: &
    dlinspace
use maria_kinds_mod,      only: &
    WP => DP
use maria_constants_mod,  only: &
    ZERO => D_ZERO,             &
    HALF => D_HALF,             &
    EPS => D_MACHTOL
use maria_comparison_mod, only: &
    safe_eq
implicit none (type, external)

    call test_sanity()
    call test_quick()
    call test()

contains
    subroutine test_sanity()
        integer  :: info
        real(WP) :: foo(1), step
        
        call dlinspace(-1, foo, HALF, ZERO, info, step=step)
        ASSERT(info == -1)

        call dlinspace(2, foo, HALF, HALF - EPS, info)
        ASSERT(info == -4)

        call dlinspace(2, foo, HALF, HALF, info)
        ASSERT(info == -4)
    end subroutine test_sanity

    subroutine test_quick()
        integer  :: info
        real(WP) :: a, b, step, x(1)
        
        a = -HALF
        b = HALF
        call dlinspace(0, x, a, b, info, step=step)
        ASSERT(info == 0)
        ASSERT(safe_eq(step, ZERO))

        call dlinspace(1, x, a, b, info, step=step)
        ASSERT(info == 0)
        ASSERT(safe_eq(x(1), a))
        ASSERT(safe_eq(step, ZERO))

        call dlinspace(1, x, a, b, info, include_b=.false., step=step)
        ASSERT(info == 0)
        ASSERT(safe_eq(x(1), a))
        ASSERT(safe_eq(step, ZERO))
    end subroutine test_quick

    subroutine test()
        integer               :: n, info, i
        real(WP)              :: a, b, step
        real(WP), allocatable :: x(:)

        n = 10
        allocate(x(n + 1))
        x = ZERO

        a = HALF
        b = 3*HALF

        call dlinspace(n, x, a, b, info, step=step)
        ASSERT(info == 0)
        ASSERT(safe_eq(step, (b - a) / (n - 1)))
        do i = 1, n
            ASSERT(safe_eq(x(i), a + (i - 1) * step))
        end do
        ASSERT(safe_eq(x(n+1), ZERO))

        call dlinspace(n, x, a, b, info, include_b=.false., step=step)
        ASSERT(info == 0)
        ASSERT(safe_eq(step, (b - a) / n))
        do i = 1, n
            ASSERT(safe_eq(x(i), a + (i - 1) * step))
        end do
        ASSERT(safe_eq(x(n+1), ZERO))
   end subroutine test
end program dlinspace_test
