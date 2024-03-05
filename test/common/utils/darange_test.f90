program darange_test
#include "maria_assert_mod.h"
use maria_utils_mod,      only: &
    darange
use maria_kinds_mod,      only: &
    WP => DP
use maria_constants_mod,  only: &
    ZERO => D_ZERO,             &
    HALF => D_HALF
use maria_comparison_mod, only: &
    safe_eq
implicit none (type, external)

    call test_sanity()
    call test_quick()
    call test()

contains
    subroutine test_sanity()
        integer  :: info
        real(WP) :: foo(1)
        
        call darange(-1, foo, HALF, ZERO, info)
        ASSERT(info == -1)
    end subroutine test_sanity

    subroutine test_quick()
        integer  :: info
        real(WP) :: foo(1)
        
        call darange(0, foo, HALF, ZERO, info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test()
        integer               :: n, info, i
        real(WP)              :: a, step
        real(WP), allocatable :: x(:)

        n = 10
        allocate(x(n + 1))
        x = ZERO

        a = HALF
        step = -3*HALF

        call darange(n, x, a, step, info)
        ASSERT(info == 0)
        do i = 1, n
            ASSERT(safe_eq(x(i), a + (i - 1) * step))
        end do
        ASSERT(safe_eq(x(n+1), ZERO))
    end subroutine test
end program darange_test
