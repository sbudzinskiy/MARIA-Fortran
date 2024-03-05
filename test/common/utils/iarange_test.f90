program iarange_test
#include "maria_assert_mod.h"
use maria_utils_mod, only: &
    iarange
implicit none (type, external)

    call test_sanity()
    call test_quick()
    call test()

contains
    subroutine test_sanity()
        integer :: info, foo(1)
        
        call iarange(-1, foo, 1, 1, info)
        ASSERT(info == -1)
    end subroutine test_sanity

    subroutine test_quick()
        integer :: info, foo(1)
        
        call iarange(0, foo, 1, 1, info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test()
        integer              :: n, info, i
        integer              :: a, step
        integer, allocatable :: x(:)

        n = 10
        allocate(x(n + 1))
        x = 3

        a = 2
        step = -3

        call iarange(n, x, a, step, info)
        ASSERT(info == 0)
        do i = 1, n
            ASSERT(x(i) == a + (i - 1) * step)
        end do
        ASSERT(x(n+1) == 3)
    end subroutine test
end program iarange_test
