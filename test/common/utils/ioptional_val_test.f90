program ioptional_val_test
#include "maria_assert_mod.h"
use maria_utils_mod, only: &
    ioptional_val
implicit none (type, external)

    call test()

contains
    subroutine test()
        ASSERT(foo(1, 2) == 2)

        ASSERT(foo(1) == 1)
    end subroutine test

    integer function foo &
    (a, b)
        integer, intent(in)           :: a
        integer, intent(in), optional :: b

        foo = ioptional_val(a, b)
    end function foo
end program ioptional_val_test
