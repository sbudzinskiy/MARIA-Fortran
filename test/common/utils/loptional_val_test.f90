program loptional_val_test
#include "maria_assert_mod.h"
use maria_utils_mod, only: &
    loptional_val
implicit none (type, external)

    call test()

contains
    subroutine test()
        ASSERT(.not. foo(.true., .false.))

        ASSERT(foo(.true.))
    end subroutine test

    logical function foo &
    (a, b)
        logical, intent(in)           :: a
        logical, intent(in), optional :: b

        foo = loptional_val(a, b)
    end function foo
end program loptional_val_test
