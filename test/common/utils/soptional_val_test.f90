program soptional_val_test
#include "maria_assert_mod.h"
use maria_utils_mod,      only: &
    soptional_val
use maria_kinds_mod,      only: &
    WP => SP
use maria_constants_mod,  only: &
    ZERO => S_ZERO,             &
    ONE => S_ONE
use maria_comparison_mod, only: &
    safe_eq
implicit none (type, external)

    call test()

contains
    subroutine test()
        ASSERT(safe_eq(foo(ZERO, ONE), ONE))

        ASSERT(safe_eq(foo(ZERO), ZERO))
    end subroutine test

    real(WP) function foo &
    (a, b)
        real(WP), intent(in)           :: a
        real(WP), intent(in), optional :: b

        foo = soptional_val(a, b)
    end function foo
end program soptional_val_test
