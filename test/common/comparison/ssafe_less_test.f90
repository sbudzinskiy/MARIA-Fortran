program ssafe_less_test
#include "maria_assert_mod.h"
use maria_comparison_mod, only: &
    ssafe_less
use maria_constants_mod,  only: &
    ONE => S_ONE,               &
    EPS => S_MACHTOL,           &
    HALF => S_HALF
implicit none (type, external)

    call test()

contains
    subroutine test()
        ASSERT(.not. ssafe_less(ONE, ONE - EPS))

        ASSERT(.not. ssafe_less(ONE, ONE + HALF * EPS))

        ASSERT(ssafe_less(ONE, ONE + EPS))
    end subroutine test
end program ssafe_less_test
