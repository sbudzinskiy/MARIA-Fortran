program ssafe_eq_test
#include "maria_assert_mod.h"
use maria_comparison_mod, only: &
    ssafe_eq
use maria_constants_mod,  only: &
    ONE => S_ONE,               &
    EPS => S_MACHTOL,           &
    HALF => S_HALF
implicit none (type, external)

    call test()

contains
    subroutine test()
        ASSERT(.not. ssafe_eq(ONE, ONE - EPS))

        ASSERT(ssafe_eq(ONE, ONE - HALF * EPS))

        ASSERT(ssafe_eq(ONE, ONE + HALF * EPS))

        ASSERT(.not. ssafe_eq(ONE, ONE + EPS))
    end subroutine test
end program ssafe_eq_test
