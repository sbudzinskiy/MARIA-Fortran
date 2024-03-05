program dsafe_less_test
#include "maria_assert_mod.h"
use maria_comparison_mod, only: &
    dsafe_less
use maria_constants_mod,  only: &
    ONE => D_ONE,               &
    EPS => D_MACHTOL,           &
    HALF => D_HALF
implicit none (type, external)

    call test()

contains
    subroutine test()
        ASSERT(.not. dsafe_less(ONE, ONE - EPS))

        ASSERT(.not. dsafe_less(ONE, ONE + HALF * EPS))

        ASSERT(dsafe_less(ONE, ONE + EPS))
    end subroutine test
end program dsafe_less_test
