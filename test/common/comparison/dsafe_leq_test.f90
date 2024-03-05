program dsafe_leq_test
#include "maria_assert_mod.h"
use maria_comparison_mod, only: &
    dsafe_leq
use maria_constants_mod,  only: &
    ONE => D_ONE,               &
    EPS => D_MACHTOL,           &
    HALF => D_HALF
implicit none (type, external)

    call test()

contains
    subroutine test()
        ASSERT(.not. dsafe_leq(ONE, ONE - EPS))

        ASSERT(dsafe_leq(ONE, ONE - HALF * EPS))

        ASSERT(dsafe_leq(ONE, ONE + EPS))
    end subroutine test
end program dsafe_leq_test
