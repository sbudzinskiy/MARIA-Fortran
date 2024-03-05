program sarg_is_bad_test
#include "maria_assert_mod.h"
use maria_constants_mod, only: &
    ONE => S_ONE,              &
    EPS => S_MACHTOL
use maria_argcheck_mod,  only: &
    sarg_is_bad,               &
    BAD_IF_LESS,               &
    BAD_IF_MORE,               &
    BAD_IF_SAME
implicit none (type, external)

    call test()

contains
    subroutine test()
        ASSERT(sarg_is_bad(BAD_IF_LESS, ONE, ONE + EPS))

        ASSERT(.not. sarg_is_bad(BAD_IF_LESS, ONE, ONE - EPS))

        ASSERT(sarg_is_bad(BAD_IF_MORE, ONE, ONE - EPS))

        ASSERT(.not. sarg_is_bad(BAD_IF_MORE, ONE, ONE + EPS))

        ASSERT(sarg_is_bad(BAD_IF_SAME, ONE, ONE + EPS/2))

        ASSERT(sarg_is_bad(BAD_IF_SAME, ONE, ONE - EPS/2))

        ASSERT(.not. sarg_is_bad(BAD_IF_SAME, ONE, ONE + EPS))

        ASSERT(sarg_is_bad('q', ONE, ONE))
    end subroutine test
end program sarg_is_bad_test
