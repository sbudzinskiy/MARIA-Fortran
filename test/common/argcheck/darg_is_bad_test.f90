program darg_is_bad_test
#include "maria_assert_mod.h"
use maria_constants_mod, only: &
    ONE => D_ONE,              &
    EPS => D_MACHTOL
use maria_argcheck_mod,  only: &
    darg_is_bad,               &
    BAD_IF_LESS,               &
    BAD_IF_MORE,               &
    BAD_IF_SAME
implicit none (type, external)

    call test()

contains
    subroutine test()
        ASSERT(darg_is_bad(BAD_IF_LESS, ONE, ONE + EPS))

        ASSERT(.not. darg_is_bad(BAD_IF_LESS, ONE, ONE - EPS))

        ASSERT(darg_is_bad(BAD_IF_MORE, ONE, ONE - EPS))

        ASSERT(.not. darg_is_bad(BAD_IF_MORE, ONE, ONE + EPS))

        ASSERT(darg_is_bad(BAD_IF_SAME, ONE, ONE + EPS/2))

        ASSERT(darg_is_bad(BAD_IF_SAME, ONE, ONE - EPS/2))

        ASSERT(.not. darg_is_bad(BAD_IF_SAME, ONE, ONE + EPS))

        ASSERT(darg_is_bad('q', ONE, ONE))
    end subroutine test
end program darg_is_bad_test
