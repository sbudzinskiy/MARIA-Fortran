program iarg_is_bad_test
#include "maria_assert_mod.h"
use maria_argcheck_mod, only: &
    iarg_is_bad,              &
    BAD_IF_LESS,              &
    BAD_IF_MORE,              &
    BAD_IF_SAME
implicit none (type, external)

    call test()

contains
    subroutine test()
        ASSERT(iarg_is_bad(BAD_IF_LESS, 1, 2))

        ASSERT(.not. iarg_is_bad(BAD_IF_LESS, 2, 1))

        ASSERT(iarg_is_bad(BAD_IF_MORE, 2, 1))

        ASSERT(.not. iarg_is_bad(BAD_IF_MORE, 1, 2))

        ASSERT(iarg_is_bad(BAD_IF_SAME, 1, 1))

        ASSERT(.not. iarg_is_bad(BAD_IF_SAME, 1, 2))

        ASSERT(iarg_is_bad('q', 1, 1))
    end subroutine test
end program iarg_is_bad_test
