program imedian_of_3_test
#include "maria_assert_mod.h"
use maria_utils_mod, only: &
    imedian_of_3
implicit none (type, external)

    call test()

contains
    subroutine test()
        ASSERT(imedian_of_3(1, 2, 3) == 2)

        ASSERT(imedian_of_3(1, 3, 2) == 2)

        ASSERT(imedian_of_3(2, 1, 3) == 2)

        ASSERT(imedian_of_3(2, 3, 1) == 2)

        ASSERT(imedian_of_3(3, 1, 2) == 2)

        ASSERT(imedian_of_3(3, 2, 1) == 2)
    end subroutine test
end program imedian_of_3_test
