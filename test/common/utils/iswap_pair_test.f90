program iswap_pair_test
#include "maria_assert_mod.h"
use maria_utils_mod, only: &
    iswap_pair
implicit none (type, external)

    call test()

contains
    subroutine test()
        integer  :: a, b

        a = 1
        b = 2
        call iswap_pair(a, b)
        ASSERT(a == 2)
        ASSERT(b == 1)
    end subroutine test
end program iswap_pair_test
