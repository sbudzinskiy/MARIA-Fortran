program sswap_pair_test
#include "maria_assert_mod.h"
use maria_utils_mod,      only: &
    sswap_pair
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
        real(WP) :: a, b

        a = ZERO
        b = ONE
        call sswap_pair(a, b)
        ASSERT(safe_eq(a, ONE))
        ASSERT(safe_eq(b, ZERO))
    end subroutine test
end program sswap_pair_test
