program sroundup_lwork_test
#include "maria_assert_mod.h"
use maria_la_core_mod, only: &
    sroundup_lwork
use maria_kinds_mod,   only: &
    WP => SP
implicit none (type, external)

    call test()

contains
    subroutine test()
        integer  :: n
        real(WP) :: x

        n = 16777217
        x = real(n, WP)
        ASSERT(nint(x) < n)
        x = sroundup_lwork(n)
        ASSERT(nint(x) >= n)
    end subroutine test
end program sroundup_lwork_test
