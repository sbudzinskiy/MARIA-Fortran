program dare_close_test
#include "maria_assert_mod.h"
use maria_utils_mod, only: &
    dare_close
use maria_constants_mod, only: &
    HALF => D_HALF,            &
    EPS => D_MACHTOL
implicit none (type, external)

    call test_sanity()
    call test()

contains
    subroutine test_sanity()
        integer :: info

        ASSERT(.not. dare_close(HALF, HALF, info, atol=-EPS))
        ASSERT(info == -4)

        ASSERT(.not. dare_close(HALF, HALF, info, rtol=-EPS))
        ASSERT(info == -5)
    end subroutine test_sanity

    subroutine test()
        integer :: info

        ASSERT(dare_close(HALF, HALF, info))
        ASSERT(info == 0)

        ASSERT(.not. dare_close(HALF, HALF + EPS, info))
        ASSERT(info == 0)

        ASSERT(.not. dare_close(HALF, HALF + EPS, info, atol=EPS))
        ASSERT(info == 0)

        ASSERT(dare_close(HALF, HALF + EPS, info, atol=2*EPS))
        ASSERT(info == 0)

        ASSERT(.not. dare_close(HALF, HALF + EPS, info, rtol=2*EPS))
        ASSERT(info == 0)

        ASSERT(dare_close(HALF, HALF + EPS, info, rtol=4*EPS))
        ASSERT(info == 0)
    end subroutine test
end program dare_close_test
