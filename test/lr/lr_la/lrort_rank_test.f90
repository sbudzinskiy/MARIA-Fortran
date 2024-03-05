program lrort_rank_test
#include "maria_assert_mod.h"
use maria_lr_la_mod, only: &
    lrort_rank
implicit none (type, external)

    call test_sanity()
    call test()

contains
    subroutine test_sanity()
        integer :: info
        integer :: newr

        call lrort_rank('x', -1, -1, -1, newr, info)
        ASSERT(info == -1)

        call lrort_rank('l', -1, -1, -1, newr, info)
        ASSERT(info == -2)

        call lrort_rank('l', 1, -1, -1, newr, info)
        ASSERT(info == -3)

        call lrort_rank('l', 1, 1, -1, newr, info)
        ASSERT(info == -4)
    end subroutine test_sanity

    subroutine test()
        integer :: m, n, r, newr, info

        m = 20
        n = 10
        r = 15
        call lrort_rank('l', m, n, r, newr, info)
        ASSERT(info == 0)
        ASSERT(newr == min(m,r))
        call lrort_rank('r', m, n, r, newr, info)
        ASSERT(info == 0)
        ASSERT(newr == min(n,r))
    end subroutine test
end program lrort_rank_test
