program dgesdd_q_test
#include "maria_assert_mod.h"
use maria_la_core_mod, only: &
    dgesdd_q
use maria_kinds_mod,   only: &
    WP => DP
implicit none (type, external)

    call test_sanity()
    call test_wrkwrk()

contains
    subroutine test_sanity()
        integer  :: info, m, n
        integer  :: ifoo(1)
        real(WP) :: foo(1)
        
        m = -1
        n = -1
        call dgesdd_q('x', m, n, foo, 0, foo, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -1)
 
        call dgesdd_q('s', m, n, foo, 0, foo, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -2)
 
        m = 2
        call dgesdd_q('s', m, n, foo, 0, foo, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -3)
 
        n = 3
        call dgesdd_q('s', m, n, foo, 0, foo, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -5)
 
        call dgesdd_q('s', m, n, foo, m, foo, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -8)
 
        call dgesdd_q('s', m, n, foo, m, foo, foo, m, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -10)
    end subroutine test_sanity

    subroutine test_wrkwrk()
        integer  :: info, m, n, ldA, ldU, ldVT, ifoo(1), lwork, liwork
        real(WP) :: foo(1)

        m = 4
        n = 5
        ldA = m + 1
        ldU = m + 1
        ldVT = n + 1

        call dgesdd_q('a', m, n, foo, ldA, foo, foo, ldU, foo, ldVT, foo, -1, ifoo, 0, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)
        ASSERT(liwork == 8 * min(m,n))
        call dgesdd_q('a', m, n, foo, ldA, foo, foo, ldU, foo, ldVT, foo, lwork-1, ifoo, liwork-1, info)
        ASSERT(info == -12)
        call dgesdd_q('a', m, n, foo, ldA, foo, foo, ldU, foo, ldVT, foo, lwork, ifoo, liwork-1, info)
        ASSERT(info == -14) 
    end subroutine test_wrkwrk
end program dgesdd_q_test
