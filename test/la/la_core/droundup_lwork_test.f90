program droundup_lwork_test
#include "maria_assert_mod.h"
use maria_la_core_mod,             only: &
    droundup_lwork
use maria_kinds_mod,               only: &
    WP => DP
use, intrinsic :: iso_fortran_env, only: &
    INT64 
implicit none (type, external)

    call test()

contains
    subroutine test()
        integer  :: n
        real(WP) :: x
        
    !-- 32-bit integers are always rounded correctly        
        if (kind(n) == INT64) then
        !-- n = 9007199254740993
            n = 2**(8 * kind(n) - 11) + 1
            x = real(n, WP)
            ASSERT(nint(x) < n)
            x = droundup_lwork(n)
            ASSERT(nint(x) >= n)
        end if
    end subroutine test
end program droundup_lwork_test
