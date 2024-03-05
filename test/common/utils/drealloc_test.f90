program drealloc_test
#include "maria_assert_mod.h"
use maria_utils_mod,      only: &
    darange,                    &
    drealloc
use maria_kinds_mod,      only: &
    WP => DP
use maria_constants_mod,  only: &
    ONE => D_ONE
use maria_comparison_mod, only: &
    safe_eq
implicit none (type, external)

    call test_sanity()
    call test()

contains
    subroutine test_sanity()
        integer  :: ierr
        real(WP), allocatable :: foo(:)
        
        call drealloc(.true., foo, 0, ierr)
        ASSERT( ierr == -3 )
        ASSERT( .not. allocated(foo) )

        call drealloc(.true., foo, 1, ierr, 1)
        ASSERT( ierr == -5 )
        ASSERT( .not. allocated(foo) )
    end subroutine test_sanity

    subroutine test()
        integer               :: n, ierr, i
        real(WP), allocatable :: x(:)

        n = 10
        call drealloc(.true., x, n, ierr)
        ASSERT( allocated(x) )
        ASSERT( size(x) == n )

        call drealloc(.false., x, 2*n, ierr)
        ASSERT( allocated(x) )
        ASSERT( size(x) == 2*n )

        call darange(2*n, x, ONE, ONE, ierr)
        call drealloc(.true., x, 3*n, ierr, factor=2)
        ASSERT( allocated(x) )
        ASSERT( size(x) == 4*n )
        do i = 1, 2*n
            ASSERT( safe_eq(x(i), i*ONE) )
        end do
    end subroutine test
end program drealloc_test
