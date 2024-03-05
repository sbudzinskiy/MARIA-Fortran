program ttrank_test
#include "maria_assert_mod.h"
#ifdef MARIA_MKL
use maria_prng_mkl_mod,     only: &
    prng => prng_mkl
#else
use maria_prng_builtin_mod, only: &
    prng => prng_builtin
#endif
use maria_tt_utils_mod,        only: &
    ttrank
implicit none (type, external)

    character(100) :: seed_str
    integer        :: seed

    call get_command_argument(1, seed_str)
    read(seed_str, *) seed

    call test_sanity()
    call test(seed)

contains
    subroutine test_sanity()
        integer :: ifoo(1), ierr
        
        call ttrank(1, ifoo, 0, ifoo(1), ifoo(1), ierr)
        ASSERT(ierr == -1)

        call ttrank(2, ifoo, 0, ifoo(1), ifoo(1), ierr)
        ASSERT(ierr == -3)

        call ttrank(2, ifoo, 3, ifoo(1), ifoo(1), ierr)
        ASSERT(ierr == -3)
    end subroutine test_sanity

    subroutine test(seed)
        integer, intent(in) :: seed

        integer              :: d, k, rl, rr, ierr, r0, rd
        type(prng)           :: rng       
        integer, allocatable :: r(:)
        
        call rng%init(seed, ierr)
        ASSERT(ierr == 0)

        d = 5
        allocate(r(d-1))
        call rng%iuniform(d-1, r, 1, 10, ierr)

        k = 2
        call ttrank(d, r, k, rl, rr, ierr)
        ASSERT( ierr == 0 )
        ASSERT( rl == r(k-1) )
        ASSERT( rr == r(k) )

        k = 1
        call ttrank(d, r, k, rl, rr, ierr)
        ASSERT( ierr == 0 )
        ASSERT( rl == 1 )
        ASSERT( rr == r(k) )

        k = d
        call ttrank(d, r, k, rl, rr, ierr)
        ASSERT( ierr == 0 )
        ASSERT( rl == r(k-1) )
        ASSERT( rr == 1 )

        k = 1
        r0 = 2
        call ttrank(d, r, k, rl, rr, ierr, r0=r0)
        ASSERT( ierr == 0 )
        ASSERT( rl == r0 )
        ASSERT( rr == r(k) )

        k = d
        rd = 2
        call ttrank(d, r, k, rl, rr, ierr, rd=rd)
        ASSERT( ierr == 0 )
        ASSERT( rl == r(k-1) )
        ASSERT( rr == rd )
    end subroutine test
end program ttrank_test
