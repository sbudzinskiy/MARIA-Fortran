program i2mi_test
#include "maria_assert_mod.h"
#ifdef MARIA_MKL
use maria_prng_mkl_mod,      only: &
    prng => prng_mkl
#else
use maria_prng_builtin_mod,  only: &
    prng => prng_builtin
#endif
use maria_access_tensor_mod, only: &
    mi2i, i2mi
implicit none (type, external)

    character(100) :: seed_str
    integer        :: seed

    call get_command_argument(1, seed_str)
    read(seed_str, *) seed

    call test_sanity()
    call test(seed)

contains
    subroutine test_sanity()
        integer :: ierr, n(2), ifoo(1)

        call i2mi(1, n, 0, ifoo, ierr)
        ASSERT(ierr == -1)

        n(1) = 0
        n(2) = 5
        call i2mi(2, n, 0, ifoo, ierr)
        ASSERT( ierr == -2 )

        n(1) = 3
        call i2mi(2, n, 0, ifoo, ierr)
        ASSERT( ierr == -3 )

        call i2mi(2, n, n(1)*n(2)+1, ifoo, ierr)
        ASSERT( ierr == -3 )
    end subroutine test_sanity

    subroutine test(seed)
        integer,    intent(in) :: seed

        integer :: ierr, d, i, ifoo(1)
        type(prng)  :: rng
        integer,    allocatable :: n(:), mi(:)

        call rng%init(seed, ierr)
        ASSERT( ierr == 0 )

        d = 5
        allocate(n(d), mi(d))
        call rng%iuniform(d, n, 1, 10, ierr)
        call rng%iuniform(1, ifoo, 1, product(n), ierr)
        call i2mi(d, n, ifoo(1), mi, ierr)
        call mi2i(d, n, mi, i, ierr)
        ASSERT( ifoo(1) == i )
    end subroutine test
end program i2mi_test
