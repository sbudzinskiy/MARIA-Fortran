program mi2i_test
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
        integer :: ierr, n(2), mi(2), ifoo

        call mi2i(1, n, mi, ifoo, ierr)
        ASSERT(ierr == -1)

        n(1) = 0
        n(2) = 5
        call mi2i(2, n, mi, ifoo, ierr)
        ASSERT( ierr == -2 )

        n(1) = 3
        mi(1) = 0
        mi(2) = n(2) + 1
        call mi2i(2, n, mi, ifoo, ierr)
        ASSERT( ierr == -3 )

        mi(1) = 1
        call mi2i(2, n, mi, ifoo, ierr)
        ASSERT( ierr == -3 )
    end subroutine test_sanity

    subroutine test(seed)
        integer,    intent(in) :: seed

        integer :: ierr, d, k, i
        type(prng)  :: rng
        integer,    allocatable :: n(:), mi(:), mimi(:)

        call rng%init(seed, ierr)
        ASSERT( ierr == 0 )

        d = 5
        allocate(n(d), mi(d), mimi(d))
        call rng%iuniform(d, n, 2, 10, ierr)
        do k = 1, d
            call rng%iuniform(1, mi(k:), 1, n(k), ierr)
        end do

        call mi2i(d, n, mi, i, ierr)
        call i2mi(d, n, i, mimi, ierr)
        do k = 1, d
            ASSERT( mi(k) == mimi(k) )
        end do
        deallocate(n, mi, mimi)

        d = 3
        allocate(n(d), mi(d))
        call rng%iuniform(d, n, 2, 10, ierr)
        do k = 1, d
            call rng%iuniform(1, mi(k:), 1, n(k), ierr)
        end do
        call mi2i(d, n, mi, i, ierr)
        i = i - mi(1) - (mi(2)-1)*n(1) - (mi(3)-1)*n(1)*n(2)
        ASSERT( i == 0 )
        deallocate(n, mi)
    end subroutine test
end program mi2i_test
