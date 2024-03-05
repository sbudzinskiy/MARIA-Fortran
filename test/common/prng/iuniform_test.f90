program iuniform_test
#include "maria_assert_mod.h"
#ifdef MARIA_MKL
use maria_prng_mkl_mod,     only: &
    prng => prng_mkl
#else
use maria_prng_builtin_mod, only: &
    prng => prng_builtin
#endif
implicit none (type, external)

    character(100) :: seed_str
    integer        :: seed

    call get_command_argument(1, seed_str)
    read(seed_str, *) seed

    call test_sanity(seed)
    call test_quick(seed)
    call test(seed)

contains
    subroutine test_sanity(seed)
        integer, intent(in) :: seed

        integer    :: info
        type(prng) :: rng
        integer    :: foo(1)

        call rng%init(seed, info)
        ASSERT(info == 0)

        call rng%iuniform(-1, foo, 0, -1, info)
        ASSERT(info == -2)
 
        call rng%iuniform(10, foo, 0, -1, info)
        ASSERT(info == -5)

        call rng%iuniform(10, foo, 0, 0, info)
        ASSERT(info == -5)

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test_sanity

    subroutine test_quick(seed)
        integer, intent(in) :: seed

        integer    :: info, foo(1)
        type(prng) :: rng

        call rng%init(seed, info)
        ASSERT(info == 0)

        call rng%iuniform(0, foo, 0, 1, info)
        ASSERT(info == 0)

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer              :: info, n, i
        integer              :: a, b
        integer, allocatable :: x(:)
        type(prng)           :: rng

        call rng%init(seed, info)
        ASSERT(info == 0)

        n = 1000
        a = -1
        b = 3

        allocate(x(n + 1))
        x = b + 1

        call rng%iuniform(n, x, a, b, info)
        ASSERT(info == 0)

    !-- Check range ------------------------------------------------------------
        ASSERT(x(n + 1) == b + 1)

    !-- Check coverage ---------------------------------------------------------
        ASSERT(all(x(1:n) >= a .and. x(1:n) <= b))
        do i = a, b
            ASSERT(any(x(1:n) == i))
        end do

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test   
end program iuniform_test
