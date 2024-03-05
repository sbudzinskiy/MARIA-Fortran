program ipermute_test
#include "maria_assert_mod.h"
#ifdef MARIA_MKL
use maria_prng_mkl_mod,     only: &
    prng => prng_mkl
#else
use maria_prng_builtin_mod, only: &
    prng => prng_builtin
#endif
use maria_utils_mod,        only: &
    ipermute,                     &
    arange
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

        integer    :: info, foo(1)
        type(prng) :: rng

        call rng%init(seed, info)
        ASSERT(info == 0)

        call ipermute(rng, -1, foo, info)
        ASSERT(info == -2)

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test_sanity

    subroutine test_quick(seed)
        integer, intent(in) :: seed

        integer    :: info, foo(1)
        type(prng) :: rng

        call rng%init(seed, info)
        ASSERT(info == 0)

        call ipermute(rng, 0, foo, info)
        ASSERT(info == 0)

        call ipermute(rng, 1, foo, info)
        ASSERT(info == 0)

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer              :: info, n, i
        type(prng)           :: rng       
        integer, allocatable :: x(:), y(:)
 
        call rng%init(seed, info)
        ASSERT(info == 0)

        n = 10
        allocate(x(n + 1), y(n))
        call arange(n, x, 1, 1, info)
        ASSERT(info == 0)
        x(n + 1) = 0
        call arange(n, y, 1, 1, info)
        ASSERT(info == 0)

        call ipermute(rng, n, x, info)
        ASSERT(info == 0)
        ASSERT(any(x(1:n) /= y))
        do i = 1, n
            ASSERT(any(x(1:n) == i))
        end do
        ASSERT(x(n + 1) == 0)

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test
end program ipermute_test
