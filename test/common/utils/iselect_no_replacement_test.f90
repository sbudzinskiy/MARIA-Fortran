program iselect_no_replacement_test
#include "maria_assert_mod.h"
#ifdef MARIA_MKL
use maria_prng_mkl_mod,     only: &
    prng => prng_mkl
#else
use maria_prng_builtin_mod, only: &
    prng => prng_builtin
#endif
use maria_utils_mod,        only: &
    iselect_no_replacement,       &
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

        call iselect_no_replacement(rng, -1, foo, -1, info)
        ASSERT(info == -2)

        call iselect_no_replacement(rng, 10, foo, -1, info)
        ASSERT(info == -4)

        call iselect_no_replacement(rng, 10, foo, 10, info)
        ASSERT(info == -4)

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test_sanity

    subroutine test_quick(seed)
        integer, intent(in) :: seed

        integer    :: info, foo(1)
        type(prng) :: rng

        call rng%init(seed, info)
        ASSERT(info == 0)

        call iselect_no_replacement(rng, 0, foo, 0, info)
        ASSERT(info == 0)

        call iselect_no_replacement(rng, 1, foo, 0, info)
        ASSERT(info == 0)

        call iselect_no_replacement(rng, 2, foo, 0, info)
        ASSERT(info == 0)

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer              :: info, n, i, k
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

        k = 5
        call iselect_no_replacement(rng, n, x, k, info)
        ASSERT(info == 0)
        ASSERT(any(x(1:k) /= y(1:k)))
        do i = 1, n
            ASSERT(any(x(1:n) == i))
        end do
        ASSERT(x(n + 1) == 0)

        k = n - 1
        call iselect_no_replacement(rng, n, x, k, info)
        ASSERT(info == 0)
        ASSERT(any(x(1:k) /= y(1:k)))
        do i = 1, n
            ASSERT(any(x(1:n) == i))
        end do
        ASSERT(x(n + 1) == 0)

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test
end program iselect_no_replacement_test
