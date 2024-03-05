program isort_test
#include "maria_assert_mod.h"
#ifdef MARIA_MKL
use maria_prng_mkl_mod,     only: &
    prng => prng_mkl
#else
use maria_prng_builtin_mod, only: &
    prng => prng_builtin
#endif
use maria_utils_mod,        only: &
    isort,                        &
    arange,                       &
    permute
implicit none (type, external)

    character(100) :: seed_str
    integer        :: seed

    call get_command_argument(1, seed_str)
    read(seed_str, *) seed

    call test_sanity()
    call test_quick()
    call test(seed)

contains
    subroutine test_sanity()
        integer :: foo(1), info
        
        call isort('x', 0, foo, info)
        ASSERT(info == -1)

        call isort('i', -1, foo, info)
        ASSERT(info == -2)
    end subroutine test_sanity

    subroutine test_quick()
        integer :: foo(1), info
        
        call isort('i', 0, foo, info)
        ASSERT(info == 0)

        call isort('i', 1, foo, info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer              :: info, n, i
        type(prng)           :: rng       
        integer, allocatable :: x(:)
        
        call rng%init(seed, info)
        ASSERT(info == 0)

        n = 115
        allocate(x(n + 1))
        x(n + 1) = -1

        call arange(n, x, 1, 1, info)
        ASSERT(info == 0)

        call permute(rng, n, x, info)
        ASSERT(info == 0)
        call isort('i', n, x, info)
        ASSERT(info == 0)
        do i = 1, n
            ASSERT(x(i) == i)
        end do
        ASSERT(x(n + 1) == -1)

        call permute(rng, n, x, info)
        ASSERT(info == 0)
        call isort('d', n, x, info)
        ASSERT(info == 0)
        do i = 1, n
            ASSERT(x(i) == n - i + 1)
        end do
        ASSERT(x(n + 1) == -1)

        x = 1
        x(1) = 0
        x(2) = 2
        call permute(rng, n, x, info)
        ASSERT(info == 0)
        call isort('i', n, x, info)
        ASSERT(info == 0)
        ASSERT(x(1) == 0)
        ASSERT(x(n) == 2)

        call permute(rng, n, x, info)
        ASSERT(info == 0)
        call isort('d', n, x, info)
        ASSERT(info == 0)
        ASSERT(x(n) == 0)
        ASSERT(x(1) == 2)

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test
end program isort_test
