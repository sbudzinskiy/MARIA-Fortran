program dmatval2slc_test
#include "maria_assert_mod.h"
#ifdef MARIA_MKL
use maria_prng_mkl_mod,      only: &
    prng => prng_mkl
#else
use maria_prng_builtin_mod,  only: &
    prng => prng_builtin
#endif
use maria_kinds_mod,         only: &
    WP => DP
use maria_constants_mod,     only: &
    ONE => D_ONE
use maria_comparison_mod,    only: &
    safe_eq
use maria_access_matrix_mod, only: &
    dmatval,                       &
    dmatval2slc,                   &
    dmatval_hilbert,               &
    dmatval_lotkin
implicit none (type, external)

    character(100) :: seed_str
    integer        :: seed

    call get_command_argument(1, seed_str)
    read(seed_str, *) seed

    call test_sanity()
    call test()

contains
    subroutine test_sanity()
        integer                     :: info
        real(WP)                    :: foo(1)
        procedure(dmatval), pointer :: fun

        call dmatval2slc(fun, 0, 0, 0, 0, foo, 0, info)
        ASSERT(info == -2)

        call dmatval2slc(fun, 1, 0, 0, 0, foo, 0, info)
        ASSERT(info == -3)

        call dmatval2slc(fun, 1, 2, 0, 0, foo, 0, info)
        ASSERT(info == -4)

        call dmatval2slc(fun, 1, 2, 3, 0, foo, 0, info)
        ASSERT(info == -4)

        call dmatval2slc(fun, 1, 2, 1, 0, foo, 0, info)
        ASSERT(info == -5)

        call dmatval2slc(fun, 1, 2, 1, 2, foo, 0, info)
        ASSERT(info == -5)

        call dmatval2slc(fun, 1, 2, 2, 0, foo, 0, info)
        ASSERT(info == -5)

        call dmatval2slc(fun, 1, 2, 2, 3, foo, 0, info)
        ASSERT(info == -5)
    end subroutine test_sanity

    subroutine test()
        integer                     :: m, n, mode, ind, incx, info, i
        procedure(dmatval), pointer :: fun
        real(WP), allocatable       :: x(:)

    !-- hilbert, incx > 0
        fun => dmatval_hilbert
        m = 4
        n = 5
        incx = 2
        allocate(x(max(m,n) * abs(incx)))
        mode = 1
        ind = 2
        call dmatval2slc(fun, m, n, mode, ind, x, incx, info)
        do i = 1, n
            ASSERT(safe_eq(x(1 + (i-1)*incx), ONE / (ind + i - 1)))
        end do
        mode = 2
        ind = 3
        call dmatval2slc(fun, m, n, mode, ind, x, incx, info)
        do i = 1, m
            ASSERT(safe_eq(x(1 + (i-1)*incx), ONE / (ind + i - 1)))
        end do
        deallocate(x)

    !-- lotkin, incx < 0
        fun => dmatval_lotkin
        m = 4
        n = 5
        incx = -2
        allocate(x(max(m,n) * abs(incx)))
        mode = 1
        ind = 1
        call dmatval2slc(fun, m, n, mode, ind, x, incx, info)
        do i = 1, n
            ASSERT(safe_eq(x(1 + (i-n)*incx), ONE))
        end do
        mode = 1
        ind = 3
        call dmatval2slc(fun, m, n, mode, ind, x, incx, info)
        do i = 1, n
            ASSERT(safe_eq(x(1 + (i-n)*incx), ONE / (ind + i - 1)))
        end do
        mode = 2
        ind = 3
        call dmatval2slc(fun, m, n, mode, ind, x, incx, info)
        ASSERT(safe_eq(x(1 + (1-m)*incx), ONE))
        do i = 2, m
            ASSERT(safe_eq(x(1 + (i-m)*incx), ONE / (ind + i - 1)))
        end do
    end subroutine test
end program dmatval2slc_test
