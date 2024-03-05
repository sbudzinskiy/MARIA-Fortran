program dtenval2fib_test
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
use maria_access_tensor_mod, only: &
    dtenval,                       &
    dtenval2fib,                   &
    dtenval_hilbert
implicit none (type, external)

    character(100) :: seed_str
    integer        :: seed

    call get_command_argument(1, seed_str)
    read(seed_str, *) seed

    call test_sanity()
    call test()

contains
    subroutine test_sanity()
        integer                     :: n(3), ierr, ifoo(1)
        real(WP)                    :: foo(1)
        procedure(dtenval), pointer :: fun

        call dtenval2fib(fun, 1, n, 0, ifoo, ifoo, foo, 0, ifoo, 0, ierr)
        ASSERT( ierr == -2 )

        n = 2
        n(1) = -1
        call dtenval2fib(fun, 3, n, 0, ifoo, ifoo, foo, 0, ifoo, 0, ierr)
        ASSERT( ierr == -3 )

        n(1) = 2
        call dtenval2fib(fun, 3, n, 0, ifoo, ifoo, foo, 0, ifoo, 0, ierr)
        ASSERT( ierr == -4 )

        call dtenval2fib(fun, 3, n, 4, ifoo, ifoo, foo, 0, ifoo, 0, ierr)
        ASSERT( ierr == -4 )
    end subroutine test_sanity

    subroutine test()
        integer                     :: i, d, k, incx, ierr, liwork, ifoo(1)
        real(WP)                    :: foo(1)
        procedure(dtenval), pointer :: fun
        integer,    allocatable     :: n(:), li(:), ri(:), iwork(:)
        real(WP),   allocatable     :: x(:)

        fun => dtenval_hilbert

        d = 5
        
        allocate( n(d), li(d), ri(d) )
        n(1) = 5
        n(2) = 7
        n(3) = 4
        n(4) = 13
        n(5) = 3

    !-- k = 1, incx > 0
        k = 1
        incx = 2

        call dtenval2fib(fun, d, n, k, ifoo, ifoo, foo, incx, ifoo, -1, ierr)
        liwork = ifoo(1)

        allocate( x(n(k)*abs(incx)), iwork(liwork) )

        ri(1) = 4
        ri(2) = 2
        ri(3) = 5
        ri(4) = 3

        call dtenval2fib(fun, d, n, k, li, ri, x, incx, iwork, liwork, ierr)
        ASSERT( ierr == 0 )

        do i = 1, n(k)
            ASSERT( safe_eq(x(1 + (i-1)*incx), ONE / (sum(li(1:k-1)) + i + sum(ri(1:d-k)) + d - 1)) )
        end do

        deallocate( x, iwork )

    !-- k = d, incx < 0
        k = d
        incx = -3

        call dtenval2fib(fun, d, n, k, ifoo, ifoo, foo, incx, ifoo, -1, ierr)
        liwork = ifoo(1)

        allocate( x(n(k)*abs(incx)), iwork(liwork) )

        li(1) = 3
        li(2) = 6
        li(3) = 2
        li(4) = 8

        call dtenval2fib(fun, d, n, k, li, ri, x, incx, iwork, liwork, ierr)
        ASSERT( ierr == 0 )

        do i = 1, n(k)
            ASSERT( safe_eq(x(1 + (i-n(k))*incx), ONE / (sum(li(1:k-1)) + i + sum(ri(1:d-k)) + d - 1)) )
        end do

        deallocate( x, iwork )

    !-- 1 < k < d, incx < 0
        k = 3
        incx = -2

        call dtenval2fib(fun, d, n, k, ifoo, ifoo, foo, incx, ifoo, -1, ierr)
        liwork = ifoo(1)

        allocate( x(n(k)*abs(incx)), iwork(liwork) )

        ri(1) = 9
        ri(2) = 2

        call dtenval2fib(fun, d, n, k, li, ri, x, incx, iwork, liwork, ierr)
        ASSERT( ierr == 0 )

        do i = 1, n(k)
            ASSERT( safe_eq(x(1 + (i-n(k))*incx), ONE / (sum(li(1:k-1)) + i + sum(ri(1:d-k)) + d - 1)) )
        end do
    end subroutine test
end program dtenval2fib_test
