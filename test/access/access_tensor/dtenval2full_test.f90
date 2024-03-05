program dtenval2full_test
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
use maria_access_tensor_mod, only:  &
    i2mi,                           &
    dtenval,                        &
    dtenval2full,                   &
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

        call dtenval2full(fun, 1, n, foo, ifoo, 0, ierr)
        ASSERT( ierr == -2 )

        n = 2
        n(1) = -1
        call dtenval2full(fun, 3, n, foo, ifoo, 0, ierr)
        ASSERT( ierr == -3 )
    end subroutine test_sanity

    subroutine test()
        integer                     :: i, d, numel, ierr, liwork, ifoo(1)
        real(WP)                    :: foo(1)
        procedure(dtenval), pointer :: fun
        integer,    allocatable     :: n(:), mi(:), iwork(:)
        real(WP),   allocatable     :: A(:)

        fun => dtenval_hilbert

        d = 5
        
        allocate( n(d), mi(d) )
        n(1) = 5
        n(2) = 7
        n(3) = 4
        n(4) = 13
        n(5) = 3
        numel = product(n(1:d))

        call dtenval2full(fun, d, n, foo, ifoo, -1, ierr)
        liwork = ifoo(1)

        allocate( A(numel), iwork(liwork) )

        call dtenval2full(fun, d, n, A, iwork, liwork, ierr)
        do i = 1, numel
            call i2mi(d, n, i, mi, ierr)
            ASSERT( safe_eq(A(i), ONE / (sum(mi(1:d)) + d - 1)) )
        end do
    end subroutine test
end program dtenval2full_test
