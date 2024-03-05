program smatval2full_test
#include "maria_assert_mod.h"
#ifdef MARIA_MKL
use maria_prng_mkl_mod,      only: &
    prng => prng_mkl
#else
use maria_prng_builtin_mod,  only: &
    prng => prng_builtin
#endif
use maria_kinds_mod,         only: &
    WP => SP
use maria_constants_mod,     only: &
    ONE => S_ONE
use maria_comparison_mod,    only: &
    safe_eq
use maria_access_matrix_mod, only: &
    smatval,                       &
    smatval2full,                  &
    smatval_lotkin
implicit none (type, external)

    character(100) :: seed_str
    integer        :: seed

    call get_command_argument(1, seed_str)
    read(seed_str, *) seed

    call test_sanity()
    call test_quick()
    call test()

contains
    subroutine test_sanity()
        integer                     :: info
        real(WP)                    :: foo(1)
        procedure(smatval), pointer :: fun

        call smatval2full(fun, 'x', -1, -1, foo, 0, info)
        ASSERT(info == -2)

        call smatval2full(fun, 'n', -1, -1, foo, 0, info)
        ASSERT(info == -3)

        call smatval2full(fun, 'n', 0, -1, foo, 0, info)
        ASSERT(info == -4)

        call smatval2full(fun, 'n', 0, 0, foo, 0, info)
        ASSERT(info == -6)

        call smatval2full(fun, 'n', 3, 2, foo, 2, info)
        ASSERT(info == -6)

        call smatval2full(fun, 't', 2, 3, foo, 2, info)
        ASSERT(info == -6)
    end subroutine test_sanity
 
    subroutine test_quick()
        integer                     :: info
        real(WP)                    :: foo(1)
        procedure(smatval), pointer :: fun

        call smatval2full(fun, 'n', 0, 1, foo, 1, info)
        ASSERT(info == 0)

        call smatval2full(fun, 'n', 1, 0, foo, 1, info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test()
        character(1)                :: trans
        integer                     :: m, n, lda, info, i, j
        procedure(smatval), pointer :: fun
        real(WP), allocatable       :: A(:)

        fun => smatval_lotkin

    !-- trans = 'n'
        trans = 'n'
        m = 4
        n = 5
        lda = m + 1
        allocate(A(lda*n))
        call smatval2full(fun, trans, m, n, A, lda, info)
        do j = 1, n
            ASSERT(safe_eq(A(1 + (j-1)*lda), ONE))
            do i = 2, m
                ASSERT(safe_eq(A(i + (j-1)*lda), ONE / (i + j - 1)))
            end do
        end do
        deallocate(A)

    !-- trans = 't'
        trans = 't'
        m = 4
        n = 5
        lda = n + 1
        allocate(A(lda*m))
        call smatval2full(fun, trans, m, n, A, lda, info)
        do j = 1, n
            ASSERT(safe_eq(A(j), ONE))
            do i = 2, m
                ASSERT(safe_eq(A(j + (i-1)*lda), ONE / (i + j - 1)))
            end do
        end do
    end subroutine test
end program smatval2full_test
