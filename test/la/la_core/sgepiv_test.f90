program sgepiv_test
#include "maria_assert_mod.h"
#ifdef MARIA_MKL
use maria_prng_mkl_mod,     only: &
    prng => prng_mkl
#else
use maria_prng_builtin_mod, only: &
    prng => prng_builtin
#endif
use maria_la_core_mod,      only: &
    sgepiv,                       &
    slacpy
use maria_kinds_mod,        only: &
    WP => SP
use maria_constants_mod,    only: &
    ZERO => S_ZERO,               &
    ONE => S_ONE
use maria_utils_mod,        only: &
    arange
use maria_la_utils_mod,     only: &
    all_close
implicit none (type, external)

    character(100) :: seed_str
    integer        :: seed

    call get_command_argument(1, seed_str)
    read(seed_str, *) seed

    call test_sanity()
    call test(seed)
    call test_runtime_err()

contains
    subroutine test_sanity()
        integer  :: info, ifoo(1)
        real(WP) :: foo(1)
        
        call sgepiv('x', 'x', 0, 0, foo, -1, -1, -1, ifoo, info)
        ASSERT(info == -1)

        call sgepiv('r', 'x', 0, 0, foo, -1, -1, -1, ifoo, info)
        ASSERT(info == -2)

        call sgepiv('r', 'f', 0, 0, foo, -1, -1, -1, ifoo, info)
        ASSERT(info == -3)

        call sgepiv('r', 'f', 1, 0, foo, -1, -1, -1, ifoo, info)
        ASSERT(info == -4)

        call sgepiv('r', 'f', 1, 2, foo, -1, -1, -1, ifoo, info)
        ASSERT(info == -6)

        call sgepiv('r', 'f', 1, 2, foo, 1, 0, -1, ifoo, info)
        ASSERT(info == -7)

        call sgepiv('r', 'f', 1, 2, foo, 1, 2, -1, ifoo, info)
        ASSERT(info == -7)

        call sgepiv('c', 'f', 1, 2, foo, 1, 3, -1, ifoo, info)
        ASSERT(info == -7)

        call sgepiv('r', 'f', 1, 2, foo, 1, 1, 0, ifoo, info)
        ASSERT(info == -8)

        call sgepiv('r', 'f', 1, 2, foo, 1, 1, 2, ifoo, info)
        ASSERT(info == -8)

        call sgepiv('c', 'f', 1, 2, foo, 1, 1, 3, ifoo, info)
        ASSERT(info == -8)
    end subroutine test_sanity

    subroutine test(seed)
        integer, intent(in) :: seed

        character(1)          :: what, dir
        integer               :: info, m, n, ldA, from, to, i
        type(prng)            :: rng
        integer,  allocatable :: ipiv(:), order(:)
        real(WP), allocatable :: A(:), tmp(:)

        call rng%init(seed, info)
        ASSERT(info == 0)

    !-- what = 'r', dir = 'f' --------------------------------------------------
        what = 'r'
        dir = 'f'
        m = 5
        n = 6
        lda = m + 1
        from = 1
        to = 4
        allocate(A(ldA*n), tmp(ldA*n), ipiv(to-from+1), order(m))
        call rng%snormal(lda*n, A, ZERO, ONE, info)
        call slacpy('a', lda, n, A, lda, tmp, ldA)
        call arange(m, order, 1, 1, info)
        ipiv(1) = 4
        ipiv(2) = 5
        ipiv(3) = 3
        ipiv(4) = 3
        call sgepiv(what, dir, m, n, A, ldA, from, to, ipiv, info, order)
        ASSERT(info == 0)
        ASSERT(order(1) == 4)
        ASSERT(order(2) == 5)
        ASSERT(order(3) == 1)
        ASSERT(order(4) == 3)
        ASSERT(order(5) == 2)
        do i = 1, m
            ASSERT(all_close(n, A(i:), ldA, tmp(order(i):), ldA, info, rtol=1e-5_WP))
        end do
        do i = m+1, ldA
            ASSERT(all_close(n, A(i:), ldA, tmp(i:), ldA, info, rtol=1e-5_WP))
        end do

    !-- what = 'r', dir = 'b' --------------------------------------------------
        dir = 'b'
        call sgepiv(what, dir, m, n, A, ldA, from, to, ipiv, info, order)
        ASSERT(info == 0)
        do i = 1, m
            ASSERT(order(i) == i)
            ASSERT(all_close(n, A(i:), ldA, tmp(order(i):), ldA, info, rtol=1e-5_WP))
        end do
        do i = m+1, ldA
            ASSERT(all_close(n, A(i:), ldA, tmp(i:), ldA, info, rtol=1e-5_WP))
        end do
        deallocate(A, tmp, ipiv, order)

    !-- what = 'c', dir = 'f' --------------------------------------------------
        what = 'c'
        dir = 'f'
        m = 6
        n = 5
        lda = m + 1
        from = 2
        to = 5
        allocate(A(ldA*n), tmp(ldA*n), ipiv(to-from+1), order(n))
        call rng%snormal(lda*n, A, ZERO, ONE, info)
        call slacpy('a', lda, n, A, lda, tmp, ldA)
        call arange(n, order, 1, 1, info)
        ipiv(1) = 5
        ipiv(2) = 1
        ipiv(3) = 3
        ipiv(4) = 4
        call sgepiv(what, dir, m, n, A, ldA, from, to, ipiv, info, order)
        ASSERT(info == 0)
        ASSERT(order(1) == 3)
        ASSERT(order(2) == 5)
        ASSERT(order(3) == 4)
        ASSERT(order(4) == 2)
        ASSERT(order(5) == 1)
        do i = 1, n
            ASSERT(all_close(m, A(1 + (i-1)*ldA:), 1, tmp(1 + (order(i)-1)*ldA:), 1, info, rtol=1e-5_WP))
            ASSERT(all_close(lda-m, A(m+1 + (i-1)*ldA:), 1, tmp(m+1 + (i-1)*ldA:), 1, info, rtol=1e-5_WP))
        end do

    !-- what = 'c', dir = 'b' --------------------------------------------------
        dir = 'b'
        call sgepiv(what, dir, m, n, A, ldA, from, to, ipiv, info, order)
        ASSERT(info == 0)
        do i = 1, n
            ASSERT(order(i) == i)
            ASSERT(all_close(lda, A(1 + (i-1)*ldA:), 1, tmp(1 + (i-1)*ldA:), 1, info, rtol=1e-5_WP))
        end do

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test

    subroutine test_runtime_err()
        character(1) :: what, dir
        integer      :: info, m, n, ldA, from, to, ipiv(2)
        real(WP)     :: foo(1)
        type(prng)   :: rng
        real(WP), allocatable :: A(:)

        call rng%init(seed, info)
        ASSERT(info == 0)

    !-- what = 'r', dir = 'f' --------------------------------------------------
        what = 'r'
        dir = 'f'
        m = 5
        n = 6
        lda = m + 1
        from = 3
        to = 4
        allocate(A(lda*n))
        ipiv(1) = 0
        ipiv(2) = m + 1
        call sgepiv(what, dir, m, n, foo, ldA, from, to, ipiv, info)
        ASSERT(info == 1)
        ipiv(1) = 1
        call sgepiv(what, dir, m, n, foo, ldA, from, to, ipiv, info)
        ASSERT(info == 2)
        deallocate(A)

    !-- what = 'c', dir = 'f' --------------------------------------------------
        what = 'c'
        dir = 'f'
        m = 5
        n = 6
        lda = n + 1
        from = 3
        to = 4
        allocate(A(lda*n))
        ipiv(1) = 0
        ipiv(2) = n + 1
        call sgepiv(what, dir, m, n, foo, ldA, from, to, ipiv, info)
        ASSERT(info == 1)
        ipiv(1) = 1
        call sgepiv(what, dir, m, n, foo, ldA, from, to, ipiv, info)
        ASSERT(info == 2)
        deallocate(A)
    end subroutine test_runtime_err
end program sgepiv_test
