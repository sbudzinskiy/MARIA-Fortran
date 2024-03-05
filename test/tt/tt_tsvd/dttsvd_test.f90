program dttsvd_test
#include "maria_assert_mod.h"
#ifdef MARIA_MKL
use maria_prng_mkl_mod,     only: &
    prng => prng_mkl
#else
use maria_prng_builtin_mod, only: &
    prng => prng_builtin
#endif
use maria_kinds_mod,        only: &
    WP => DP
use maria_arr_mod,          only: &
    AR => darr
use maria_constants_mod, only: &
    ONE => D_ONE
use maria_comparison_mod, only: &
    safe_eq, safe_leq
use maria_la_utils_mod, only: &
    dall_const
use maria_la_core_mod, only: &
    daxpy, &
    dcopy, &
    dnrm2
use maria_access_tensor_mod,    only: &
    mi2i
use maria_tt_utils_mod,        only: &
    dtt2full
use maria_tt_tsvd_mod,  only: &
    dttsvd
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
        integer :: ierr, n(3), maxr(2), ifoo(1)
        real(WP) :: foo(1)
        type(AR) :: cores(1)

        n = 1
        maxr = 1

        call dttsvd('x', 1, n, foo, ifoo, cores, foo, 0, ifoo, 0, ierr, maxr=maxr, rtolf=-ONE)
        ASSERT( ierr == -1 )

        call dttsvd('l', 1, n, foo, ifoo, cores, foo, 0, ifoo, 0, ierr, maxr=maxr, rtolf=-ONE)
        ASSERT( ierr == -2 )

        n(1) = -1
        call dttsvd('l', 3, n, foo, ifoo, cores, foo, 0, ifoo, 0, ierr, maxr=maxr, rtolf=-ONE)
        ASSERT( ierr == -3 )

        n(1) = 1
        maxr(1) = -1
        call dttsvd('l', 3, n, foo, ifoo, cores, foo, 0, ifoo, 0, ierr, maxr=maxr, rtolf=-ONE)
        ASSERT( ierr == -12 )

        maxr(1) = 1
        call dttsvd('l', 3, n, foo, ifoo, cores, foo, 0, ifoo, 0, ierr, maxr=maxr, rtolf=-ONE)
        ASSERT( ierr == -13 )
    end subroutine test_sanity

    subroutine test_quick()
        integer :: ierr, n(3), maxr(2), r(2), ifoo(1)
        real(WP) :: foo(1)
        type(AR) :: cores(1)

        n = 2
        maxr = 3

        n(1) = 0
        call dttsvd('l', 3, n, foo, ifoo, cores, foo, 0, ifoo, 0, ierr, maxr=maxr)
        ASSERT( ierr == 0 )

        n(1) = 2
        maxr(1) = 0
        call dttsvd('l', 3, n, foo, r, cores, foo, 0, ifoo, 0, ierr, maxr=maxr)
        ASSERT( ierr == 0 )
        ASSERT( r(1) == 0 .and. r(2) == 0 )
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer              :: d, ierr, lwork, liwork, numel, i, j, k, l, ind, ifoo(1), mi3(4)
        real(WP)    :: foo(1), rtolf, rerrf, aerrf, nrmA
        type(prng) :: rng
        integer, allocatable :: n(:), r(:), maxr(:), iwork(:)
        real(WP), allocatable :: A(:), work(:), tmp(:)
        type(AR), allocatable :: cores(:)

        call rng%init(seed, ierr)
        ASSERT( ierr == 0 )

    !-- i1 + i2 + ... + id
        d = 4
        allocate( n(d), r(d-1), maxr(d-1), cores(d) )
        n(1) = 10
        n(2) = 5
        n(3) = 8
        n(4) = 10
        numel = product(n)
        rtolf = 1e-12_WP

        call dttsvd('l', d, n, foo, r, cores, foo, -1, ifoo, -1, ierr, rtolf=rtolf)
        lwork = int(foo(1))
        liwork = ifoo(1)

        allocate( A(numel), tmp(numel), work(lwork), iwork(liwork) )        
        do i = 1, n(1)
            do j = 1, n(2)
                do k = 1, n(3)
                    do l = 1, n(4)
                        mi3(1) = i
                        mi3(2) = j
                        mi3(3) = k
                        mi3(4) = l
                        call mi2i(4, n, mi3, ind, ierr)
                        A(ind) = ONE * (i + j + k + l)
                    end do
                end do
            end do
        end do
        call dcopy(numel, A, 1, tmp, 1)

        call dttsvd('l', d, n, tmp, r, cores, work, lwork, iwork, liwork, ierr, rtolf=rtolf, rerrf=rerrf)
        ASSERT( ierr == 0 )
        ASSERT( all(r == 2) )
        ASSERT( safe_leq(rerrf, rtolf) )
        call dtt2full(d, n, r, cores, tmp, work, lwork, ierr)
        call daxpy(numel, -ONE, A, 1, tmp, 1)
        nrmA = dnrm2(numel, A, 1)
        aerrf = dnrm2(numel, tmp, 1)
        ASSERT( safe_leq(rerrf, aerrf / nrmA) )

        maxr(1) = 1
        maxr(2) = 2
        maxr(3) = 2

        call dcopy(numel, A, 1, tmp, 1)
        call dttsvd('r', d, n, tmp, r, cores, work, lwork, iwork, liwork, ierr, maxr=maxr, rerrf=rerrf)
        ASSERT( ierr == 0 )
        ASSERT( all(r <= maxr) )
        call dtt2full(d, n, r, cores, tmp, work, lwork, ierr)
        call daxpy(numel, -ONE, A, 1, tmp, 1)
        aerrf = dnrm2(numel, tmp, 1)
        ASSERT( safe_leq(rerrf, aerrf / nrmA) )

        do k = 1, d
            deallocate( cores(k)%arr )
        end do
        deallocate( A, tmp, work, iwork, n, r, cores )
    end subroutine test
end program dttsvd_test
