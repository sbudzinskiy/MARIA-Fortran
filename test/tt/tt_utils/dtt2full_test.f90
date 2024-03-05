program dtt2full_test
#include "maria_assert_mod.h"
use maria_kinds_mod,        only: &
    WP => DP
use maria_arr_mod,          only: &
    AR => darr
use maria_constants_mod, only: &
    ZERO => D_ZERO, &
    ONE => D_ONE
use maria_comparison_mod, only: &
    safe_eq
use maria_la_utils_mod, only: &
    dall_const
use maria_access_tensor_mod,    only: &
    mi2i
use maria_tt_utils_mod,        only: &
    dtt2full
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
        integer :: ierr, n(3), r(2)
        real(WP) :: foo(1)
        type(AR) :: cores(1)

        n = 1
        r = 1

        call dtt2full(1, n, r, cores, foo, foo, 0, ierr)
        ASSERT( ierr == -1 )

        n(1) = -1
        call dtt2full(3, n, r, cores, foo, foo, 0, ierr)
        ASSERT( ierr == -2 )

        n(1) = 1
        r(1) = -1
        call dtt2full(3, n, r, cores, foo, foo, 0, ierr)
        ASSERT( ierr == -3 )
    end subroutine test_sanity

    subroutine test_quick()
        integer :: ierr, n(3), r(2), numel
        real(WP) :: foo(1), A(2*2*2)
        type(AR) :: cores(1)

        n = 2
        r = 3
        numel = 2*2*2

        n(1) = 0
        call dtt2full(3, n, r, cores, A, foo, 0, ierr)
        ASSERT( ierr == 0 )

        r(2) = 0
        n(1) = 2
        call dtt2full(3, n, r, cores, A, foo, 0, ierr)
        ASSERT( ierr == 0 )
        ASSERT( dall_const(numel, A, 1, ZERO, ierr) )
    end subroutine test_quick

    subroutine test()
        integer              :: d, ierr, lwork, numel, k, i, j, l, ind, mi3(4), n3(4)
        real(WP)    :: foo(1)
        integer, allocatable :: n(:), r(:)
        real(WP), allocatable :: A(:), work(:)
        type(AR), allocatable :: cores(:)

        d = 4
        allocate( n(d), r(d-1), cores(d) )
        n(1) = 2
        n(2) = 3
        n(3) = 4
        n(4) = 5

    !-- i1 + i2 + ... + id
        r = 2
        numel = product(n)
        call dtt2full(d, n, r, cores, foo, foo, -1, ierr)
        lwork = int(foo(1))

        allocate( A(numel), work(lwork) )
        allocate( cores(1)%arr(n(1)*r(1)) )
        do k = 2, d-1
            allocate( cores(k)%arr(r(k-1)*n(k)*r(k)) )
        end do
        allocate( cores(d)%arr(r(d-1)*n(d)) )

        ! First core
        n3(1) = 1
        n3(2) = n(1)
        n3(3) = r(1)
        do j = 1, n(1)
            do l = 1, r(1)
                mi3(1) = 1
                mi3(2) = j
                mi3(3) = l
                call mi2i(3, n3, mi3, ind, ierr)
                if ( l == 1 ) then
                    cores(1)%arr(ind) = j
                else
                    cores(1)%arr(ind) = 1
                end if
            end do
        end do
        
        ! Middle cores
        do k = 2, d-1
            n3(1) = r(k-1)
            n3(2) = n(k)
            n3(3) = r(k)
            do i = 1, r(k-1)
                do j = 1, n(k)
                    do l = 1, r(k)
                        mi3(1) = i
                        mi3(2) = j
                        mi3(3) = l
                        call mi2i(3, n3, mi3, ind, ierr) 
                        if ( i == 1 .and. l == 1 ) then
                            cores(k)%arr(ind) = 1
                        elseif ( i == 1 .and. l == 2 ) then
                            cores(k)%arr(ind) = 0
                        elseif ( i == 2 .and. l == 1 ) then
                            cores(k)%arr(ind) = j
                        else
                            cores(k)%arr(ind) = 1
                        end if
                    end do
                end do
            end do
        end do
        
        ! Last core
        n3(1) = r(d-1)
        n3(2) = n(d)
        n3(3) = 1
        do i = 1, r(d-1)
            do j = 1, n(d)
                mi3(1) = i
                mi3(2) = j
                mi3(3) = 1
                call mi2i(3, n3, mi3, ind, ierr) 
                if ( i == 1 ) then
                    cores(d)%arr(ind) = 1
                else
                    cores(d)%arr(ind) = j
                end if
            end do
        end do
        
        call dtt2full(d, n, r, cores, A, work, lwork, ierr)
        ASSERT( ierr == 0 )

        do i = 1, n(1)
            do j = 1, n(2)
                do k = 1, n(3)
                    do l = 1, n(4)
                        mi3(1) = i
                        mi3(2) = j
                        mi3(3) = k
                        mi3(4) = l
                        call mi2i(4, n, mi3, ind, ierr)
                        ASSERT( safe_eq(ONE*(i + j + k + l), A(ind)) )
                    end do
                end do
            end do
        end do
    end subroutine test
end program dtt2full_test
