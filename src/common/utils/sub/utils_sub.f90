!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the implementation of the [[maria_utils_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Implements the [[maria_utils_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
submodule (maria_utils_mod) maria_utils_sub
implicit none (type, external)

integer, parameter :: INSERTION_SORT_SIZE = 16
integer, parameter :: QUICKSORT_PARTITION_STACK_SIZE = 32

contains
    !------------------------------------------------------------------------------------------------------------------------

    module function loptional_val &
    (default_val, val)
    !-- Input/output arguments -------------------------------------------------
        logical, intent(in)           :: default_val
        logical, intent(in), optional :: val
        logical                       :: loptional_val

    !-- Executable section -----------------------------------------------------
        if (present(val)) then
            loptional_val = val
        else
            loptional_val = default_val
        end if
    end function loptional_val

    !------------------------------------------------------------------------------------------------------------------------

    module function ioptional_val &
    (default_val, val)
    !-- Input/output arguments -------------------------------------------------
        integer, intent(in)           :: default_val
        integer, intent(in), optional :: val
        integer                       :: ioptional_val

    !-- Executable section -----------------------------------------------------
        if (present(val)) then
            ioptional_val = val
        else
            ioptional_val = default_val
        end if
    end function ioptional_val

    !------------------------------------------------------------------------------------------------------------------------

    module function soptional_val &
    (default_val, val)
    use maria_kinds_mod, only: &
        WP => SP
    !-- Input/output arguments -------------------------------------------------
        real(WP), intent(in)           :: default_val
        real(WP), intent(in), optional :: val
        real(WP)                       :: soptional_val

    !-- Executable section -----------------------------------------------------
        if (present(val)) then
            soptional_val = val
        else
            soptional_val = default_val
        end if
    end function soptional_val

    !------------------------------------------------------------------------------------------------------------------------

    module function doptional_val &
    (default_val, val)
    use maria_kinds_mod, only: &
        WP => DP
    !-- Input/output arguments -------------------------------------------------
        real(WP), intent(in)           :: default_val
        real(WP), intent(in), optional :: val
        real(WP)                       :: doptional_val

    !-- Executable section -----------------------------------------------------
        if (present(val)) then
            doptional_val = val
        else
            doptional_val = default_val
        end if
    end function doptional_val

    !------------------------------------------------------------------------------------------------------------------------

    module function sare_close &
    (a, b, info, atol, rtol)
    use maria_kinds_mod,      only: &
        WP => SP
    use maria_constants_mod,  only: &
        EPS => S_MACHTOL,           &
        ZERO => S_ZERO
    use maria_comparison_mod, only: &
        safe_less
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        real(WP), intent(in)           :: a
        real(WP), intent(in)           :: b
        integer,  intent(out)          :: info
        real(WP), intent(in), optional :: atol
        real(WP), intent(in), optional :: rtol
        logical                        :: sare_close

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SARE_CLOSE'
        real(WP)                :: atol_, rtol_, ref

    !-- Process optional arguments and their default values --------------------
        atol_ = optional_val(EPS, atol)
        rtol_ = optional_val(ZERO, rtol)

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, atol_, ZERO)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, rtol_, ZERO)) then
            info = -5
            exit sanity
        end if
    end block sanity

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            sare_close = .false.
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        ref = max(rtol_ * max(abs(a), abs(b)), atol_)
        sare_close = safe_less(abs(a - b), ref)
    end function sare_close

    !------------------------------------------------------------------------------------------------------------------------

    module function dare_close &
    (a, b, info, atol, rtol)
    use maria_kinds_mod,      only: &
        WP => DP
    use maria_constants_mod,  only: &
        EPS => D_MACHTOL,           &
        ZERO => D_ZERO
    use maria_comparison_mod, only: &
        safe_less
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        real(WP), intent(in)           :: a
        real(WP), intent(in)           :: b
        integer,  intent(out)          :: info
        real(WP), intent(in), optional :: atol
        real(WP), intent(in), optional :: rtol
        logical                        :: dare_close

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DARE_CLOSE'
        real(WP)                :: atol_, rtol_, ref

    !-- Process optional arguments and their default values --------------------
        atol_ = optional_val(EPS, atol)
        rtol_ = optional_val(ZERO, rtol)

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, atol_, ZERO)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, rtol_, ZERO)) then
            info = -5
            exit sanity
        end if
    end block sanity

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            dare_close = .false.
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        ref = max(rtol_ * max(abs(a), abs(b)), atol_)
        dare_close = safe_less(abs(a - b), ref)
    end function dare_close

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine iswap_pair &
    (a, b)
    !-- Input/output arguments -------------------------------------------------
        integer, intent(inout) :: a
        integer, intent(inout) :: b

    !-- Inner variables --------------------------------------------------------
        integer :: tmp

    !-- Executable section -----------------------------------------------------
        tmp = a
        a = b
        b = tmp
    end subroutine iswap_pair

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine sswap_pair &
    (a, b)
    use maria_kinds_mod, only: &
        WP => SP
    !-- Input/output arguments -------------------------------------------------
        real(WP), intent(inout) :: a
        real(WP), intent(inout) :: b

    !-- Inner variables --------------------------------------------------------
        real(WP) :: tmp

    !-- Executable section -----------------------------------------------------
        tmp = a
        a = b
        b = tmp
    end subroutine sswap_pair

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dswap_pair &
    (a, b)
    use maria_kinds_mod, only: &
        WP => DP
    !-- Input/output arguments -------------------------------------------------
        real(WP), intent(inout) :: a
        real(WP), intent(inout) :: b

    !-- Inner variables --------------------------------------------------------
        real(WP) :: tmp

    !-- Executable section -----------------------------------------------------
        tmp = a
        a = b
        b = tmp
    end subroutine dswap_pair

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine slinspace &
    (n, x, a, b, info, include_b, step)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        ZERO => S_ZERO
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS,               &
        BAD_IF_SAME
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)              :: n
        real(WP), intent(out), contiguous :: x(:)
        real(WP), intent(in)              :: a
        real(WP), intent(in)              :: b
        integer,  intent(out)             :: info
        logical,  intent(in),  optional   :: include_b
        real(WP), intent(out), optional   :: step

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SLINSPACE'
        logical                 :: include_b_
        integer                 :: i
        real(WP)                :: h

    !-- Process optional arguments and their default values --------------------
        include_b_ = optional_val(.true., include_b)

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -1
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, b, a)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_SAME, b, a)) then
            info = -4
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (present(step)) then
            step = ZERO
        endif
        if (n == 0) return
        if (n == 1) then
            x(1) = a
            return
        end if
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (include_b_) then
            h = (b - a) / (n - 1)
        else
            h = (b - a) / n
        end if
        
        do i = 1, n
            x(i) = a + h * (i - 1)
        end do

        if (present(step)) then
            step = h
        end if
    end subroutine slinspace

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dlinspace &
    (n, x, a, b, info, include_b, step)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        ZERO => D_ZERO
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS,               &
        BAD_IF_SAME
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)              :: n
        real(WP), intent(out), contiguous :: x(:)
        real(WP), intent(in)              :: a
        real(WP), intent(in)              :: b
        integer,  intent(out)             :: info
        logical,  intent(in),  optional   :: include_b
        real(WP), intent(out), optional   :: step

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DLINSPACE'
        logical                 :: include_b_
        integer                 :: i
        real(WP)                :: h

    !-- Process optional arguments and their default values --------------------
        include_b_ = optional_val(.true., include_b)

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -1
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, b, a)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_SAME, b, a)) then
            info = -4
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (present(step)) then
            step = ZERO
        endif
        if (n == 0) return
        if (n == 1) then
            x(1) = a
            return
        end if
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------   
        if (include_b_) then
            h = (b - a) / (n - 1)
        else
            h = (b - a) / n
        end if
        
        do i = 1, n
            x(i) = a + h * (i - 1)
        end do

        if (present(step)) then
            step = h
        end if
    end subroutine dlinspace

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine iarange &
    (n, x, a, step, info)
    use maria_argcheck_mod, only: &
        arg_is_bad,               &
        BAD_IF_LESS
    use maria_reports_mod,  only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer, intent(in)              :: n
        integer, intent(out), contiguous :: x(:)
        integer, intent(in)              :: a
        integer, intent(in)              :: step
        integer, intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'IARANGE'
        integer                 :: i

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -1
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (n == 0) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------         
        do i = 1, n
            x(i) = a + step * (i - 1)
        end do
    end subroutine iarange

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine sarange &
    (n, x, a, step, info)
    use maria_kinds_mod,    only: &
        WP => SP
    use maria_argcheck_mod, only: &
        arg_is_bad,               &
        BAD_IF_LESS
    use maria_reports_mod,  only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)              :: n
        real(WP), intent(out), contiguous :: x(:)
        real(WP), intent(in)              :: a
        real(WP), intent(in)              :: step
        integer,  intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SARANGE'
        integer                 :: i

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -1
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (n == 0) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------               
        do i = 1, n
            x(i) = a + step * (i - 1)
        end do
    end subroutine sarange

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine darange &
    (n, x, a, step, info)
    use maria_kinds_mod,    only: &
        WP => DP
    use maria_argcheck_mod, only: &
        arg_is_bad,               &
        BAD_IF_LESS
    use maria_reports_mod,  only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)              :: n
        real(WP), intent(out), contiguous :: x(:)
        real(WP), intent(in)              :: a
        real(WP), intent(in)              :: step
        integer,  intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DARANGE'
        integer                 :: i

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -1
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (n == 0) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------               
        do i = 1, n
            x(i) = a + step * (i - 1)
        end do
    end subroutine darange

    !------------------------------------------------------------------------------------------------------------------------

    module function imedian_of_3 &
    (a, b, c)
    !-- Input/output arguments -------------------------------------------------
        integer, intent(in) :: a
        integer, intent(in) :: b
        integer, intent(in) :: c 
        integer             :: imedian_of_3

    !-- Executable section -----------------------------------------------------               
        if (a < b) then
            if (c < a) then
                imedian_of_3 = a
            else if (c < b) then
                imedian_of_3 = c
            else
                imedian_of_3 = b
            end if
        else
            if (c < b) then
                imedian_of_3 = b
            else if (c < a) then
                imedian_of_3 = c
            else
                imedian_of_3 = a
            end if
        end if
    end function imedian_of_3

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine isort &
    (id, n, x, info)
    use maria_argcheck_mod, only: &
        arg_is_bad,               &
        BAD_IF_LESS
    use maria_reports_mod,  only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        character(1), intent(in)                :: id
        integer,      intent(in)                :: n
        integer,      intent(inout), contiguous :: x(:)
        integer,      intent(out)               :: info
    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'ISORT'
        logical                 :: increasing
        integer                 :: i, j, dir, left, right, mid, top, pivot, tmp, &
                                   partition_size 
        integer                 :: partition_stack(2, QUICKSORT_PARTITION_STACK_SIZE)

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        select case (id)
            case ('i', 'I')
                increasing = .true.
            case ('d', 'D')
                increasing = .false.
            case default
                increasing = .true.
                info = -1
                exit sanity
        end select
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -2
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (n == 0 .or. n == 1) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if
 
    !-- Executable section -----------------------------------------------------                  
        if (increasing) then
            dir = 1
        else
            dir = -1
        end if

        top = 1
        partition_stack(1, top) = 1
        partition_stack(2, top) = n
        partition_loop: do while (top > 0)
        !-- Pick the next partition from stack ---------------------------------
            left = partition_stack(1, top)
            right = partition_stack(2, top)
            top = top - 1

            partition_size = right - left + 1
            if (partition_size <= INSERTION_SORT_SIZE) then
            !-- Apply insertion sort to d(left:right) --------------------------
                forward: do i = left+1, right
                    backward: do j = i, left+1, -1
                        tmp = dir * (x(j) - x(j-1))
                        if (tmp < 0) then
                            call swap_pair(x(j), x(j-1))
                        else
                            exit backward
                        end if
                    end do backward
                end do forward
            else
            !-- Partition d(left:right) ----------------------------------------
                mid = left + partition_size / 2
                pivot = median_of_3(x(left), x(mid), x(right))
                
                i = left - 1
                j = right + 1
                order_wrt_pivot: do
                    do
                        j = j - 1
                        tmp = dir * (x(j) - pivot)
                        if (tmp <= 0) exit
                    end do

                    do
                        i = i + 1
                        tmp = dir * (x(i) - pivot) 
                        if (tmp >= 0) exit
                    end do

                    if (i < j) then
                        call swap_pair(x(i), x(j))
                    else
                        exit order_wrt_pivot
                    end if
                end do order_wrt_pivot

            !-- Stack the smaller partition on top -----------------------------
                if (j - left > right - j - 1) then
                    top = top + 1
                    partition_stack(1, top) = left
                    partition_stack(2, top) = j
                    top = top + 1
                    partition_stack(1, top) = j + 1
                    partition_stack(2, top) = right
                else
                    top = top + 1
                    partition_stack(1, top) = j + 1
                    partition_stack(2, top) = right
                    top = top + 1
                    partition_stack(1, top) = left
                    partition_stack(2, top) = j
                end if 
            end if 
        end do partition_loop
    end subroutine isort

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine iselect_no_replacement &
    (rng, n, x, k, info)
    use maria_utils_mod,    only: &
        swap_pair
    use maria_argcheck_mod, only: &
        arg_is_bad,               &
        BAD_IF_LESS,              &
        BAD_IF_MORE
    use maria_reports_mod,  only: &
        report_bad_arg
    use maria_prng_mod,     only: &
        prng
    !-- Input/output arguments -------------------------------------------------
        class(prng), intent(in)                :: rng        
        integer,     intent(in)                :: n
        integer,     intent(inout), contiguous :: x(:)
        integer,     intent(in)                :: k
        integer,     intent(out)               :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'ISELECT_NO_REPLACEMENT'
        integer                 :: i, next(1)

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, k, 0)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, k, max(n-1,0))) then
            info = -4
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (n == 0 .or. n == 1 .or. k == 0) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        do i = 1, k
            call rng%iuniform(1, next, i, n, info)
            if (i /= next(1)) then
                call swap_pair(x(i), x(next(1)))
            end if
        end do
    end subroutine iselect_no_replacement

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine ipermute &
    (rng, n, x, info)
    use maria_argcheck_mod, only: &
        arg_is_bad,               &
        BAD_IF_LESS
    use maria_reports_mod,  only: &
        report_bad_arg
    use maria_prng_mod,     only: &
        prng
    !-- Input/output arguments -------------------------------------------------
        class(prng), intent(in)                :: rng       
        integer,     intent(in)                :: n
        integer,     intent(inout), contiguous :: x(:)
        integer,     intent(out)               :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'IPERMUTE'

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -2
            exit sanity
        end if
    end block sanity
 
    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (n == 0 .or. n == 1) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        call select_no_replacement(rng, n, x, n-1, info)
    end subroutine ipermute

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine srealloc &
    (copy, x, n, ierr, factor)
    use maria_kinds_mod,    only:   &
        WP => SP
    use maria_argcheck_mod, only: &
        arg_is_bad,               &
        BAD_IF_LESS
    use maria_reports_mod,  only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        logical,    intent(in   )               :: copy
        real(WP),   intent(inout),  allocatable :: x(:)
        integer,    intent(in   )               :: n
        integer,    intent(  out)               :: ierr
        integer,    intent(in   ),  optional    :: factor

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SREALLOC'
        integer :: new_size
        real(WP),   allocatable :: tmp(:)

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        ierr = 0
        if (arg_is_bad(BAD_IF_LESS, n, 1)) then
            ierr = -3
            exit sanity
        end if
        if (present(factor)) then
            if (arg_is_bad(BAD_IF_LESS, factor, 2)) then
                ierr = -5
                exit sanity
            end if
        end if
    end block sanity
 
    !-- Report bad input -------------------------------------------------------
        if (ierr /= 0) then
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------
    if ( .not. allocated(x) ) then
        allocate(x(n), stat=ierr)
        return
    end if
    
    if ( size(x) >= n ) return

    if ( present(factor) ) then
        new_size = size(x)
        do
            new_size = factor * new_size
            if ( new_size >= n ) exit
        end do
    else
        new_size = n
    end if

    if ( copy ) then
        allocate( tmp(new_size), stat=ierr )
        tmp(1:size(x)) = x
        call move_alloc(tmp, x)
    else
        deallocate( x )
        allocate( x(new_size), stat=ierr )
    end if
    end subroutine srealloc

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine drealloc &
    (copy, x, n, ierr, factor)
    use maria_kinds_mod,    only:   &
        WP => DP
    use maria_argcheck_mod, only: &
        arg_is_bad,               &
        BAD_IF_LESS
    use maria_reports_mod,  only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        logical,    intent(in   )               :: copy
        real(WP),   intent(inout),  allocatable :: x(:)
        integer,    intent(in   )               :: n
        integer,    intent(  out)               :: ierr
        integer,    intent(in   ),  optional    :: factor

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DREALLOC'
        integer :: new_size
        real(WP),   allocatable :: tmp(:)

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        ierr = 0
        if (arg_is_bad(BAD_IF_LESS, n, 1)) then
            ierr = -3
            exit sanity
        end if
        if (present(factor)) then
            if (arg_is_bad(BAD_IF_LESS, factor, 2)) then
                ierr = -5
                exit sanity
            end if
        end if
    end block sanity
 
    !-- Report bad input -------------------------------------------------------
        if (ierr /= 0) then
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------
    if ( .not. allocated(x) ) then
        allocate(x(n), stat=ierr)
        return
    end if
    
    if ( size(x) >= n ) return

    if ( present(factor) ) then
        new_size = size(x)
        do
            new_size = factor * new_size
            if ( new_size >= n ) exit
        end do
    else
        new_size = n
    end if

    if ( copy ) then
        allocate( tmp(new_size), stat=ierr )
        tmp(1:size(x)) = x
        call move_alloc(tmp, x)
    else
        deallocate( x )
        allocate( x(new_size), stat=ierr )
    end if
    end subroutine drealloc

    !------------------------------------------------------------------------------------------------------------------------
end submodule maria_utils_sub
