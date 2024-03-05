!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the implementation of the [[maria_access_tensor_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Implements the [[maria_access_tensor_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------

submodule (maria_access_tensor_mod) maria_access_tensor_sub
implicit none (type, external)

contains
    !------------------------------------------------------------------------------------------------------------------------

    module subroutine mi2i  &
    (d, n, mi, i, ierr)
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS,               &
        BAD_IF_MORE
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer,    intent(in   )               :: d
        integer,    intent(in   ),  contiguous  :: n(:)
        integer,    intent(in   ),  contiguous  :: mi(:)
        integer,    intent(  out)               :: i
        integer,    intent(  out)               :: ierr

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "MI2I"
        integer                 :: k

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        ierr = 0
        if (arg_is_bad(BAD_IF_LESS, d, 2)) then
            ierr = -1
            exit sanity
        end if
        do k = 1, d
            if (arg_is_bad(BAD_IF_LESS, n(k), 1)) then
                ierr = -2
                exit sanity
            end if
        end do
        do k = 1, d
            if (arg_is_bad(BAD_IF_LESS, mi(k), 1)) then
                ierr = -3
                exit sanity
            end if
            if (arg_is_bad(BAD_IF_MORE, mi(k), n(k))) then
                ierr = -3
                exit sanity
            end if
        end do
    end block sanity

    !-- Report bad input -------------------------------------------------------
        if (ierr /= 0) then
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------
        i = (mi(d) - 1) * n(d-1)
        do k = d-1, 2, -1
            i = (i + mi(k) - 1) * n(k-1)
        end do
        i = i + mi(1)
    end subroutine mi2i

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine i2mi  &
    (d, n, i, mi, ierr)
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS,               &
        BAD_IF_MORE
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer,    intent(in   )               :: d
        integer,    intent(in   ),  contiguous  :: n(:)
        integer,    intent(in   )               :: i
        integer,    intent(  out),  contiguous  :: mi(:)
        integer,    intent(  out)               :: ierr

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "I2MI"
        integer                 :: k, tmp

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        ierr = 0
        if (arg_is_bad(BAD_IF_LESS, d, 2)) then
            ierr = -1
            exit sanity
        end if
        do k = 1, d
            if (arg_is_bad(BAD_IF_LESS, n(k), 1)) then
                ierr = -2
                exit sanity
            end if
        end do
        if (arg_is_bad(BAD_IF_LESS, i, 1)) then
            ierr = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, i, product(n(1:d)))) then
            ierr = -3
            exit sanity
        end if
    end block sanity

    !-- Report bad input -------------------------------------------------------
        if (ierr /= 0) then
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------
    tmp = i - 1
    do k = 1, d
        mi(k) = mod(tmp, n(k)) + 1
        tmp = tmp / n(k)
    end do
    end subroutine i2mi

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine stenval2fib &
    (tenval, d, n, k, li, ri, x, incx, iwork, liwork, ierr)
    !-- Kinds, types and constants ---------------------------------------------
    use maria_kinds_mod,        only:   &
        WP => SP
    !-- Auxiliary subroutines --------------------------------------------------
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Arguments --------------------------------------------------------------
        procedure(stenval), intent(in   ),  pointer     :: tenval           !  1
        integer,            intent(in   )               :: d                !  2
        integer,            intent(in   ),  contiguous  :: n(:)             !  3
        integer,            intent(in   )               :: k                !  4
        integer,            intent(in   ),  contiguous  :: li(:)            !  5
        integer,            intent(in   ),  contiguous  :: ri(:)            !  6
        real(WP),           intent(  out),  contiguous  :: x(:)             !  7
        integer,            intent(in   )               :: incx             !  8
        integer,            intent(  out),  contiguous  :: iwork(:)         !  9
        integer,            intent(in   )               :: liwork           ! 10
        integer,            intent(  out)               :: ierr             ! 11
    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'STENVAL2FIB'
    !-- Variables --------------------------------------------------------------
        logical ::  lw_query, bad_arg(11)
        integer ::  i, j, posx, req_liw
    
    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(2) = ( d < 2 )
        do j = 1, d
            bad_arg(3) = bad_arg(3) .or. ( n(j) < 0 )
        end do
        bad_arg(4) = ( k < 1 ) .or. ( k > d )
        
    !-- Quick return if possible -----------------------------------------------
    quick: block
        if ( any(bad_arg) ) exit quick

        if ( any(n(1:d) == 0) ) return
    end block quick

    !-- Estimate workspace -----------------------------------------------------
    workspace: block
        if ( any(bad_arg) ) exit workspace

        req_liw = d

        lw_query = ( liwork == -1 )
        if (lw_query) then
            iwork(1) = req_liw
            return
        end if

        bad_arg(10) = ( liwork < req_liw )
    end block workspace


    !-- Report bad input -------------------------------------------------------
        if ( any(bad_arg) ) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------    
        if ( incx >= 0 ) then
            posx = 1
        else
            posx = 1 - (n(k)-1) * incx
        end if

        associate( mi => iwork )
            do j = 1, n(k)
                do i = 1, k-1
                    mi(i) = li(i)
                end do
                mi(k) = j
                do i = 1, d-k
                    mi(k+i) = ri(i)
                end do
                x(posx) = tenval(d, n, mi, ierr)
                posx = posx + incx
            end do
        end associate
    end subroutine stenval2fib

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dtenval2fib &
    (tenval, d, n, k, li, ri, x, incx, iwork, liwork, ierr)
    !-- Kinds, types and constants ---------------------------------------------
    use maria_kinds_mod,        only:   &
        WP => DP
    !-- Auxiliary subroutines --------------------------------------------------
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Arguments --------------------------------------------------------------
        procedure(dtenval), intent(in   ),  pointer     :: tenval           !  1
        integer,            intent(in   )               :: d                !  2
        integer,            intent(in   ),  contiguous  :: n(:)             !  3
        integer,            intent(in   )               :: k                !  4
        integer,            intent(in   ),  contiguous  :: li(:)            !  5
        integer,            intent(in   ),  contiguous  :: ri(:)            !  6
        real(WP),           intent(  out),  contiguous  :: x(:)             !  7
        integer,            intent(in   )               :: incx             !  8
        integer,            intent(  out),  contiguous  :: iwork(:)         !  9
        integer,            intent(in   )               :: liwork           ! 10
        integer,            intent(  out)               :: ierr             ! 11
    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'DTENVAL2FIB'
    !-- Variables --------------------------------------------------------------
        logical ::  lw_query, bad_arg(11)
        integer ::  i, j, posx, req_liw
    
    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(2) = ( d < 2 )
        do j = 1, d
            bad_arg(3) = bad_arg(3) .or. ( n(j) < 0 )
        end do
        bad_arg(4) = ( k < 1 ) .or. ( k > d )
        
    !-- Quick return if possible -----------------------------------------------
    quick: block
        if ( any(bad_arg) ) exit quick

        if ( any(n(1:d) == 0) ) return
    end block quick

    !-- Estimate workspace -----------------------------------------------------
    workspace: block
        if ( any(bad_arg) ) exit workspace

        req_liw = d

        lw_query = ( liwork == -1 )
        if (lw_query) then
            iwork(1) = req_liw
            return
        end if

        bad_arg(10) = ( liwork < req_liw )
    end block workspace


    !-- Report bad input -------------------------------------------------------
        if ( any(bad_arg) ) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------    
        if ( incx >= 0 ) then
            posx = 1
        else
            posx = 1 - (n(k)-1) * incx
        end if

        associate( mi => iwork )
            do j = 1, n(k)
                do i = 1, k-1
                    mi(i) = li(i)
                end do
                mi(k) = j
                do i = 1, d-k
                    mi(k+i) = ri(i)
                end do
                x(posx) = tenval(d, n, mi, ierr)
                posx = posx + incx
            end do
        end associate
    end subroutine dtenval2fib

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine stenval2full &
    (tenval, d, n, A, iwork, liwork, ierr)
    !-- Kinds, types and constants ---------------------------------------------
    use maria_kinds_mod,        only:   &
        WP => SP
    !-- Auxiliary subroutines --------------------------------------------------
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Arguments --------------------------------------------------------------
        procedure(stenval), intent(in   ),  pointer     :: tenval           !  1
        integer,            intent(in   )               :: d                !  2
        integer,            intent(in   ),  contiguous  :: n(:)             !  3
        real(WP),           intent(  out),  contiguous  :: A(:)             !  4
        integer,            intent(  out),  contiguous  :: iwork(:)         !  5
        integer,            intent(in   )               :: liwork           !  6
        integer,            intent(  out)               :: ierr             !  7
    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'STENVAL2FULL'
    !-- Variables --------------------------------------------------------------
        logical ::  lw_query, bad_arg(7)
        integer ::  j, numel, req_liw
    
    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(2) = ( d < 2 )
        do j = 1, d
            bad_arg(3) = bad_arg(3) .or. ( n(j) < 0 )
        end do
        
    !-- Quick return if possible -----------------------------------------------
    quick: block
        if ( any(bad_arg) ) exit quick

        if ( any(n(1:d) == 0) ) return
    end block quick

    !-- Estimate workspace -----------------------------------------------------
    workspace: block
        if ( any(bad_arg) ) exit workspace

        req_liw = d

        lw_query = ( liwork == -1 )
        if (lw_query) then
            iwork(1) = req_liw
            return
        end if

        bad_arg(6) = ( liwork < req_liw )
    end block workspace


    !-- Report bad input -------------------------------------------------------
        if ( any(bad_arg) ) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------    
        numel = product(n(1:d))
        associate( mi => iwork )
            do j = 1, numel
                call i2mi(d, n, j, mi, ierr)
                A(j) = tenval(d, n, mi, ierr)        
            end do
        end associate
    end subroutine stenval2full

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dtenval2full &
    (tenval, d, n, A, iwork, liwork, ierr)
    !-- Kinds, types and constants ---------------------------------------------
    use maria_kinds_mod,        only:   &
        WP => DP
    !-- Auxiliary subroutines --------------------------------------------------
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Arguments --------------------------------------------------------------
        procedure(dtenval), intent(in   ),  pointer     :: tenval           !  1
        integer,            intent(in   )               :: d                !  2
        integer,            intent(in   ),  contiguous  :: n(:)             !  3
        real(WP),           intent(  out),  contiguous  :: A(:)             !  4
        integer,            intent(  out),  contiguous  :: iwork(:)         !  5
        integer,            intent(in   )               :: liwork           !  6
        integer,            intent(  out)               :: ierr             !  7
    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'DTENVAL2FULL'
    !-- Variables --------------------------------------------------------------
        logical ::  lw_query, bad_arg(7)
        integer ::  j, numel, req_liw
    
    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(2) = ( d < 2 )
        do j = 1, d
            bad_arg(3) = bad_arg(3) .or. ( n(j) < 0 )
        end do
        
    !-- Quick return if possible -----------------------------------------------
    quick: block
        if ( any(bad_arg) ) exit quick

        if ( any(n(1:d) == 0) ) return
    end block quick

    !-- Estimate workspace -----------------------------------------------------
    workspace: block
        if ( any(bad_arg) ) exit workspace

        req_liw = d

        lw_query = ( liwork == -1 )
        if (lw_query) then
            iwork(1) = req_liw
            return
        end if

        bad_arg(6) = ( liwork < req_liw )
    end block workspace


    !-- Report bad input -------------------------------------------------------
        if ( any(bad_arg) ) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------    
        numel = product(n(1:d))
        associate( mi => iwork )
            do j = 1, numel
                call i2mi(d, n, j, mi, ierr)
                A(j) = tenval(d, n, mi, ierr)        
            end do
        end associate
    end subroutine dtenval2full

    !------------------------------------------------------------------------------------------------------------------------

    module function stenval_hilbert &
    (d, n, mi, ierr)
    !-- Kinds, types and constants ---------------------------------------------
    use maria_kinds_mod,        only:   &
        WP => SP
    use maria_constants_mod,    only:   &
        ZERO => S_ZERO,                 &
        ONE => S_ONE
    !-- Auxiliary subroutines --------------------------------------------------
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Arguments --------------------------------------------------------------
        integer,    intent(in   )               :: d                        !  1
        integer,    intent(in   ),  contiguous  :: n(:)                     !  2
        integer,    intent(in   ),  contiguous  :: mi(:)                    !  3
        integer,    intent(  out)               :: ierr                     !  4
        real(WP)                                :: stenval_hilbert
    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'STENVAL_HILBERT'
    !-- Variables --------------------------------------------------------------
        logical ::  bad_arg(4)
        integer ::  j

    !-- Default value ----------------------------------------------------------
        stenval_hilbert = ZERO

    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(1) = ( d < 1 )
        do j = 1, d
            bad_arg(2) = bad_arg(2) .or. ( n(j) < 0 )
        end do
        do j = 1, d
            bad_arg(3) = bad_arg(3) .or. ( mi(j) < 0 ) .or. ( mi(j) > n(j) )
        end do
        
    !-- Quick return if possible -----------------------------------------------
    quick: block
        if ( any(bad_arg) ) exit quick

        if ( any(n(1:d) == 0) ) return
    end block quick

    !-- Report bad input -------------------------------------------------------
        if ( any(bad_arg) ) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------    
        stenval_hilbert = ONE / (sum(mi(1:d)) + d - 1)
    end function stenval_hilbert

    !------------------------------------------------------------------------------------------------------------------------

    module function dtenval_hilbert &
    (d, n, mi, ierr)
    !-- Kinds, types and constants ---------------------------------------------
    use maria_kinds_mod,        only:   &
        WP => DP
    use maria_constants_mod,    only:   &
        ZERO => D_ZERO,                 &
        ONE => D_ONE
    !-- Auxiliary subroutines --------------------------------------------------
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Arguments --------------------------------------------------------------
        integer,    intent(in   )               :: d                        !  1
        integer,    intent(in   ),  contiguous  :: n(:)                     !  2
        integer,    intent(in   ),  contiguous  :: mi(:)                    !  3
        integer,    intent(  out)               :: ierr                     !  4
        real(WP)                                :: dtenval_hilbert
    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'DTENVAL_HILBERT'
    !-- Variables --------------------------------------------------------------
        logical ::  bad_arg(4)
        integer ::  j

    !-- Default value ----------------------------------------------------------
        dtenval_hilbert = ZERO

    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(1) = ( d < 1 )
        do j = 1, d
            bad_arg(2) = bad_arg(2) .or. ( n(j) < 0 )
        end do
        do j = 1, d
            bad_arg(3) = bad_arg(3) .or. ( mi(j) < 0 ) .or. ( mi(j) > n(j) )
        end do
        
    !-- Quick return if possible -----------------------------------------------
    quick: block
        if ( any(bad_arg) ) exit quick

        if ( any(n(1:d) == 0) ) return
    end block quick

    !-- Report bad input -------------------------------------------------------
        if ( any(bad_arg) ) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------    
        dtenval_hilbert = ONE / (sum(mi(1:d)) + d - 1)
    end function dtenval_hilbert

    !------------------------------------------------------------------------------------------------------------------------

    module function stenval_sum &
    (d, n, mi, ierr)
    !-- Kinds, types and constants ---------------------------------------------
    use maria_kinds_mod,        only:   &
        WP => SP
    use maria_constants_mod,    only:   &
        ZERO => S_ZERO,                 &
        ONE => S_ONE
    !-- Auxiliary subroutines --------------------------------------------------
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Arguments --------------------------------------------------------------
        integer,    intent(in   )               :: d                        !  1
        integer,    intent(in   ),  contiguous  :: n(:)                     !  2
        integer,    intent(in   ),  contiguous  :: mi(:)                    !  3
        integer,    intent(  out)               :: ierr                     !  4
        real(WP)                                :: stenval_sum
    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'STENVAL_SUM'
    !-- Variables --------------------------------------------------------------
        logical ::  bad_arg(4)
        integer ::  j

    !-- Default value ----------------------------------------------------------
        stenval_sum = ZERO

    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(1) = ( d < 1 )
        do j = 1, d
            bad_arg(2) = bad_arg(2) .or. ( n(j) < 0 )
        end do
        do j = 1, d
            bad_arg(3) = bad_arg(3) .or. ( mi(j) < 0 ) .or. ( mi(j) > n(j) )
        end do
        
    !-- Quick return if possible -----------------------------------------------
    quick: block
        if ( any(bad_arg) ) exit quick

        if ( any(n(1:d) == 0) ) return
    end block quick

    !-- Report bad input -------------------------------------------------------
        if ( any(bad_arg) ) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------    
        stenval_sum = ONE * sum(mi(1:d))
    end function stenval_sum

    !------------------------------------------------------------------------------------------------------------------------

    module function dtenval_sum &
    (d, n, mi, ierr)
    !-- Kinds, types and constants ---------------------------------------------
    use maria_kinds_mod,        only:   &
        WP => DP
    use maria_constants_mod,    only:   &
        ZERO => D_ZERO,                 &
        ONE => D_ONE
    !-- Auxiliary subroutines --------------------------------------------------
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Arguments --------------------------------------------------------------
        integer,    intent(in   )               :: d                        !  1
        integer,    intent(in   ),  contiguous  :: n(:)                     !  2
        integer,    intent(in   ),  contiguous  :: mi(:)                    !  3
        integer,    intent(  out)               :: ierr                     !  4
        real(WP)                                :: dtenval_sum
    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'DTENVAL_SUM'
    !-- Variables --------------------------------------------------------------
        logical ::  bad_arg(4)
        integer ::  j

    !-- Default value ----------------------------------------------------------
        dtenval_sum = ZERO

    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(1) = ( d < 1 )
        do j = 1, d
            bad_arg(2) = bad_arg(2) .or. ( n(j) < 0 )
        end do
        do j = 1, d
            bad_arg(3) = bad_arg(3) .or. ( mi(j) < 0 ) .or. ( mi(j) > n(j) )
        end do
        
    !-- Quick return if possible -----------------------------------------------
    quick: block
        if ( any(bad_arg) ) exit quick

        if ( any(n(1:d) == 0) ) return
    end block quick

    !-- Report bad input -------------------------------------------------------
        if ( any(bad_arg) ) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------    
        dtenval_sum = ONE * sum(mi(1:d))
    end function dtenval_sum

    !------------------------------------------------------------------------------------------------------------------------
end submodule maria_access_tensor_sub
