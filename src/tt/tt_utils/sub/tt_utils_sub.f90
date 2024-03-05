!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the implementation of the [[maria_tt_utils_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Implements the [[maria_tt_utils_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
submodule (maria_tt_utils_mod) maria_tt_utils_sub
implicit none (type, external)

contains
    !------------------------------------------------------------------------------------------------------------------------
    module subroutine ttrank &
    (d, r, k, rl, rr, ierr, r0, rd)
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS,                &
        BAD_IF_MORE
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer,    intent(in   )               :: d
        integer,    intent(in   ),  contiguous  :: r(:)
        integer,    intent(in   )               :: k
        integer,    intent(  out)               :: rl
        integer,    intent(  out)               :: rr
        integer,    intent(  out)               :: ierr
        integer,    intent(in   ),  optional    :: r0
        integer,    intent(in   ),  optional    :: rd

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'TTRANK'
    
    !-- Sanity check -----------------------------------------------------------
    sanity: block
        ierr = 0
        if (arg_is_bad(BAD_IF_LESS, d, 2)) then
            ierr = -1
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, k, 1)) then
            ierr = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, k, d)) then
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
        if ( k > 1 .and. k < d ) then
            rl = r(k-1)
            rr = r(k)
            return
        end if

        if ( k == 1 ) then
            rl = 1
            if (present(r0)) rl = r0
            rr = r(1)
            return
        end if

        if ( k == d ) then
            rl = r(d-1)
            rr = 1
            if (present(rd)) rr = rd
            return
        end if
    end subroutine ttrank

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine stt2full &
    (d, n, r, cores, A, work, lwork, ierr)
    !-- Kinds, types and constants ---------------------------------------------
    use maria_kinds_mod,        only:   &
        WP      =>  SP
    use maria_arr_mod,          only:   & 
        AR      =>  sarr
    use maria_constants_mod,    only:   &
        ZERO    =>  S_ZERO,             &
        ONE     =>  S_ONE
    !-- Auxiliary subroutines --------------------------------------------------
    use maria_reports_mod,  only:   &
        report_bad_arg
    use maria_la_core_mod,  only:   &
        scopy,                      &
        slaset,                     &
        sroundup_lwork
    !-- Computational subroutines ----------------------------------------------
    use maria_la_core_mod,  only:   &
        sgemm
    !-- Arguments --------------------------------------------------------------
        integer,    intent(in   )               :: d                        !  1
        integer,    intent(in   ),  contiguous  :: n(:)                     !  2
        integer,    intent(in   ),  contiguous  :: r(:)                     !  3
        type(AR),   intent(in   )               :: cores(:)                 !  4
        real(WP),   intent(  out),  contiguous  :: A(:)                     !  5
        real(WP),   intent(  out),  contiguous  :: work(:)                  !  6
        integer,    intent(in   )               :: lwork                    !  7
        integer,    intent(  out)               :: ierr                     !  8
    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'STT2FULL'
    !-- Variables --------------------------------------------------------------
        logical :: lw_query, bad_arg(8)
        integer :: k, rl, nk, rr, prodn, req_lw
    
    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(1) = ( d < 2 )
        do k = 1, d
            bad_arg(2) = bad_arg(2) .or. ( n(k) < 0 )
        end do
        do k = 1, d-1
            bad_arg(3) = bad_arg(3) .or. ( r(k) < 0 )
        end do

    !-- Quick return if possible -----------------------------------------------
    quick: block
        if ( any(bad_arg) ) exit quick

        prodn = product(n(1:d))
        if ( prodn == 0 ) return

        if ( any(r(1:d-1) == 0) ) then
            call slaset('a', 1, prodn, ZERO, ZERO, A, 1)
            return
        end if
    end block quick

    !-- Estimate workspace -----------------------------------------------------
    workspace: block
        if ( any(bad_arg) ) exit workspace

        if ( mod(d,2) == 0 ) then
            prodn = n(d)
            req_lw = r(d-1) * prodn
            do k = d-2, 2, -2
                prodn = prodn * n(k) * n(k+1)
                req_lw = max(req_lw, r(k-1) * prodn)
            end do
        else
            prodn = 1
            req_lw = prodn
            do k = d-1, 2, -2
                prodn = prodn * n(k) * n(k+1)
                req_lw = max(req_lw, r(k-1) * prodn)
            end do
        end if

        lw_query = (lwork == -1)
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            return
        end if

        bad_arg(7) = ( lwork < req_lw )
    end block workspace

    !-- Report bad input -------------------------------------------------------
        if ( any(bad_arg) ) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------
        prodn = n(d)
        if ( mod(d,2) == 0 ) then
            call scopy(r(d-1)*n(d), cores(d)%arr, 1, work, 1)
        else
            call scopy(r(d-1)*n(d), cores(d)%arr, 1, A, 1)
        end if
        do k = d-1, 2, -1
            rl = r(k-1)
            nk = n(k)
            rr = r(k)
            if ( mod(k,2) == 0 ) then
                call sgemm('n', 'n', rl*nk, prodn, rr, ONE, cores(k)%arr, rl*nk, A, rr, ZERO, work, rl*nk)
            else
                call sgemm('n', 'n', rl*nk, prodn, rr, ONE, cores(k)%arr, rl*nk, work, rr, ZERO, A, rl*nk)
            end if
            prodn = prodn * nk
        end do
        nk = n(1)
        rr = r(1)
        call sgemm('n', 'n', nk, prodn, rr, ONE, cores(1)%arr, nk, work, rr, ZERO, A, nk)
    end subroutine stt2full

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dtt2full &
    (d, n, r, cores, A, work, lwork, ierr)
    !-- Kinds, types and constants ---------------------------------------------
    use maria_kinds_mod,        only:   &
        WP      =>  DP
    use maria_arr_mod,          only:   & 
        AR      =>  darr
    use maria_constants_mod,    only:   &
        ZERO    =>  D_ZERO,             &
        ONE     =>  D_ONE
    !-- Auxiliary subroutines --------------------------------------------------
    use maria_reports_mod,  only:   &
        report_bad_arg
    use maria_la_core_mod,  only:   &
        dcopy,                      &
        dlaset,                     &
        droundup_lwork
    !-- Computational subroutines ----------------------------------------------
    use maria_la_core_mod,  only:   &
        dgemm
    !-- Arguments --------------------------------------------------------------
        integer,    intent(in   )               :: d                        !  1
        integer,    intent(in   ),  contiguous  :: n(:)                     !  2
        integer,    intent(in   ),  contiguous  :: r(:)                     !  3
        type(AR),   intent(in   )               :: cores(:)                 !  4
        real(WP),   intent(  out),  contiguous  :: A(:)                     !  5
        real(WP),   intent(  out),  contiguous  :: work(:)                  !  6
        integer,    intent(in   )               :: lwork                    !  7
        integer,    intent(  out)               :: ierr                     !  8
    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'DTT2FULL'
    !-- Variables --------------------------------------------------------------
        logical :: lw_query, bad_arg(8)
        integer :: k, rl, nk, rr, prodn, req_lw
    
    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(1) = ( d < 2 )
        do k = 1, d
            bad_arg(2) = bad_arg(2) .or. ( n(k) < 0 )
        end do
        do k = 1, d-1
            bad_arg(3) = bad_arg(3) .or. ( r(k) < 0 )
        end do

    !-- Quick return if possible -----------------------------------------------
    quick: block
        if ( any(bad_arg) ) exit quick

        prodn = product(n(1:d))
        if ( prodn == 0 ) return

        if ( any(r(1:d-1) == 0) ) then
            call dlaset('a', 1, prodn, ZERO, ZERO, A, 1)
            return
        end if
    end block quick

    !-- Estimate workspace -----------------------------------------------------
    workspace: block
        if ( any(bad_arg) ) exit workspace

        if ( mod(d,2) == 0 ) then
            prodn = n(d)
            req_lw = r(d-1) * prodn
            do k = d-2, 2, -2
                prodn = prodn * n(k) * n(k+1)
                req_lw = max(req_lw, r(k-1) * prodn)
            end do
        else
            prodn = 1
            req_lw = prodn
            do k = d-1, 2, -2
                prodn = prodn * n(k) * n(k+1)
                req_lw = max(req_lw, r(k-1) * prodn)
            end do
        end if

        lw_query = (lwork == -1)
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            return
        end if

        bad_arg(7) = ( lwork < req_lw )
    end block workspace

    !-- Report bad input -------------------------------------------------------
        if ( any(bad_arg) ) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------
        prodn = n(d)
        if ( mod(d,2) == 0 ) then
            call dcopy(r(d-1)*n(d), cores(d)%arr, 1, work, 1)
        else
            call dcopy(r(d-1)*n(d), cores(d)%arr, 1, A, 1)
        end if
        do k = d-1, 2, -1
            rl = r(k-1)
            nk = n(k)
            rr = r(k)
            if ( mod(k,2) == 0 ) then
                call dgemm('n', 'n', rl*nk, prodn, rr, ONE, cores(k)%arr, rl*nk, A, rr, ZERO, work, rl*nk)
            else
                call dgemm('n', 'n', rl*nk, prodn, rr, ONE, cores(k)%arr, rl*nk, work, rr, ZERO, A, rl*nk)
            end if
            prodn = prodn * nk
        end do
        nk = n(1)
        rr = r(1)
        call dgemm('n', 'n', nk, prodn, rr, ONE, cores(1)%arr, nk, work, rr, ZERO, A, nk)
    end subroutine dtt2full

    !------------------------------------------------------------------------------------------------------------------------
end submodule maria_tt_utils_sub
