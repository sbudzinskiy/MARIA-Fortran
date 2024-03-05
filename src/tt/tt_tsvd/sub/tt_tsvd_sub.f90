!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the implementation of the [[maria_tt_tsvd_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Implements the [[maria_tt_tsvd_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
submodule (maria_tt_tsvd_mod) maria_tt_tsvd_sub
implicit none (type, external)

contains
    !------------------------------------------------------------------------------------------------------------------------

    module subroutine sttsvd &
    (ort, d, n, A, r, cores, work, lwork, iwork, liwork, ierr, maxr, rtolf, rerrf)
    !-- Kinds, types and constants ---------------------------------------------
    use maria_kinds_mod,        only:   &
        WP      =>  SP
    use maria_arr_mod,          only:   & 
        AR      =>  sarr
    use maria_constants_mod,    only:   &
        ZERO    =>  S_ZERO,             &
        ONE     =>  S_ONE,              &
        EPS     =>  S_MACHTOL
    !-- Auxiliary subroutines --------------------------------------------------
    use maria_reports_mod,      only:   &
        report_bad_arg
    use maria_comparison_mod,   only:   &
        safe_less
    use maria_utils_mod,        only:   &
        srealloc
    use maria_la_core_mod,      only:   &
        slacpy,                         &
        sroundup_lwork
    !-- Computational subroutines ----------------------------------------------
    use maria_la_core_mod,  only:   &
        sgesdd_q,                   &
        sdgmm
    use maria_lr_tsvd_mod,  only:   &
        schop
    !-- Arguments --------------------------------------------------------------
        character(1),   intent(in   )               :: ort                  !  1
        integer,        intent(in   )               :: d                    !  2
        integer,        intent(in   ),  contiguous  :: n(:)                 !  3
        real(WP),       intent(inout),  contiguous  :: A(:)                 !  4
        integer,        intent(  out),  contiguous  :: r(:)                 !  5
        type(AR),       intent(  out)               :: cores(:)             !  6
        real(WP),       intent(  out),  contiguous  :: work(:)              !  7
        integer,        intent(in   )               :: lwork                !  8
        integer,        intent(  out),  contiguous  :: iwork(:)             !  9
        integer,        intent(in   )               :: liwork               ! 10
        integer,        intent(  out)               :: ierr                 ! 11
    !-- Optional arguments -----------------------------------------------------
        integer,    intent(in   ),  contiguous, optional    :: maxr(:)      ! 12
        real(WP),   intent(in   ),              optional    :: rtolf        ! 13
        real(WP),   intent(  out),              optional    :: rerrf        ! 14
    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'STTSVD'
    !-- Variables --------------------------------------------------------------
        logical     ::  lw_query, left_ort, chop_rank, chop_rtolf, bad_arg(14)
        integer     ::  k, rl, rr, prodn, lsize, rsize, minsize, ifoo(1), &
                        req_lw, req_liw, i_U, mem_U, i_S, mem_S, i_VT, mem_VT, &
                        i_wrk, lwrk, lw_proc, liw_proc, i_actual_rerrf, mem_actual_rerrf
        real(WP)    ::  rtolf_cur, target_rerrf, foo(1)
    
    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(1) = ( ort /= 'l' ) .and. ( ort /= 'r' )
        bad_arg(2) = ( d < 2 )
        do k = 1, d
            bad_arg(3) = bad_arg(3) .or. ( n(k) < 0 )
        end do
        if ( present(maxr) ) then
            do k = 1, d-1
                bad_arg(12) = bad_arg(12) .or. ( maxr(k) < 0 )
            end do
        end if
        if ( present(rtolf) ) then
            bad_arg(13) = ( safe_less(rtolf, ZERO) )
        end if

        left_ort = ( ort == 'l' )
        chop_rank = present(maxr)
        chop_rtolf = present(rtolf)

    !-- Quick return if possible -----------------------------------------------
    quick: block
        if ( any(bad_arg) ) exit quick

        prodn = product(n(1:d))
        if ( prodn == 0 ) return
        
        if ( chop_rank ) then
            if ( any(maxr(1:d-1) == 0) ) then
                r = 0
                return
            end if
        end if
    end block quick

    !-- Estimate workspace -----------------------------------------------------
    workspace: block
        if ( any(bad_arg) ) exit workspace

        ! Storage: U, S, VT, actual_rerrf
        mem_U = 0
        mem_S = 0
        mem_VT = 0
        mem_actual_rerrf = d-1
        lw_proc = 0
        liw_proc = 0

        if ( left_ort ) then
            rl = 1
            rsize = prodn
            do k = 1, d-1
                lsize = rl * n(k)
                rsize = rsize / n(k)
                minsize = min(lsize, rsize)
                mem_U = max(mem_U, lsize * minsize)
                mem_S = max(mem_S, minsize)
                mem_VT = max(mem_VT, minsize * rsize)
                call sgesdd_q('s', lsize, rsize, foo, lsize, foo, foo, lsize, foo, minsize, foo, -1, ifoo, -1, ierr)
                lw_proc = max(lw_proc, int(foo(1)))
                liw_proc = max(liw_proc, ifoo(1))
                rl = minsize
            end do
        else
            lsize = prodn
            rr = 1
            do k = d, 2, -1
                lsize = lsize / n(k)
                rsize = rr * n(k)
                minsize = min(lsize, rsize)
                mem_U = max(mem_U, lsize * minsize)
                mem_S = max(mem_S, minsize)
                mem_VT = max(mem_VT, minsize * rsize)
                call sgesdd_q('s', lsize, rsize, foo, lsize, foo, foo, lsize, foo, minsize, foo, -1, ifoo, -1, ierr)
                lw_proc = max(lw_proc, int(foo(1)))
                liw_proc = max(liw_proc, ifoo(1))
                rr = minsize
            end do
        end if

        req_lw = mem_U + mem_S + mem_VT + mem_actual_rerrf + lw_proc
        req_liw = liw_proc

        lw_query = ( lwork == -1 ) .or. ( liwork == -1 )
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        bad_arg(8) = ( lwork < req_lw )
        bad_arg(10) = ( liwork < req_liw )
    end block workspace

    !-- Report bad input -------------------------------------------------------
        if ( any(bad_arg) ) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if ( chop_rtolf ) then
            rtolf_cur = rtolf
        else
            rtolf_cur = ZERO
        end if

        ! Slice workspace |...|...|....|..............|.....|
        !                   U   S   VT   actual_rerrf   wrk
        i_U = 1
        i_S = i_U + mem_U
        i_VT = i_S + mem_S
        i_actual_rerrf = i_VT + mem_VT
        i_wrk = i_actual_rerrf + mem_actual_rerrf
        lwrk = lwork - i_wrk + 1
        associate &
        (U => work(i_U:), S => work(i_S:), VT => work(i_VT:), actual_rerrf => work(i_actual_rerrf:), wrk => work(i_wrk:))
            ! Compute the approximation
            if ( left_ort ) then
                rl = 1
                rsize = prodn
                do k = 1, d-1
                    lsize = rl * n(k)
                    rsize = rsize / n(k)
                    minsize = min(lsize, rsize)
                    call sgesdd_q('s', lsize, rsize, A, lsize, S, U, lsize, VT, minsize, wrk, lwrk, iwork, liwork, ierr)

                    associate ( cur_rerrf => actual_rerrf(k) )
                        target_rerrf = rtolf_cur / sqrt(ONE*d-k)
                        if ( chop_rank ) then
                            r(k) = schop(minsize, S, ierr, maxr=maxr(k), rtolf=target_rerrf, rerrf=cur_rerrf)
                        else
                            r(k) = schop(minsize, S, ierr, rtolf=target_rerrf, rerrf=cur_rerrf)
                        end if
                        rtolf_cur = (rtolf_cur + cur_rerrf) / (1 - cur_rerrf**2) * (rtolf_cur - cur_rerrf)
                        rtolf_cur = sqrt(max(ZERO, rtolf_cur))
                    end associate
                    
                    call srealloc(.false., cores(k)%arr, lsize*r(k), ierr, factor=2)
                    call slacpy('a', lsize, r(k), U, lsize, cores(k)%arr, lsize)
                    call sdgmm('l', r(k), rsize, VT, minsize, S, 1, ierr)
                    if ( k < d-1 ) then
                        call slacpy('a', r(k), rsize, VT, minsize, A, r(k))
                    end if
                    rl = r(k)
                end do
                call srealloc(.false., cores(d)%arr, r(d-1)*n(d), ierr, factor=2)
                call slacpy('a', r(d-1), n(d), VT, minsize, cores(d)%arr, r(d-1))
            else
                lsize = prodn
                rr = 1
                do k = d, 2, -1
                    lsize = lsize / n(k)
                    rsize = rr * n(k)
                    minsize = min(lsize, rsize)
                    call sgesdd_q('s', lsize, rsize, A, lsize, S, U, lsize, VT, minsize, wrk, lwrk, iwork, liwork, ierr)

                    associate ( cur_rerrf => actual_rerrf(d-k+1) )
                        target_rerrf = rtolf_cur / sqrt(ONE*k-1)
                        if ( chop_rank ) then
                            r(k-1) = schop(minsize, S, ierr, maxr=maxr(k-1), rtolf=target_rerrf, rerrf=cur_rerrf)
                        else
                            r(k-1) = schop(minsize, S, ierr, rtolf=target_rerrf, rerrf=cur_rerrf)
                        end if
                        rtolf_cur = (rtolf_cur + cur_rerrf) / (1 - cur_rerrf**2) * (rtolf_cur - cur_rerrf)
                        rtolf_cur = sqrt(max(ZERO, rtolf_cur))
                    end associate
                    
                    call srealloc(.false., cores(k)%arr, r(k-1)*rsize, ierr, factor=2)
                    call slacpy('a', r(k-1), rsize, VT, minsize, cores(k)%arr, r(k-1))
                    call sdgmm('r', lsize, r(k-1), U, lsize, S, 1, ierr)
                    if ( k > 2 ) then
                        call slacpy('a', lsize, r(k-1), U, lsize, A, lsize)
                    end if
                    rr = r(k-1)
                end do
                call srealloc(.false., cores(1)%arr, n(1)*r(1), ierr, factor=2)
                call slacpy('a', n(1), r(1), U, n(1), cores(1)%arr, n(1))
            end if

            ! Sum the errors
            if ( present(rerrf) ) then
                rerrf = actual_rerrf(d-1)**2
                do k = d-2, 1, -1
                    rerrf = rerrf + actual_rerrf(k)**2 * (1 - rerrf)
                end do
                rerrf = sqrt(rerrf)
            end if
        end associate
    end subroutine sttsvd

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dttsvd &
    (ort, d, n, A, r, cores, work, lwork, iwork, liwork, ierr, maxr, rtolf, rerrf)
    !-- Kinds, types and constants ---------------------------------------------
    use maria_kinds_mod,        only:   &
        WP      =>  DP
    use maria_arr_mod,          only:   & 
        AR      =>  darr
    use maria_constants_mod,    only:   &
        ZERO    =>  D_ZERO,             &
        ONE     =>  D_ONE,              &
        EPS     =>  D_MACHTOL
    !-- Auxiliary subroutines --------------------------------------------------
    use maria_reports_mod,      only:   &
        report_bad_arg
    use maria_comparison_mod,   only:   &
        safe_less
    use maria_utils_mod,        only:   &
        drealloc
    use maria_la_core_mod,      only:   &
        dlacpy,                         &
        droundup_lwork
    !-- Computational subroutines ----------------------------------------------
    use maria_la_core_mod,  only:   &
        dgesdd_q,                   &
        ddgmm
    use maria_lr_tsvd_mod,  only:   &
        dchop
    !-- Arguments --------------------------------------------------------------
        character(1),   intent(in   )               :: ort                  !  1
        integer,        intent(in   )               :: d                    !  2
        integer,        intent(in   ),  contiguous  :: n(:)                 !  3
        real(WP),       intent(inout),  contiguous  :: A(:)                 !  4
        integer,        intent(  out),  contiguous  :: r(:)                 !  5
        type(AR),       intent(  out)               :: cores(:)             !  6
        real(WP),       intent(  out),  contiguous  :: work(:)              !  7
        integer,        intent(in   )               :: lwork                !  8
        integer,        intent(  out),  contiguous  :: iwork(:)             !  9
        integer,        intent(in   )               :: liwork               ! 10
        integer,        intent(  out)               :: ierr                 ! 11
    !-- Optional arguments -----------------------------------------------------
        integer,    intent(in   ),  contiguous, optional    :: maxr(:)      ! 12
        real(WP),   intent(in   ),              optional    :: rtolf        ! 13
        real(WP),   intent(  out),              optional    :: rerrf        ! 14
    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'DTTSVD'
    !-- Variables --------------------------------------------------------------
        logical     ::  lw_query, left_ort, chop_rank, chop_rtolf, bad_arg(14)
        integer     ::  k, rl, rr, prodn, lsize, rsize, minsize, ifoo(1), &
                        req_lw, req_liw, i_U, mem_U, i_S, mem_S, i_VT, mem_VT, &
                        i_wrk, lwrk, lw_proc, liw_proc, i_actual_rerrf, mem_actual_rerrf
        real(WP)    ::  rtolf_cur, target_rerrf, foo(1)
    
    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(1) = ( ort /= 'l' ) .and. ( ort /= 'r' )
        bad_arg(2) = ( d < 2 )
        do k = 1, d
            bad_arg(3) = bad_arg(3) .or. ( n(k) < 0 )
        end do
        if ( present(maxr) ) then
            do k = 1, d-1
                bad_arg(12) = bad_arg(12) .or. ( maxr(k) < 0 )
            end do
        end if
        if ( present(rtolf) ) then
            bad_arg(13) = ( safe_less(rtolf, ZERO) )
        end if

        left_ort = ( ort == 'l' )
        chop_rank = present(maxr)
        chop_rtolf = present(rtolf)

    !-- Quick return if possible -----------------------------------------------
    quick: block
        if ( any(bad_arg) ) exit quick

        prodn = product(n(1:d))
        if ( prodn == 0 ) return
        
        if ( chop_rank ) then
            if ( any(maxr(1:d-1) == 0) ) then
                r = 0
                return
            end if
        end if
    end block quick

    !-- Estimate workspace -----------------------------------------------------
    workspace: block
        if ( any(bad_arg) ) exit workspace

        ! Storage: U, S, VT, actual_rerrf
        mem_U = 0
        mem_S = 0
        mem_VT = 0
        mem_actual_rerrf = d-1
        lw_proc = 0
        liw_proc = 0

        if ( left_ort ) then
            rl = 1
            rsize = prodn
            do k = 1, d-1
                lsize = rl * n(k)
                rsize = rsize / n(k)
                minsize = min(lsize, rsize)
                mem_U = max(mem_U, lsize * minsize)
                mem_S = max(mem_S, minsize)
                mem_VT = max(mem_VT, minsize * rsize)
                call dgesdd_q('s', lsize, rsize, foo, lsize, foo, foo, lsize, foo, minsize, foo, -1, ifoo, -1, ierr)
                lw_proc = max(lw_proc, int(foo(1)))
                liw_proc = max(liw_proc, ifoo(1))
                rl = minsize
            end do
        else
            lsize = prodn
            rr = 1
            do k = d, 2, -1
                lsize = lsize / n(k)
                rsize = rr * n(k)
                minsize = min(lsize, rsize)
                mem_U = max(mem_U, lsize * minsize)
                mem_S = max(mem_S, minsize)
                mem_VT = max(mem_VT, minsize * rsize)
                call dgesdd_q('s', lsize, rsize, foo, lsize, foo, foo, lsize, foo, minsize, foo, -1, ifoo, -1, ierr)
                lw_proc = max(lw_proc, int(foo(1)))
                liw_proc = max(liw_proc, ifoo(1))
                rr = minsize
            end do
        end if

        req_lw = mem_U + mem_S + mem_VT + mem_actual_rerrf + lw_proc
        req_liw = liw_proc

        lw_query = ( lwork == -1 ) .or. ( liwork == -1 )
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        bad_arg(8) = ( lwork < req_lw )
        bad_arg(10) = ( liwork < req_liw )
    end block workspace

    !-- Report bad input -------------------------------------------------------
        if ( any(bad_arg) ) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if ( chop_rtolf ) then
            rtolf_cur = rtolf
        else
            rtolf_cur = ZERO
        end if

        ! Slice workspace |...|...|....|..............|.....|
        !                   U   S   VT   actual_rerrf   wrk
        i_U = 1
        i_S = i_U + mem_U
        i_VT = i_S + mem_S
        i_actual_rerrf = i_VT + mem_VT
        i_wrk = i_actual_rerrf + mem_actual_rerrf
        lwrk = lwork - i_wrk + 1
        associate &
        (U => work(i_U:), S => work(i_S:), VT => work(i_VT:), actual_rerrf => work(i_actual_rerrf:), wrk => work(i_wrk:))
            ! Compute the approximation
            if ( left_ort ) then
                rl = 1
                rsize = prodn
                do k = 1, d-1
                    lsize = rl * n(k)
                    rsize = rsize / n(k)
                    minsize = min(lsize, rsize)
                    call dgesdd_q('s', lsize, rsize, A, lsize, S, U, lsize, VT, minsize, wrk, lwrk, iwork, liwork, ierr)

                    associate ( cur_rerrf => actual_rerrf(k) )
                        target_rerrf = rtolf_cur / sqrt(ONE*d-k)
                        if ( chop_rank ) then
                            r(k) = dchop(minsize, S, ierr, maxr=maxr(k), rtolf=target_rerrf, rerrf=cur_rerrf)
                        else
                            r(k) = dchop(minsize, S, ierr, rtolf=target_rerrf, rerrf=cur_rerrf)
                        end if
                        rtolf_cur = (rtolf_cur + cur_rerrf) / (1 - cur_rerrf**2) * (rtolf_cur - cur_rerrf)
                        rtolf_cur = sqrt(max(ZERO, rtolf_cur))
                    end associate
                    
                    call drealloc(.false., cores(k)%arr, lsize*r(k), ierr, factor=2)
                    call dlacpy('a', lsize, r(k), U, lsize, cores(k)%arr, lsize)
                    call ddgmm('l', r(k), rsize, VT, minsize, S, 1, ierr)
                    if ( k < d-1 ) then
                        call dlacpy('a', r(k), rsize, VT, minsize, A, r(k))
                    end if
                    rl = r(k)
                end do
                call drealloc(.false., cores(d)%arr, r(d-1)*n(d), ierr, factor=2)
                call dlacpy('a', r(d-1), n(d), VT, minsize, cores(d)%arr, r(d-1))
            else
                lsize = prodn
                rr = 1
                do k = d, 2, -1
                    lsize = lsize / n(k)
                    rsize = rr * n(k)
                    minsize = min(lsize, rsize)
                    call dgesdd_q('s', lsize, rsize, A, lsize, S, U, lsize, VT, minsize, wrk, lwrk, iwork, liwork, ierr)

                    associate ( cur_rerrf => actual_rerrf(d-k+1) )
                        target_rerrf = rtolf_cur / sqrt(ONE*k-1)
                        if ( chop_rank ) then
                            r(k-1) = dchop(minsize, S, ierr, maxr=maxr(k-1), rtolf=target_rerrf, rerrf=cur_rerrf)
                        else
                            r(k-1) = dchop(minsize, S, ierr, rtolf=target_rerrf, rerrf=cur_rerrf)
                        end if
                        rtolf_cur = (rtolf_cur + cur_rerrf) / (1 - cur_rerrf**2) * (rtolf_cur - cur_rerrf)
                        rtolf_cur = sqrt(max(ZERO, rtolf_cur))
                    end associate
                    
                    call drealloc(.false., cores(k)%arr, r(k-1)*rsize, ierr, factor=2)
                    call dlacpy('a', r(k-1), rsize, VT, minsize, cores(k)%arr, r(k-1))
                    call ddgmm('r', lsize, r(k-1), U, lsize, S, 1, ierr)
                    if ( k > 2 ) then
                        call dlacpy('a', lsize, r(k-1), U, lsize, A, lsize)
                    end if
                    rr = r(k-1)
                end do
                call drealloc(.false., cores(1)%arr, n(1)*r(1), ierr, factor=2)
                call dlacpy('a', n(1), r(1), U, n(1), cores(1)%arr, n(1))
            end if

            ! Sum the errors
            if ( present(rerrf) ) then
                rerrf = actual_rerrf(d-1)**2
                do k = d-2, 1, -1
                    rerrf = rerrf + actual_rerrf(k)**2 * (1 - rerrf)
                end do
                rerrf = sqrt(rerrf)
            end if
        end associate
    end subroutine dttsvd

    !------------------------------------------------------------------------------------------------------------------------
end submodule maria_tt_tsvd_sub
