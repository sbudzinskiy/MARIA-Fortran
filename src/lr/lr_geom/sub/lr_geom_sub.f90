!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the implementation of the [[maria_lr_geom_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Implements the [[maria_lr_geom_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
submodule (maria_lr_geom_mod) maria_lr_geom_sub
implicit none (type, external)

contains
    !------------------------------------------------------------------------------------------------------------------------

    module subroutine slrproj_tangent &
    (m, n, r, U, ldu, VT, ldvt, B2AB, B2BA, pU, ldpu, pVT, ldpvt, UTAV, ldutav, ierr)
    !-- Kinds, types and constants ---------------------------------------------
    use maria_kinds_mod,        only:   &
        WP => SP
    use maria_constants_mod,    only:   &
        ZERO => S_ZERO,                 &
        ONE => S_ONE
    use maria_la_core_mod,      only:   &
        MM => smatmul
    !-- Auxiliary subroutines --------------------------------------------------
    use maria_reports_mod,    only: &
        report_bad_arg
    use maria_la_core_mod,  only:   &
        sroundup_lwork
    !-- Computational subroutines ----------------------------------------------
    use maria_la_core_mod,  only:   &
        sgemm
    !-- Arguments --------------------------------------------------------------
        integer,        intent(in   )               :: m                    !  1
        integer,        intent(in   )               :: n                    !  2
        integer,        intent(in   )               :: r                    !  3
        real(WP),       intent(in   ),  contiguous  :: U(:)                 !  4
        integer,        intent(in   )               :: ldu                  !  5
        real(WP),       intent(in   ),  contiguous  :: VT(:)                !  6
        integer,        intent(in   )               :: ldvt                 !  7
        procedure(MM),  intent(in   ),  pointer     :: B2AB                 !  8
        procedure(MM),  intent(in   ),  pointer     :: B2BA                 !  9
        real(WP),       intent(  out),  contiguous  :: pU(:)                ! 10
        integer,        intent(in   )               :: ldpu                 ! 11
        real(WP),       intent(  out),  contiguous  :: pVT(:)               ! 12
        integer,        intent(in   )               :: ldpvt                ! 13
        real(WP),       intent(  out),  contiguous  :: UTAV(:)              ! 14
        integer,        intent(in   )               :: ldutav               ! 15
        integer,        intent(  out)               :: ierr                 ! 16
    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'SLRPROJ_TANGENT'
    !-- Variables --------------------------------------------------------------
        logical ::  bad_arg(16)
    
    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(1) = ( m < 0 )
        bad_arg(2) = ( n < 0 )
        bad_arg(3) = ( r < 0 ) .or. ( r > min(m,n) )
        bad_arg(5) = ( ldu < max(1,m) )
        bad_arg(7) = ( ldvt < max(1,r) )
        bad_arg(11) = ( ldpu < max(1,m) )
        bad_arg(13) = ( ldpvt < max(1,r) )
        bad_arg(15) = ( ldutav < max(1,r) )

    !-- Quick return if possible -----------------------------------------------
    quick: block
        if ( any(bad_arg) ) exit quick

        if ( m == 0 .or. n == 0 .or. r == 0 ) return
    end block quick

    !-- Report bad input -------------------------------------------------------
        if ( any(bad_arg) ) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------    
        call B2BA('t', r, n, m, ONE, U, ldu, ZERO, pVT, ldpvt, ierr)
        call B2AB('t', m, r, n, ONE, VT, ldvt, ZERO, pU, ldpu, ierr)
        if ( m < n ) then
            call sgemm('t', 'n', r, r, m, ONE, U, ldu, pU, ldpu, ZERO, UTAV, ldutav)
        else
            call sgemm('n', 't', r, r, n, ONE, pVT, ldpvt, VT, ldvt, ZERO, UTAV, ldutav)
        end if
        call sgemm('n', 'n', m, r, r, -ONE, U, ldu, UTAV, ldutav, ONE, pU, ldpu)
        call sgemm('n', 'n', r, n, r, -ONE, UTAV, ldutav, VT, ldvt, ONE, pVT, ldpvt)
    end subroutine slrproj_tangent 

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dlrproj_tangent &
    (m, n, r, U, ldu, VT, ldvt, B2AB, B2BA, pU, ldpu, pVT, ldpvt, UTAV, ldutav, ierr)
    !-- Kinds, types and constants ---------------------------------------------
    use maria_kinds_mod,        only:   &
        WP => DP
    use maria_constants_mod,    only:   &
        ZERO => D_ZERO,                 &
        ONE => D_ONE
    use maria_la_core_mod,      only:   &
        MM => dmatmul
    !-- Auxiliary subroutines --------------------------------------------------
    use maria_reports_mod,    only: &
        report_bad_arg
    use maria_la_core_mod,  only:   &
        droundup_lwork
    !-- Computational subroutines ----------------------------------------------
    use maria_la_core_mod,  only:   &
        dgemm
    !-- Arguments --------------------------------------------------------------
        integer,        intent(in   )               :: m                    !  1
        integer,        intent(in   )               :: n                    !  2
        integer,        intent(in   )               :: r                    !  3
        real(WP),       intent(in   ),  contiguous  :: U(:)                 !  4
        integer,        intent(in   )               :: ldu                  !  5
        real(WP),       intent(in   ),  contiguous  :: VT(:)                !  6
        integer,        intent(in   )               :: ldvt                 !  7
        procedure(MM),  intent(in   ),  pointer     :: B2AB                 !  8
        procedure(MM),  intent(in   ),  pointer     :: B2BA                 !  9
        real(WP),       intent(  out),  contiguous  :: pU(:)                ! 10
        integer,        intent(in   )               :: ldpu                 ! 11
        real(WP),       intent(  out),  contiguous  :: pVT(:)               ! 12
        integer,        intent(in   )               :: ldpvt                ! 13
        real(WP),       intent(  out),  contiguous  :: UTAV(:)              ! 14
        integer,        intent(in   )               :: ldutav               ! 15
        integer,        intent(  out)               :: ierr                 ! 16
    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'DLRPROJ_TANGENT'
    !-- Variables --------------------------------------------------------------
        logical ::  bad_arg(16)
    
    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(1) = ( m < 0 )
        bad_arg(2) = ( n < 0 )
        bad_arg(3) = ( r < 0 ) .or. ( r > min(m,n) )
        bad_arg(5) = ( ldu < max(1,m) )
        bad_arg(7) = ( ldvt < max(1,r) )
        bad_arg(11) = ( ldpu < max(1,m) )
        bad_arg(13) = ( ldpvt < max(1,r) )
        bad_arg(15) = ( ldutav < max(1,r) )

    !-- Quick return if possible -----------------------------------------------
    quick: block
        if ( any(bad_arg) ) exit quick

        if ( m == 0 .or. n == 0 .or. r == 0 ) return
    end block quick

    !-- Report bad input -------------------------------------------------------
        if ( any(bad_arg) ) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------    
        call B2BA('t', r, n, m, ONE, U, ldu, ZERO, pVT, ldpvt, ierr)
        call B2AB('t', m, r, n, ONE, VT, ldvt, ZERO, pU, ldpu, ierr)
        if ( m < n ) then
            call dgemm('t', 'n', r, r, m, ONE, U, ldu, pU, ldpu, ZERO, UTAV, ldutav)
        else
            call dgemm('n', 't', r, r, n, ONE, pVT, ldpvt, VT, ldvt, ZERO, UTAV, ldutav)
        end if
        call dgemm('n', 'n', m, r, r, -ONE, U, ldu, UTAV, ldutav, ONE, pU, ldpu)
        call dgemm('n', 'n', r, n, r, -ONE, UTAV, ldutav, VT, ldvt, ONE, pVT, ldpvt)
    end subroutine dlrproj_tangent 

    !------------------------------------------------------------------------------------------------------------------------

    module function slrdotf_tangent &
    (m, n, r, pU1, ldpu1, pVT1, ldpvt1, C1, ldc1, pU2, ldpu2, pVT2, ldpvt2, C2, ldc2, ierr)
    !-- Kinds, types and constants ---------------------------------------------
    use maria_kinds_mod,        only:   &
        WP => SP
    use maria_constants_mod,    only:   &
        ZERO    =>  S_ZERO
    !-- Auxiliary subroutines --------------------------------------------------
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Computational subroutines ----------------------------------------------
    use maria_la_core_mod,  only:   &
        sgedotf
    !-- Arguments --------------------------------------------------------------
        integer,        intent(in   )               :: m                    !  1
        integer,        intent(in   )               :: n                    !  2
        integer,        intent(in   )               :: r                    !  3
        real(WP),       intent(in   ),  contiguous  :: pU1(:)               !  4
        integer,        intent(in   )               :: ldpu1                !  5
        real(WP),       intent(in   ),  contiguous  :: pVT1(:)              !  6
        integer,        intent(in   )               :: ldpvt1               !  7
        real(WP),       intent(in   ),  contiguous  :: C1(:)                !  8
        integer,        intent(in   )               :: ldc1                 !  9
        real(WP),       intent(in   ),  contiguous  :: pU2(:)               ! 10
        integer,        intent(in   )               :: ldpu2                ! 11
        real(WP),       intent(in   ),  contiguous  :: pVT2(:)              ! 12
        integer,        intent(in   )               :: ldpvt2               ! 13
        real(WP),       intent(in   ),  contiguous  :: C2(:)                ! 14
        integer,        intent(in   )               :: ldc2                 ! 15
        integer,        intent(  out)               :: ierr                 ! 16
        real(WP)                                    :: slrdotf_tangent
    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'SLRDOTF_TANGENT'
    !-- Variables --------------------------------------------------------------
        logical ::  bad_arg(16)
 
    !-- Default value ----------------------------------------------------------
        slrdotf_tangent = ZERO
   
    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(1) = ( m < 0 )
        bad_arg(2) = ( n < 0 )
        bad_arg(3) = ( r < 0 ) .or. ( r > min(m,n) )
        bad_arg(5) = ( ldpu1 < max(1,m) )
        bad_arg(7) = ( ldpvt1 < max(1,r) )
        bad_arg(9) = ( ldc1 < max(1,r) )
        bad_arg(11) = ( ldpu2 < max(1,m) )
        bad_arg(13) = ( ldpvt2 < max(1,r) )
        bad_arg(15) = ( ldc2 < max(1,r) )

    !-- Quick return if possible -----------------------------------------------
    quick: block
        if ( any(bad_arg) ) exit quick

        if ( m == 0 .or. n == 0 .or. r == 0 ) return
    end block quick

    !-- Report bad input -------------------------------------------------------
        if ( any(bad_arg) ) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------
        slrdotf_tangent = sgedotf('n', m, r, pU1, ldpu1, pU2, ldpu2, ierr)  &
            + sgedotf('n', r, n, pVT1, ldpvt1, pVT2, ldpvt2, ierr)          &
            + sgedotf('n', r, r, C1, ldc1, C2, ldc2, ierr)
    end function slrdotf_tangent 

    !------------------------------------------------------------------------------------------------------------------------

    module function dlrdotf_tangent &
    (m, n, r, pU1, ldpu1, pVT1, ldpvt1, C1, ldc1, pU2, ldpu2, pVT2, ldpvt2, C2, ldc2, ierr)
    !-- Kinds, types and constants ---------------------------------------------
    use maria_kinds_mod,        only:   &
        WP => DP
    use maria_constants_mod,    only:   &
        ZERO    =>  D_ZERO
    !-- Auxiliary subroutines --------------------------------------------------
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Computational subroutines ----------------------------------------------
    use maria_la_core_mod,  only:   &
        dgedotf
    !-- Arguments --------------------------------------------------------------
        integer,        intent(in   )               :: m                    !  1
        integer,        intent(in   )               :: n                    !  2
        integer,        intent(in   )               :: r                    !  3
        real(WP),       intent(in   ),  contiguous  :: pU1(:)               !  4
        integer,        intent(in   )               :: ldpu1                !  5
        real(WP),       intent(in   ),  contiguous  :: pVT1(:)              !  6
        integer,        intent(in   )               :: ldpvt1               !  7
        real(WP),       intent(in   ),  contiguous  :: C1(:)                !  8
        integer,        intent(in   )               :: ldc1                 !  9
        real(WP),       intent(in   ),  contiguous  :: pU2(:)               ! 10
        integer,        intent(in   )               :: ldpu2                ! 11
        real(WP),       intent(in   ),  contiguous  :: pVT2(:)              ! 12
        integer,        intent(in   )               :: ldpvt2               ! 13
        real(WP),       intent(in   ),  contiguous  :: C2(:)                ! 14
        integer,        intent(in   )               :: ldc2                 ! 15
        integer,        intent(  out)               :: ierr                 ! 16
        real(WP)                                    :: dlrdotf_tangent
    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'DLRDOTF_TANGENT'
    !-- Variables --------------------------------------------------------------
        logical ::  bad_arg(16)
 
    !-- Default value ----------------------------------------------------------
        dlrdotf_tangent = ZERO
   
    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(1) = ( m < 0 )
        bad_arg(2) = ( n < 0 )
        bad_arg(3) = ( r < 0 ) .or. ( r > min(m,n) )
        bad_arg(5) = ( ldpu1 < max(1,m) )
        bad_arg(7) = ( ldpvt1 < max(1,r) )
        bad_arg(9) = ( ldc1 < max(1,r) )
        bad_arg(11) = ( ldpu2 < max(1,m) )
        bad_arg(13) = ( ldpvt2 < max(1,r) )
        bad_arg(15) = ( ldc2 < max(1,r) )

    !-- Quick return if possible -----------------------------------------------
    quick: block
        if ( any(bad_arg) ) exit quick

        if ( m == 0 .or. n == 0 .or. r == 0 ) return
    end block quick

    !-- Report bad input -------------------------------------------------------
        if ( any(bad_arg) ) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------
        dlrdotf_tangent = dgedotf('n', m, r, pU1, ldpu1, pU2, ldpu2, ierr)  &
            + dgedotf('n', r, n, pVT1, ldpvt1, pVT2, ldpvt2, ierr)          &
            + dgedotf('n', r, r, C1, ldc1, C2, ldc2, ierr)
    end function dlrdotf_tangent 

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine slrretr_tangent &
    (m, n, r, S, U, ldu, VT, ldvt, alpha, pU, ldpu, pVT, ldpvt, C, ldc, work, lwork, iwork, liwork, ierr)
    !-- Kinds, types and constants ---------------------------------------------
    use maria_kinds_mod,        only:   &
        WP => SP
    use maria_constants_mod,    only:   &
        ZERO    =>  S_ZERO,             &
        ONE     =>  S_ONE
    !-- Auxiliary subroutines --------------------------------------------------
    use maria_reports_mod,  only:   &
        report_bad_arg
    use maria_la_core_mod,  only:   &
        slaset,                     &
        slacpy,                     &
        scopy,                      &
        sroundup_lwork
    !-- Computational subroutines ----------------------------------------------
    use maria_la_core_mod,  only:   &
        sgescal,                    &
        saxpy,                      &
        sgemm,                      &
        sgeqrf,                     &
        sorgqr,                     &
        sgelqf,                     &
        sorglq,                     &
        sgesdd_q
    !-- Arguments --------------------------------------------------------------
        integer,        intent(in   )               :: m                    !  1
        integer,        intent(in   )               :: n                    !  2
        integer,        intent(in   )               :: r                    !  3
        real(WP),       intent(inout),  contiguous  :: S(:)                 !  4
        real(WP),       intent(inout),  contiguous  :: U(:)                 !  5
        integer,        intent(in   )               :: ldu                  !  6
        real(WP),       intent(inout),  contiguous  :: VT(:)                !  7
        integer,        intent(in   )               :: ldvt                 !  8
        real(WP),       intent(in   )               :: alpha                !  9
        real(WP),       intent(in   ),  contiguous  :: pU(:)                ! 10
        integer,        intent(in   )               :: ldpu                 ! 11
        real(WP),       intent(in   ),  contiguous  :: pVT(:)               ! 12
        integer,        intent(in   )               :: ldpvt                ! 13
        real(WP),       intent(in   ),  contiguous  :: C(:)                 ! 14
        integer,        intent(in   )               :: ldc                  ! 15
        real(WP),       intent(  out),  contiguous  :: work(:)              ! 16
        integer,        intent(in   )               :: lwork                ! 17
        integer,        intent(  out),  contiguous  :: iwork(:)             ! 18
        integer,        intent(in   )               :: liwork               ! 19
        integer,        intent(  out)               :: ierr                 ! 20
    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'SLRRETR_TANGENT'
    !-- Variables --------------------------------------------------------------
        logical     ::  lw_query, bad_arg(20)
        integer     ::  ifoo(1), lsize_core, rsize_core, msize_core, mem_UU, i_UU,  &
                        mem_VVT, i_VVT, mem_core, i_core, mem_coreU, i_coreU,       &
                        mem_coreVT, i_coreVT, mem_coreS, i_coreS, mem_tau, i_tau,   &
                        i_wrk, lwrk, lw_proc, liw_proc, req_lw, req_liw
        real(WP)    ::  foo(1)
    
    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(1) = ( m < 0 )
        bad_arg(2) = ( n < 0 )
        bad_arg(3) = ( r < 0 ) .or. ( r > min(m,n) )
        bad_arg(6) = ( ldu < max(1,m) )
        bad_arg(8) = ( ldvt < max(1,r) )
        bad_arg(11) = ( ldpu < max(1,m) )
        bad_arg(13) = ( ldpvt < max(1,r) )
        bad_arg(15) = ( ldc < max(1,r) )

    !-- Quick return if possible -----------------------------------------------
    quick: block
        if ( any(bad_arg) ) exit quick

        if ( m == 0 .or. n == 0 .or. r == 0 ) return
    end block quick

    !-- Estimate workspace -----------------------------------------------------
    workspace: block
        if ( any(bad_arg) ) exit workspace

        ! Storage: UU, VVT, core, coreU, coreVT, coreS, tau
        lsize_core = min(m, 2*r)
        rsize_core = min(2*r, n)
        msize_core = min(lsize_core, rsize_core)

        mem_UU = m * (2*r)
        mem_VVT = (2*r) * n
        mem_core = lsize_core * rsize_core
        mem_coreU = lsize_core * msize_core
        mem_coreVT = msize_core * rsize_core
        mem_coreS = msize_core
        if ( m >= 2*r ) then
            mem_tau = r
        else
            mem_tau = m
        end if
        if ( n >= 2*r ) then
            mem_tau = max(mem_tau, r)
        else
            mem_tau = max(mem_tau, n)
        end if

        ! Query procedures
        if ( m >= 2*r ) then
            call sgeqrf(m, r, foo, m, foo, foo, -1, ierr)
            lw_proc = int(foo(1))
            call sorgqr(m, r, r, foo, m, foo, foo, -1, ierr)
            lw_proc = max(lw_proc, int(foo(1)))
        else
            !todo: QR with column pivoting and select only rank(pU) columns
            call sgeqrf(m, 2*r, foo, m, foo, foo, -1, ierr)
            lw_proc = int(foo(1))
            call sorgqr(m, m, m, foo, m, foo, foo, -1, ierr)
            lw_proc = max(lw_proc, int(foo(1)))
        end if

        if ( n >= 2*r ) then
            call sgelqf(r, n, foo, 2*r, foo, foo, -1, ierr)
            lw_proc = max(lw_proc, int(foo(1)))
            call sorglq(r, n, r, foo, 2*r, foo, foo, -1, ierr)
            lw_proc = max(lw_proc, int(foo(1)))
        else
            call sgelqf(2*r, n, foo, 2*r, foo, foo, -1, ierr)
            lw_proc = max(lw_proc, int(foo(1)))
            call sorglq(n, n, n, foo, 2*r, foo, foo, -1, ierr)
            lw_proc = max(lw_proc, int(foo(1)))
        end if

        call sgesdd_q('s', lsize_core, rsize_core, foo, 2*r, foo, foo, lsize_core, foo, msize_core, foo, -1, ifoo, -1, ierr)
        lw_proc = max(lw_proc, int(foo(1)))
        liw_proc = ifoo(1)

        req_lw = mem_UU + mem_VVT + mem_core + mem_coreU + mem_coreVT + mem_coreS + mem_tau  + lw_proc
        req_liw = liw_proc

        lw_query = ( lwork == -1 ) .or. ( liwork == -1 )
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        bad_arg(17) = ( lwork < req_lw )
        bad_arg(19) = ( liwork < req_liw )
    end block workspace

    !-- Report bad input -------------------------------------------------------
        if ( any(bad_arg) ) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------
        ! Slice workspace: |....|.....|.....|......|.......|.......|........|.....|
        !                    UU   VVT   tau   core   coreS   coreU   coreVT   wrk
        i_UU = 1
        i_VVT = i_UU + mem_UU
        i_tau = i_VVT + mem_VVT
        i_core = i_tau + mem_tau
        i_coreS = i_core + mem_core
        i_coreU = i_coreS + mem_coreS
        i_coreVT = i_coreU + mem_coreU
        i_wrk = i_coreVT + mem_coreVT
        lwrk = lwork - i_wrk + 1
    associate( UU => work(i_UU:), VVT => work(i_VVT:), tau => work(i_tau:), &
            core => work(i_core:), coreS => work(i_coreS:), coreU => work(i_coreU:), &
            coreVT => work(i_coreVT:), wrk => work(i_wrk:) )
        ! Set top left (without singular values) and bottom right blocks of the core
        call slaset('a', lsize_core, rsize_core, ZERO, ZERO, core, lsize_core)
        call slacpy('a', r, r, C, ldc, core, lsize_core)

        ! Form left orthogonal matrix and bottom left block of the core
        associate( rightUU => UU(1 + r*m:), botleft_core => core(1+r:) )
            call slacpy('a', m, r, U, ldu, UU, m)
            call slacpy('a', m, r, pU, ldpU, rightUU, m)
            if ( m >= 2*r ) then
                call sgeqrf(m, r, rightUU, m, tau, wrk, lwrk, ierr)
                call slacpy('u', r, r, rightUU, m, botleft_core, lsize_core)
                call sorgqr(m, r, r, rightUU, m, tau, wrk, lwrk, ierr)
            else
                call sgeqrf(m, 2*r, UU, m, tau, wrk, lwrk, ierr)
                associate( trapez => rightUU(1+r:) )
                    call slacpy('u', m-r, r, trapez, m, botleft_core, lsize_core)
                end associate
                call sorgqr(m, m, m, UU, m, tau, wrk, lwrk, ierr)
                call slacpy('a', m, r, U, ldu, UU, m)
            end if
        end associate

        ! Form right orthogonal matrix and top right block of the core
        associate( botVVT => VVT(1 + r:), topright_core => core(1+r*lsize_core:) )
            call slacpy('a', r, n, VT, ldvt, VVT, 2*r)
            call slacpy('a', r, n, pVT, ldpvt, botVVT, 2*r)
            if ( n >= 2*r ) then
                call sgelqf(r, n, botVVT, 2*r, tau, wrk, lwrk, ierr)
                call slacpy('l', r, r, botVVT, 2*r, topright_core, lsize_core)
                call sorglq(r, n, r, botVVT, 2*r, tau, wrk, lwrk, ierr)
            else
                call sgelqf(2*r, n, VVT, 2*r, tau, wrk, lwrk, ierr)
                associate( trapez => botVVT(1 + r*(2*r):) )
                    call slacpy('l', r, n-r, trapez, 2*r, topright_core, lsize_core)
                end associate
                call sorglq(n, n, n, VVT, 2*r, tau, wrk, lwrk, ierr)
            end if
        end associate

        ! Scale and add singular values
        call sgescal(lsize_core, rsize_core, alpha, core, lsize_core, ierr)
        call saxpy(r, ONE, S, 1, core, lsize_core+1)

        ! Update S, U, and VT
        call sgesdd_q('s', lsize_core, rsize_core, core, lsize_core, &
            coreS, coreU, lsize_core, coreVT, msize_core, wrk, lwrk, iwork, liwork, ierr)
        call scopy(r, coreS, 1, S, 1)
        call sgemm('n', 'n', m, r, lsize_core, ONE, UU, m, coreU, lsize_core, ZERO, U, ldu)
        call sgemm('n', 'n', r, n, rsize_core, ONE, coreVT, msize_core, VVT, 2*r, ZERO, VT, ldvt)
    end associate
    end subroutine slrretr_tangent

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dlrretr_tangent &
    (m, n, r, S, U, ldu, VT, ldvt, alpha, pU, ldpu, pVT, ldpvt, C, ldc, work, lwork, iwork, liwork, ierr)
    !-- Kinds, types and constants ---------------------------------------------
    use maria_kinds_mod,        only:   &
        WP => DP
    use maria_constants_mod,    only:   &
        ZERO    =>  D_ZERO,             &
        ONE     =>  D_ONE
    !-- Auxiliary subroutines --------------------------------------------------
    use maria_reports_mod,  only:   &
        report_bad_arg
    use maria_la_core_mod,  only:   &
        dlaset,                     &
        dlacpy,                     &
        dcopy,                      &
        droundup_lwork
    !-- Computational subroutines ----------------------------------------------
    use maria_la_core_mod,  only:   &
        dgescal,                    &
        daxpy,                      &
        dgemm,                      &
        dgeqrf,                     &
        dorgqr,                     &
        dgelqf,                     &
        dorglq,                     &
        dgesdd_q
    !-- Arguments --------------------------------------------------------------
        integer,        intent(in   )               :: m                    !  1
        integer,        intent(in   )               :: n                    !  2
        integer,        intent(in   )               :: r                    !  3
        real(WP),       intent(inout),  contiguous  :: S(:)                 !  4
        real(WP),       intent(inout),  contiguous  :: U(:)                 !  5
        integer,        intent(in   )               :: ldu                  !  6
        real(WP),       intent(inout),  contiguous  :: VT(:)                !  7
        integer,        intent(in   )               :: ldvt                 !  8
        real(WP),       intent(in   )               :: alpha                !  9
        real(WP),       intent(in   ),  contiguous  :: pU(:)                ! 10
        integer,        intent(in   )               :: ldpu                 ! 11
        real(WP),       intent(in   ),  contiguous  :: pVT(:)               ! 12
        integer,        intent(in   )               :: ldpvt                ! 13
        real(WP),       intent(in   ),  contiguous  :: C(:)                 ! 14
        integer,        intent(in   )               :: ldc                  ! 15
        real(WP),       intent(  out),  contiguous  :: work(:)              ! 16
        integer,        intent(in   )               :: lwork                ! 17
        integer,        intent(  out),  contiguous  :: iwork(:)             ! 18
        integer,        intent(in   )               :: liwork               ! 19
        integer,        intent(  out)               :: ierr                 ! 20
    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'DLRRETR_TANGENT'
    !-- Variables --------------------------------------------------------------
        logical     ::  lw_query, bad_arg(20)
        integer     ::  ifoo(1), lsize_core, rsize_core, msize_core, mem_UU, i_UU,  &
                        mem_VVT, i_VVT, mem_core, i_core, mem_coreU, i_coreU,       &
                        mem_coreVT, i_coreVT, mem_coreS, i_coreS, mem_tau, i_tau,   &
                        i_wrk, lwrk, lw_proc, liw_proc, req_lw, req_liw
        real(WP)    ::  foo(1)
    
    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(1) = ( m < 0 )
        bad_arg(2) = ( n < 0 )
        bad_arg(3) = ( r < 0 ) .or. ( r > min(m,n) )
        bad_arg(6) = ( ldu < max(1,m) )
        bad_arg(8) = ( ldvt < max(1,r) )
        bad_arg(11) = ( ldpu < max(1,m) )
        bad_arg(13) = ( ldpvt < max(1,r) )
        bad_arg(15) = ( ldc < max(1,r) )

    !-- Quick return if possible -----------------------------------------------
    quick: block
        if ( any(bad_arg) ) exit quick

        if ( m == 0 .or. n == 0 .or. r == 0 ) return
    end block quick

    !-- Estimate workspace -----------------------------------------------------
    workspace: block
        if ( any(bad_arg) ) exit workspace

        ! Storage: UU, VVT, core, coreU, coreVT, coreS, tau
        lsize_core = min(m, 2*r)
        rsize_core = min(2*r, n)
        msize_core = min(lsize_core, rsize_core)

        mem_UU = m * (2*r)
        mem_VVT = (2*r) * n
        mem_core = lsize_core * rsize_core
        mem_coreU = lsize_core * msize_core
        mem_coreVT = msize_core * rsize_core
        mem_coreS = msize_core
        if ( m >= 2*r ) then
            mem_tau = r
        else
            mem_tau = m
        end if
        if ( n >= 2*r ) then
            mem_tau = max(mem_tau, r)
        else
            mem_tau = max(mem_tau, n)
        end if

        ! Query procedures
        if ( m >= 2*r ) then
            call dgeqrf(m, r, foo, m, foo, foo, -1, ierr)
            lw_proc = int(foo(1))
            call dorgqr(m, r, r, foo, m, foo, foo, -1, ierr)
            lw_proc = max(lw_proc, int(foo(1)))
        else
            !todo: QR with column pivoting and select only rank(pU) columns
            call dgeqrf(m, 2*r, foo, m, foo, foo, -1, ierr)
            lw_proc = int(foo(1))
            call dorgqr(m, m, m, foo, m, foo, foo, -1, ierr)
            lw_proc = max(lw_proc, int(foo(1)))
        end if

        if ( n >= 2*r ) then
            call dgelqf(r, n, foo, 2*r, foo, foo, -1, ierr)
            lw_proc = max(lw_proc, int(foo(1)))
            call dorglq(r, n, r, foo, 2*r, foo, foo, -1, ierr)
            lw_proc = max(lw_proc, int(foo(1)))
        else
            call dgelqf(2*r, n, foo, 2*r, foo, foo, -1, ierr)
            lw_proc = max(lw_proc, int(foo(1)))
            call dorglq(n, n, n, foo, 2*r, foo, foo, -1, ierr)
            lw_proc = max(lw_proc, int(foo(1)))
        end if

        call dgesdd_q('s', lsize_core, rsize_core, foo, 2*r, foo, foo, lsize_core, foo, msize_core, foo, -1, ifoo, -1, ierr)
        lw_proc = max(lw_proc, int(foo(1)))
        liw_proc = ifoo(1)

        req_lw = mem_UU + mem_VVT + mem_core + mem_coreU + mem_coreVT + mem_coreS + mem_tau  + lw_proc
        req_liw = liw_proc

        lw_query = ( lwork == -1 ) .or. ( liwork == -1 )
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        bad_arg(17) = ( lwork < req_lw )
        bad_arg(19) = ( liwork < req_liw )
    end block workspace

    !-- Report bad input -------------------------------------------------------
        if ( any(bad_arg) ) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------
        ! Slice workspace: |....|.....|.....|......|.......|.......|........|.....|
        !                    UU   VVT   tau   core   coreS   coreU   coreVT   wrk
        i_UU = 1
        i_VVT = i_UU + mem_UU
        i_tau = i_VVT + mem_VVT
        i_core = i_tau + mem_tau
        i_coreS = i_core + mem_core
        i_coreU = i_coreS + mem_coreS
        i_coreVT = i_coreU + mem_coreU
        i_wrk = i_coreVT + mem_coreVT
        lwrk = lwork - i_wrk + 1
    associate( UU => work(i_UU:), VVT => work(i_VVT:), tau => work(i_tau:), &
            core => work(i_core:), coreS => work(i_coreS:), coreU => work(i_coreU:), &
            coreVT => work(i_coreVT:), wrk => work(i_wrk:) )
        ! Set top left (without singular values) and bottom right blocks of the core
        call dlaset('a', lsize_core, rsize_core, ZERO, ZERO, core, lsize_core)
        call dlacpy('a', r, r, C, ldc, core, lsize_core)

        ! Form left orthogonal matrix and bottom left block of the core
        associate( rightUU => UU(1 + r*m:), botleft_core => core(1+r:) )
            call dlacpy('a', m, r, U, ldu, UU, m)
            call dlacpy('a', m, r, pU, ldpU, rightUU, m)
            if ( m >= 2*r ) then
                call dgeqrf(m, r, rightUU, m, tau, wrk, lwrk, ierr)
                call dlacpy('u', r, r, rightUU, m, botleft_core, lsize_core)
                call dorgqr(m, r, r, rightUU, m, tau, wrk, lwrk, ierr)
            else
                call dgeqrf(m, 2*r, UU, m, tau, wrk, lwrk, ierr)
                associate( trapez => rightUU(1+r:) )
                    call dlacpy('u', m-r, r, trapez, m, botleft_core, lsize_core)
                end associate
                call dorgqr(m, m, m, UU, m, tau, wrk, lwrk, ierr)
                call dlacpy('a', m, r, U, ldu, UU, m)
            end if
        end associate

        ! Form right orthogonal matrix and top right block of the core
        associate( botVVT => VVT(1 + r:), topright_core => core(1+r*lsize_core:) )
            call dlacpy('a', r, n, VT, ldvt, VVT, 2*r)
            call dlacpy('a', r, n, pVT, ldpvt, botVVT, 2*r)
            if ( n >= 2*r ) then
                call dgelqf(r, n, botVVT, 2*r, tau, wrk, lwrk, ierr)
                call dlacpy('l', r, r, botVVT, 2*r, topright_core, lsize_core)
                call dorglq(r, n, r, botVVT, 2*r, tau, wrk, lwrk, ierr)
            else
                call dgelqf(2*r, n, VVT, 2*r, tau, wrk, lwrk, ierr)
                associate( trapez => botVVT(1 + r*(2*r):) )
                    call dlacpy('l', r, n-r, trapez, 2*r, topright_core, lsize_core)
                end associate
                call dorglq(n, n, n, VVT, 2*r, tau, wrk, lwrk, ierr)
            end if
        end associate

        ! Scale and add singular values
        call dgescal(lsize_core, rsize_core, alpha, core, lsize_core, ierr)
        call daxpy(r, ONE, S, 1, core, lsize_core+1)

        ! Update S, U, and VT
        call dgesdd_q('s', lsize_core, rsize_core, core, lsize_core, &
            coreS, coreU, lsize_core, coreVT, msize_core, wrk, lwrk, iwork, liwork, ierr)
        call dcopy(r, coreS, 1, S, 1)
        call dgemm('n', 'n', m, r, lsize_core, ONE, UU, m, coreU, lsize_core, ZERO, U, ldu)
        call dgemm('n', 'n', r, n, rsize_core, ONE, coreVT, msize_core, VVT, 2*r, ZERO, VT, ldvt)
    end associate
    end subroutine dlrretr_tangent

    !------------------------------------------------------------------------------------------------------------------------
end submodule maria_lr_geom_sub
