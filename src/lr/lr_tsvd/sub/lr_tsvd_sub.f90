!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the implementation of the [[maria_lr_tsvd_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Implements the [[maria_lr_tsvd_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
submodule (maria_lr_tsvd_mod) maria_lr_tsvd_sub
implicit none (type, external)

contains
    !------------------------------------------------------------------------------------------------------------------------

    module function schop &
    (n, S, info, maxr, rtolf, atolf, rtol2, atol2, rerrf, aerrf, rerr2, aerr2)
    use maria_kinds_mod,      only: &
        WP => SP
    use maria_constants_mod,  only: &
        ZERO => S_ZERO
    use maria_comparison_mod, only: &
        safe_eq,                    &
        safe_less
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)              :: n
        real(WP), intent(in),  contiguous :: S(:)
        integer,  intent(out)             :: info
        integer,  intent(in),  optional   :: maxr
        real(WP), intent(in),  optional   :: rtolf
        real(WP), intent(in),  optional   :: atolf
        real(WP), intent(in),  optional   :: rtol2
        real(WP), intent(in),  optional   :: atol2
        real(WP), intent(out), optional   :: rerrf
        real(WP), intent(out), optional   :: aerrf
        real(WP), intent(out), optional   :: rerr2
        real(WP), intent(out), optional   :: aerr2
        integer                           :: schop

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SCHOP'
        logical                 :: chopmr, choprf, chopaf, chopr2, chopa2
        integer                 :: r, i
        real(WP)                :: summ, errfsq, err2, newerrfsq, newerr2
    
    !-- Default values ---------------------------------------------------------
        schop = 0
        chopmr = .false.
        choprf = .false.
        chopaf = .false.
        chopr2 = .false.
        chopa2 = .false.
        if (present(rerrf)) rerrf = ZERO
        if (present(aerrf)) aerrf = ZERO
        if (present(rerr2)) rerr2 = ZERO
        if (present(aerr2)) aerr2 = ZERO

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -1
            exit sanity
        end if
        if (present(maxr)) then
            chopmr = .true.
            if (arg_is_bad(BAD_IF_LESS, maxr, 0)) then
                info = -4
                exit sanity
            end if
        end if
        if (present(rtolf)) then
            choprf = .true.
            if (arg_is_bad(BAD_IF_LESS, rtolf, ZERO)) then
                info = -5
                exit sanity
            end if
        end if
        if (present(atolf)) then
            chopaf = .true.
            if (arg_is_bad(BAD_IF_LESS, atolf, ZERO)) then
                info = -6
                exit sanity
            end if
        end if
        if (present(rtol2)) then
            chopr2 = .true.
            if (arg_is_bad(BAD_IF_LESS, rtol2, ZERO)) then
                info = -7
                exit sanity
            end if
        end if
        if (present(atol2)) then
            chopa2 = .true.
            if (arg_is_bad(BAD_IF_LESS, atol2, ZERO)) then
                info = -8
                exit sanity
            end if
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (n == 0) return
        if (safe_eq(S(1), ZERO)) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        r = n
        errfsq = ZERO
        err2 = ZERO

        if (chopmr) r = min(r, maxr)
        do i = n, r+1, -1
            errfsq = errfsq + S(i)**2
            err2 = S(i)
        end do

        if (choprf .or. present(rerrf)) then
            summ = sum(S(1:n)**2)
        else
            summ = ZERO
        end if

        if (choprf .or. chopaf .or. chopr2 .or. chopa2) then
            reduce_rank: do
                if (r == 0) exit reduce_rank
                
                newerrfsq = errfsq + S(r)**2
                newerr2 = S(r)

                if (chopa2) then
                    if (safe_less(atol2, newerr2)) exit reduce_rank
                end if
                if (chopr2) then
                    if (safe_less(rtol2 * S(1), newerr2)) exit reduce_rank
                end if
                if (chopaf) then
                    if (safe_less(atolf, sqrt(newerrfsq))) exit reduce_rank
                end if
                if (choprf) then
                    if (safe_less(rtolf, sqrt(newerrfsq / summ))) exit reduce_rank
                end if
                
                errfsq = newerrfsq
                err2 = newerr2
                r = r - 1
            end do reduce_rank            
        end if

        schop = r
        if (present(rerrf)) rerrf = sqrt(errfsq / summ)
        if (present(rerr2)) rerr2 = err2 / S(1)
        if (present(aerrf)) aerrf = sqrt(errfsq)
        if (present(aerr2)) aerr2 = err2
    end function schop

    !------------------------------------------------------------------------------------------------------------------------

    module function dchop &
    (n, S, info, maxr, rtolf, atolf, rtol2, atol2, rerrf, aerrf, rerr2, aerr2)
    use maria_kinds_mod,      only: &
        WP => DP
    use maria_constants_mod,  only: &
        ZERO => D_ZERO
    use maria_comparison_mod, only: &
        safe_eq,                    &
        safe_less
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)              :: n
        real(WP), intent(in),  contiguous :: S(:)
        integer,  intent(out)             :: info
        integer,  intent(in),  optional   :: maxr
        real(WP), intent(in),  optional   :: rtolf
        real(WP), intent(in),  optional   :: atolf
        real(WP), intent(in),  optional   :: rtol2
        real(WP), intent(in),  optional   :: atol2
        real(WP), intent(out), optional   :: rerrf
        real(WP), intent(out), optional   :: aerrf
        real(WP), intent(out), optional   :: rerr2
        real(WP), intent(out), optional   :: aerr2
        integer                           :: dchop

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DCHOP'
        logical                 :: chopmr, choprf, chopaf, chopr2, chopa2
        integer                 :: r, i
        real(WP)                :: summ, errfsq, err2, newerrfsq, newerr2
    
    !-- Default values ---------------------------------------------------------
        dchop = 0
        chopmr = .false.
        choprf = .false.
        chopaf = .false.
        chopr2 = .false.
        chopa2 = .false.
        if (present(rerrf)) rerrf = ZERO
        if (present(aerrf)) aerrf = ZERO
        if (present(rerr2)) rerr2 = ZERO
        if (present(aerr2)) aerr2 = ZERO

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -1
            exit sanity
        end if
        if (present(maxr)) then
            chopmr = .true.
            if (arg_is_bad(BAD_IF_LESS, maxr, 0)) then
                info = -4
                exit sanity
            end if
        end if
        if (present(rtolf)) then
            choprf = .true.
            if (arg_is_bad(BAD_IF_LESS, rtolf, ZERO)) then
                info = -5
                exit sanity
            end if
        end if
        if (present(atolf)) then
            chopaf = .true.
            if (arg_is_bad(BAD_IF_LESS, atolf, ZERO)) then
                info = -6
                exit sanity
            end if
        end if
        if (present(rtol2)) then
            chopr2 = .true.
            if (arg_is_bad(BAD_IF_LESS, rtol2, ZERO)) then
                info = -7
                exit sanity
            end if
        end if
        if (present(atol2)) then
            chopa2 = .true.
            if (arg_is_bad(BAD_IF_LESS, atol2, ZERO)) then
                info = -8
                exit sanity
            end if
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (n == 0) return
        if (safe_eq(S(1), ZERO)) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        r = n
        errfsq = ZERO
        err2 = ZERO

        if (chopmr) r = min(r, maxr)
        do i = n, r+1, -1
            errfsq = errfsq + S(i)**2
            err2 = S(i)
        end do

        if (choprf .or. present(rerrf)) then
            summ = sum(S(1:n)**2)
        else
            summ = ZERO
        end if

        if (choprf .or. chopaf .or. chopr2 .or. chopa2) then
            reduce_rank: do
                if (r == 0) exit reduce_rank
                
                newerrfsq = errfsq + S(r)**2
                newerr2 = S(r)

                if (chopa2) then
                    if (safe_less(atol2, newerr2)) exit reduce_rank
                end if
                if (chopr2) then
                    if (safe_less(rtol2 * S(1), newerr2)) exit reduce_rank
                end if
                if (chopaf) then
                    if (safe_less(atolf, sqrt(newerrfsq))) exit reduce_rank
                end if
                if (choprf) then
                    if (safe_less(rtolf, sqrt(newerrfsq / summ))) exit reduce_rank
                end if
                
                errfsq = newerrfsq
                err2 = newerr2
                r = r - 1
            end do reduce_rank            
        end if

        dchop = r
        if (present(rerrf)) rerrf = sqrt(errfsq / summ)
        if (present(rerr2)) rerr2 = err2 / S(1)
        if (present(aerrf)) aerrf = sqrt(errfsq)
        if (present(aerr2)) aerr2 = err2
    end function dchop

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine sgersvd1 &
    (m, n, B2AB, B2BA, kc, csk, ldcsk, niter_power, S, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,      only: &
        WP => SP
    use maria_constants_mod,  only: &
        ZERO => S_ZERO,             &
        ONE => S_ONE
    use maria_la_core_mod,    only: &
        MM => smatmul,              &
        sgeqrf,                     &
        sorgqr,                     &
        sgelqf,                     &
        sorglq,                     &
        sgesdd_q,                   &
        sgemm,                      &
        sroundup_lwork
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS,                &
        BAD_IF_MORE
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer,       intent(in)                :: m
        integer,       intent(in)                :: n
        procedure(MM), intent(in),    pointer    :: B2AB
        procedure(MM), intent(in),    pointer    :: B2BA
        integer,       intent(in)                :: kc
        real(WP),      intent(inout), contiguous :: csk(:)
        integer,       intent(in)                :: ldcsk
        integer,       intent(in)                :: niter_power
        real(WP),      intent(out),   contiguous :: S(:)
        real(WP),      intent(out),   contiguous :: U(:)
        integer,       intent(in)                :: ldu
        real(WP),      intent(out),   contiguous :: VT(:)
        integer,       intent(in)                :: ldvt
        real(WP),      intent(out),   contiguous :: work(:)
        integer,       intent(in)                :: lwork
        integer,       intent(out),   contiguous :: iwork(:)
        integer,       intent(in)                :: liwork
        integer,       intent(out)               :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SGERSVD1'
        logical                 :: lw_query
        integer                 :: req_lw, req_liw, req_lw1, req_lw2, lw_proc, liw_proc, &
                                   mem_QTA, i_QTA, mem_tau, i_tau, mem_tmpU, i_tmpU, &
                                   lwrk, i_wrk, ifoo(1), i
        real(WP)                :: foo(1) 
    
    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, m, 0)) then
            info = -1
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, kc, 0)) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, kc, n)) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldcsk, max(1,m))) then
            info = -7
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, niter_power, 0)) then
            info = -8
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -11
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,kc))) then
            info = -13
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. kc == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) exit wrkwrk

        !-- 1. Storage: QTA, tau
        mem_QTA = kc * n
        mem_tau = kc
        !-- 1. Query procedures
        call sgeqrf(m, kc, foo, ldcsk, foo, foo, -1, info)
        lw_proc = int(foo(1))
        call sorgqr(m, kc, kc, foo, ldcsk, foo, foo, -1, info)
        lw_proc = max(lw_proc, int(foo(1)))
        call sgelqf(kc, n, foo, kc, foo, foo, -1, info)
        lw_proc = max(lw_proc, int(foo(1)))
        call sorglq(kc, n, kc, foo, kc, foo, foo, -1, info)
        lw_proc = max(lw_proc, int(foo(1)))
            
        req_lw1 = lw_proc + mem_QTA + mem_tau

        !-- 2. Storage: QTA, tmpU
        mem_QTA = kc * n
        mem_tmpU = kc * kc
        !-- 2. Query procedures
        call sgesdd_q('s', kc, n, foo, kc, foo, foo, kc, foo, ldvt, foo, -1, ifoo, -1, info)
        lw_proc = int(foo(1))
        liw_proc = ifoo(1)

        req_lw2 = lw_proc + mem_QTA + mem_tmpU

        !-- Total
        req_lw = max(req_lw1, req_lw2)
        req_liw = liw_proc

        lw_query = (lwork == -1 .or. liwork == -1)
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -15
            exit wrkwrk
        end if
        if (arg_is_bad(BAD_IF_LESS, liwork, req_liw)) then
            info = -17
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------

    !-- 1. |..QTA..|..tau..|..wrk..|
        i_QTA = 1
        i_tau = i_QTA + mem_QTA
        i_wrk = i_tau + mem_tau
        lwrk = lwork - i_wrk + 1
    associate &
    (QTA => work(i_QTA:), tau => work(i_tau:), wrk => work(i_wrk:))
        call sgeqrf(m, kc, csk, ldcsk, tau, wrk, lwrk, info)
        call sorgqr(m, kc, kc, csk, ldcsk, tau, wrk, lwrk, info)
        do i = 1, niter_power
            call B2BA('t', kc, n, m, ONE, csk, ldcsk, ZERO, QTA, kc, info)
            call sgelqf(kc, n, QTA, kc, tau, wrk, lwrk, info)
            call sorglq(kc, n, kc, QTA, kc, tau, wrk, lwrk, info)
            call B2AB('t', m, kc, n, ONE, QTA, kc, ZERO, csk, ldcsk, info)
            call sgeqrf(m, kc, csk, ldcsk, tau, wrk, lwrk, info)
            call sorgqr(m, kc, kc, csk, ldcsk, tau, wrk, lwrk, info)
        end do
    end associate

    !-- 2. |..QTA..|..tmpU..|..wrk..|
        i_QTA = 1
        i_tmpU = i_QTA + mem_QTA
        i_wrk = i_tmpU + mem_tmpU
        lwrk = lwork - i_wrk + 1
    associate &
    (QTA => work(i_QTA:), tmpU => work(i_tmpU:), wrk => work(i_wrk:))
        call B2BA('t', kc, n, m, ONE, csk, ldcsk, ZERO, QTA, kc, info)
        call sgesdd_q('s', kc, n, QTA, kc, S, tmpU, kc, VT, ldvt, wrk, lwrk, iwork, liwork, info)
        call sgemm('n', 'n', m, kc, kc, ONE, csk, ldcsk, tmpU, kc, ZERO, U, ldu)
    end associate
    end subroutine sgersvd1

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dgersvd1 &
    (m, n, B2AB, B2BA, kc, csk, ldcsk, niter_power, S, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,      only: &
        WP => DP
    use maria_constants_mod,  only: &
        ZERO => D_ZERO,             &
        ONE => D_ONE
    use maria_la_core_mod,    only: &
        MM => Dmatmul,              &
        dgeqrf,                     &
        dorgqr,                     &
        dgelqf,                     &
        dorglq,                     &
        dgesdd_q,                   &
        dgemm,                      &
        droundup_lwork
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS,                &
        BAD_IF_MORE
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer,       intent(in)                :: m
        integer,       intent(in)                :: n
        procedure(MM), intent(in),    pointer    :: B2AB
        procedure(MM), intent(in),    pointer    :: B2BA
        integer,       intent(in)                :: kc
        real(WP),      intent(inout), contiguous :: csk(:)
        integer,       intent(in)                :: ldcsk
        integer,       intent(in)                :: niter_power
        real(WP),      intent(out),   contiguous :: S(:)
        real(WP),      intent(out),   contiguous :: U(:)
        integer,       intent(in)                :: ldu
        real(WP),      intent(out),   contiguous :: VT(:)
        integer,       intent(in)                :: ldvt
        real(WP),      intent(out),   contiguous :: work(:)
        integer,       intent(in)                :: lwork
        integer,       intent(out),   contiguous :: iwork(:)
        integer,       intent(in)                :: liwork
        integer,       intent(out)               :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DGERSVD1'
        logical                 :: lw_query
        integer                 :: req_lw, req_liw, req_lw1, req_lw2, lw_proc, liw_proc, &
                                   mem_QTA, i_QTA, mem_tau, i_tau, mem_tmpU, i_tmpU, &
                                   lwrk, i_wrk, ifoo(1), i
        real(WP)                :: foo(1) 
    
    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, m, 0)) then
            info = -1
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, kc, 0)) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, kc, n)) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldcsk, max(1,m))) then
            info = -7
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, niter_power, 0)) then
            info = -8
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -11
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,kc))) then
            info = -13
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. kc == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) exit wrkwrk

        !-- 1. Storage: QTA, tau
        mem_QTA = kc * n
        mem_tau = kc
        !-- 1. Query procedures
        call dgeqrf(m, kc, foo, ldcsk, foo, foo, -1, info)
        lw_proc = int(foo(1))
        call dorgqr(m, kc, kc, foo, ldcsk, foo, foo, -1, info)
        lw_proc = max(lw_proc, int(foo(1)))
        call dgelqf(kc, n, foo, kc, foo, foo, -1, info)
        lw_proc = max(lw_proc, int(foo(1)))
        call dorglq(kc, n, kc, foo, kc, foo, foo, -1, info)
        lw_proc = max(lw_proc, int(foo(1)))
            
        req_lw1 = lw_proc + mem_QTA + mem_tau

        !-- 2. Storage: QTA, tmpU
        mem_QTA = kc * n
        mem_tmpU = kc * kc
        !-- 2. Query procedures
        call dgesdd_q('s', kc, n, foo, kc, foo, foo, kc, foo, ldvt, foo, -1, ifoo, -1, info)
        lw_proc = int(foo(1))
        liw_proc = ifoo(1)

        req_lw2 = lw_proc + mem_QTA + mem_tmpU

        !-- Total
        req_lw = max(req_lw1, req_lw2)
        req_liw = liw_proc

        lw_query = (lwork == -1 .or. liwork == -1)
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -15
            exit wrkwrk
        end if
        if (arg_is_bad(BAD_IF_LESS, liwork, req_liw)) then
            info = -17
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------

    !-- 1. |..QTA..|..tau..|..wrk..|
        i_QTA = 1
        i_tau = i_QTA + mem_QTA
        i_wrk = i_tau + mem_tau
        lwrk = lwork - i_wrk + 1
    associate &
    (QTA => work(i_QTA:), tau => work(i_tau:), wrk => work(i_wrk:))
        call dgeqrf(m, kc, csk, ldcsk, tau, wrk, lwrk, info)
        call dorgqr(m, kc, kc, csk, ldcsk, tau, wrk, lwrk, info)
        do i = 1, niter_power
            call B2BA('t', kc, n, m, ONE, csk, ldcsk, ZERO, QTA, kc, info)
            call dgelqf(kc, n, QTA, kc, tau, wrk, lwrk, info)
            call dorglq(kc, n, kc, QTA, kc, tau, wrk, lwrk, info)
            call B2AB('t', m, kc, n, ONE, QTA, kc, ZERO, csk, ldcsk, info)
            call dgeqrf(m, kc, csk, ldcsk, tau, wrk, lwrk, info)
            call dorgqr(m, kc, kc, csk, ldcsk, tau, wrk, lwrk, info)
        end do
    end associate

    !-- 2. |..QTA..|..tmpU..|..wrk..|
        i_QTA = 1
        i_tmpU = i_QTA + mem_QTA
        i_wrk = i_tmpU + mem_tmpU
        lwrk = lwork - i_wrk + 1
    associate &
    (QTA => work(i_QTA:), tmpU => work(i_tmpU:), wrk => work(i_wrk:))
        call B2BA('t', kc, n, m, ONE, csk, ldcsk, ZERO, QTA, kc, info)
        call dgesdd_q('s', kc, n, QTA, kc, S, tmpU, kc, VT, ldvt, wrk, lwrk, iwork, liwork, info)
        call dgemm('n', 'n', m, kc, kc, ONE, csk, ldcsk, tmpU, kc, ZERO, U, ldu)
    end associate
    end subroutine dgersvd1

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine sgersvd2 &
    (m, n, kc, csk, ldcsk, test_row, kr, rsk, ldrsk, S, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,      only: &
        WP => SP
    use maria_constants_mod,  only: &
        ZERO => S_ZERO,             &
        ONE => S_ONE
    use maria_la_core_mod,    only: &
        MM => smatmul,              &
        sgeqrf,                     &
        sorgqr,                     &
        sormqr,                     &
        sgesdd_q,                   &
        sgemm,                      &
        strsm,                      &
        sroundup_lwork
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS,                &
        BAD_IF_MORE
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer,       intent(in)                :: m
        integer,       intent(in)                :: n
        integer,       intent(in)                :: kc
        real(WP),      intent(inout), contiguous :: csk(:)
        integer,       intent(in)                :: ldcsk
        procedure(MM), intent(in),    pointer    :: test_row
        integer,       intent(in)                :: kr
        real(WP),      intent(inout), contiguous :: rsk(:)
        integer,       intent(in)                :: ldrsk
        real(WP),      intent(out),   contiguous :: S(:)
        real(WP),      intent(out),   contiguous :: U(:)
        integer,       intent(in)                :: ldu
        real(WP),      intent(out),   contiguous :: VT(:)
        integer,       intent(in)                :: ldvt
        real(WP),      intent(out),   contiguous :: work(:)
        integer,       intent(in)                :: lwork
        integer,       intent(out),   contiguous :: iwork(:)
        integer,       intent(in)                :: liwork
        integer,       intent(out)               :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SGERSVD2'
        logical                 :: lw_query
        integer                 :: req_lw, req_liw, req_lw1, req_lw2, lw_proc, liw_proc, &
                                   mem_B, i_B, mem_tau, i_tau, mem_tmpU, i_tmpU, &
                                   lwrk, i_wrk, ifoo(1)
        real(WP)                :: foo(1) 
    
    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, m, 0)) then
            info = -1
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, kc, 0)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, kc, n)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldcsk, max(1,m))) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, kr, kc)) then
            info = -7
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, kr, m)) then
            info = -7
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldrsk, max(1,kr))) then
            info = -9
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -12
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,kc))) then
            info = -14
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. kc == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) exit wrkwrk

        !-- 1. Storage: tau
        mem_tau = kc
        !-- 1. Query procedures
        call sgeqrf(m, kc, foo, ldcsk, foo, foo, -1, info)
        lw_proc = int(foo(1))
        call sorgqr(m, kc, kc, foo, ldcsk, foo, foo, -1, info)
        lw_proc = max(lw_proc, int(foo(1)))
            
        req_lw1 = lw_proc + mem_tau

        !-- 2. Storage: B, tau, tmpU
        mem_B = kr * kc
!        mem_tau = kc
        mem_tmpU = kc * kc
        !-- 2. Query procedures
        call sgeqrf(kr, kc, foo, kr, foo, foo, -1, info)
        lw_proc = int(foo(1))
        call sormqr('l', 't', kr, n, kc, foo, kr, foo, rsk, ldrsk, foo, -1, info)
        lw_proc = max(lw_proc, int(foo(1)))
        call sgesdd_q('s', kc, n, foo, ldrsk, foo, foo, kc, foo, ldvt, foo, -1, ifoo, -1, info)
        lw_proc = max(lw_proc, int(foo(1)))
        liw_proc = ifoo(1)

        req_lw2 = lw_proc + mem_B + mem_tau + mem_tmpU

        !-- Total
        req_lw = max(req_lw1, req_lw2)
        req_liw = liw_proc

        lw_query = (lwork == -1 .or. liwork == -1)
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -16
            exit wrkwrk
        end if
        if (arg_is_bad(BAD_IF_LESS, liwork, req_liw)) then
            info = -18
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------

    !-- 1. |..tau..|..wrk..|
        i_tau = 1
        i_wrk = i_tau + mem_tau
        lwrk = lwork - i_wrk + 1
    associate &
    (tau => work(i_tau:), wrk => work(i_wrk:))
        call sgeqrf(m, kc, csk, ldcsk, tau, wrk, lwrk, info)
        call sorgqr(m, kc, kc, csk, ldcsk, tau, wrk, lwrk, info)
    end associate

    !-- 2. |..B..|..tau..|..tmpU..|
        i_B = 1
        i_tau = i_B + mem_B
        i_tmpU = i_tau + mem_tau
        i_wrk = i_tmpU + mem_tmpU
        lwrk = lwork - i_wrk + 1
    associate &
    (B => work(i_B:), tau => work(i_tau:), tmpU => work(i_tmpU:), wrk => work(i_wrk:))
        call test_row('n', kr, kc, m, ONE, csk, ldcsk, ZERO, B, kr, info)
        call sgeqrf(kr, kc, B, kr, tau, wrk, lwrk, info)
        call sormqr('l', 't', kr, n, kc, B, kr, tau, rsk, ldrsk, wrk, lwrk, info)
        call strsm('l', 'u', 'n', 'n', kc, n, ONE, B, kr, rsk, ldrsk)
        call sgesdd_q('s', kc, n, rsk, ldrsk, S, tmpU, kc, VT, ldvt, wrk, lwrk, iwork, liwork, info)
        call sgemm('n', 'n', m, kc, kc, ONE, csk, ldcsk, tmpU, kc, ZERO, U, ldu)
    end associate
    end subroutine sgersvd2

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dgersvd2 &
    (m, n, kc, csk, ldcsk, test_row, kr, rsk, ldrsk, S, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,      only: &
        WP => DP
    use maria_constants_mod,  only: &
        ZERO => D_ZERO,             &
        ONE => D_ONE
    use maria_la_core_mod,    only: &
        MM => dmatmul,              &
        dgeqrf,                     &
        dorgqr,                     &
        dormqr,                     &
        dgesdd_q,                   &
        dgemm,                      &
        dtrsm,                      &
        droundup_lwork
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS,                &
        BAD_IF_MORE
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer,       intent(in)                :: m
        integer,       intent(in)                :: n
        integer,       intent(in)                :: kc
        real(WP),      intent(inout), contiguous :: csk(:)
        integer,       intent(in)                :: ldcsk
        procedure(MM), intent(in),    pointer    :: test_row
        integer,       intent(in)                :: kr
        real(WP),      intent(inout), contiguous :: rsk(:)
        integer,       intent(in)                :: ldrsk
        real(WP),      intent(out),   contiguous :: S(:)
        real(WP),      intent(out),   contiguous :: U(:)
        integer,       intent(in)                :: ldu
        real(WP),      intent(out),   contiguous :: VT(:)
        integer,       intent(in)                :: ldvt
        real(WP),      intent(out),   contiguous :: work(:)
        integer,       intent(in)                :: lwork
        integer,       intent(out),   contiguous :: iwork(:)
        integer,       intent(in)                :: liwork
        integer,       intent(out)               :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DGERSVD2'
        logical                 :: lw_query
        integer                 :: req_lw, req_liw, req_lw1, req_lw2, lw_proc, liw_proc, &
                                   mem_B, i_B, mem_tau, i_tau, mem_tmpU, i_tmpU, &
                                   lwrk, i_wrk, ifoo(1)
        real(WP)                :: foo(1) 
    
    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, m, 0)) then
            info = -1
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, kc, 0)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, kc, n)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldcsk, max(1,m))) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, kr, kc)) then
            info = -7
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, kr, m)) then
            info = -7
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldrsk, max(1,kr))) then
            info = -9
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -12
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,kc))) then
            info = -14
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. kc == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) exit wrkwrk

        !-- 1. Storage: tau
        mem_tau = kc
        !-- 1. Query procedures
        call dgeqrf(m, kc, foo, ldcsk, foo, foo, -1, info)
        lw_proc = int(foo(1))
        call dorgqr(m, kc, kc, foo, ldcsk, foo, foo, -1, info)
        lw_proc = max(lw_proc, int(foo(1)))
            
        req_lw1 = lw_proc + mem_tau

        !-- 2. Storage: B, tau, tmpU
        mem_B = kr * kc
!        mem_tau = kc
        mem_tmpU = kc * kc
        !-- 2. Query procedures
        call dgeqrf(kr, kc, foo, kr, foo, foo, -1, info)
        lw_proc = int(foo(1))
        call dormqr('l', 't', kr, n, kc, foo, kr, foo, rsk, ldrsk, foo, -1, info)
        lw_proc = max(lw_proc, int(foo(1)))
        call dgesdd_q('s', kc, n, foo, ldrsk, foo, foo, kc, foo, ldvt, foo, -1, ifoo, -1, info)
        lw_proc = max(lw_proc, int(foo(1)))
        liw_proc = ifoo(1)

        req_lw2 = lw_proc + mem_B + mem_tau + mem_tmpU

        !-- Total
        req_lw = max(req_lw1, req_lw2)
        req_liw = liw_proc

        lw_query = (lwork == -1 .or. liwork == -1)
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -16
            exit wrkwrk
        end if
        if (arg_is_bad(BAD_IF_LESS, liwork, req_liw)) then
            info = -18
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------

    !-- 1. |..tau..|..wrk..|
        i_tau = 1
        i_wrk = i_tau + mem_tau
        lwrk = lwork - i_wrk + 1
    associate &
    (tau => work(i_tau:), wrk => work(i_wrk:))
        call dgeqrf(m, kc, csk, ldcsk, tau, wrk, lwrk, info)
        call dorgqr(m, kc, kc, csk, ldcsk, tau, wrk, lwrk, info)
    end associate

    !-- 2. |..B..|..tau..|..tmpU..|
        i_B = 1
        i_tau = i_B + mem_B
        i_tmpU = i_tau + mem_tau
        i_wrk = i_tmpU + mem_tmpU
        lwrk = lwork - i_wrk + 1
    associate &
    (B => work(i_B:), tau => work(i_tau:), tmpU => work(i_tmpU:), wrk => work(i_wrk:))
        call test_row('n', kr, kc, m, ONE, csk, ldcsk, ZERO, B, kr, info)
        call dgeqrf(kr, kc, B, kr, tau, wrk, lwrk, info)
        call dormqr('l', 't', kr, n, kc, B, kr, tau, rsk, ldrsk, wrk, lwrk, info)
        call dtrsm('l', 'u', 'n', 'n', kc, n, ONE, B, kr, rsk, ldrsk)
        call dgesdd_q('s', kc, n, rsk, ldrsk, S, tmpU, kc, VT, ldvt, wrk, lwrk, iwork, liwork, info)
        call dgemm('n', 'n', m, kc, kc, ONE, csk, ldcsk, tmpU, kc, ZERO, U, ldu)
    end associate
    end subroutine dgersvd2

    !------------------------------------------------------------------------------------------------------------------------

end submodule maria_lr_tsvd_sub
