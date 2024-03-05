!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the implementation of the [[maria_lr_la_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Implements the [[maria_lr_la_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
submodule (maria_lr_la_mod) maria_lr_la_sub
implicit none (type, external)

contains
    !------------------------------------------------------------------------------------------------------------------------

    module function slrval &
    (m, n, i, j, r, U, ldU, VT, ldVT, info)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        ZERO => S_ZERO
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS,               &
        BAD_IF_MORE
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_la_core_mod,   only: &
        sdot
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)              :: m
        integer,  intent(in)              :: n
        integer,  intent(in)              :: i
        integer,  intent(in)              :: j
        integer,  intent(in)              :: r
        real(WP), intent(in),  contiguous :: U(:)
        integer,  intent(in)              :: ldU
        real(WP), intent(in),  contiguous :: VT(:)
        integer,  intent(in)              :: ldVT
        integer,  intent(out)             :: info
        real(WP)                          :: slrval

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SLRVAL'
    
    !-- Default values ---------------------------------------------------------
        slrval = ZERO

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, m, 1)) then
            info = -1
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, n, 1)) then
            info = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, i, 1)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, i, m)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, j, 1)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, j, n)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -7
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,r))) then
            info = -9
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (r == 0) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        slrval = sdot(r, U(i:), ldu, VT(1 + (j-1)*ldvt:), 1)
    end function slrval

    !------------------------------------------------------------------------------------------------------------------------

    module function dlrval &
    (m, n, i, j, r, U, ldU, VT, ldVT, info)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        ZERO => D_ZERO
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS,               &
        BAD_IF_MORE
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_la_core_mod,   only: &
        ddot
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)              :: m
        integer,  intent(in)              :: n
        integer,  intent(in)              :: i
        integer,  intent(in)              :: j
        integer,  intent(in)              :: r
        real(WP), intent(in),  contiguous :: U(:)
        integer,  intent(in)              :: ldU
        real(WP), intent(in),  contiguous :: VT(:)
        integer,  intent(in)              :: ldVT
        integer,  intent(out)             :: info
        real(WP)                          :: dlrval

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DLRVAL'
    
    !-- Default values ---------------------------------------------------------
        dlrval = ZERO

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, m, 1)) then
            info = -1
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, n, 1)) then
            info = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, i, 1)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, i, m)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, j, 1)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, j, n)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -7
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,r))) then
            info = -9
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (r == 0) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        associate(x => U(i:), y => VT(1 + (j-1)*ldvt:))
            dlrval = ddot(r, x, ldu, y, 1)
        end associate
    end function dlrval

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine slr2full &
    (m, n, r, U, ldU, VT, ldVT, A, lda, info)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        ZERO => S_ZERO,            &
        ONE => S_ONE
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_la_core_mod,   only: &
        sgemm,                     &
        slaset
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)              :: m
        integer,  intent(in)              :: n
        integer,  intent(in)              :: r
        real(WP), intent(in),  contiguous :: U(:)
        integer,  intent(in)              :: ldU
        real(WP), intent(in),  contiguous :: VT(:)
        integer,  intent(in)              :: ldVT
        real(WP), intent(out), contiguous :: A(:)
        integer,  intent(in)              :: ldA
        integer,  intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SLR2FULL'
    
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
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,r))) then
            info = -7
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldA, max(1,m))) then
            info = -9
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (r == 0) then
            call slaset('a', m, n, ZERO, ZERO, A, lda)
        else
            call sgemm('n', 'n', m, n, r, ONE, U, ldu, VT, ldvt, ZERO, A, lda)
        end if
    end subroutine slr2full

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dlr2full &
    (m, n, r, U, ldU, VT, ldVT, A, lda, info)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        ZERO => D_ZERO,            &
        ONE => D_ONE
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_la_core_mod,   only: &
        dgemm,                     &
        dlaset
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)              :: m
        integer,  intent(in)              :: n
        integer,  intent(in)              :: r
        real(WP), intent(in),  contiguous :: U(:)
        integer,  intent(in)              :: ldU
        real(WP), intent(in),  contiguous :: VT(:)
        integer,  intent(in)              :: ldVT
        real(WP), intent(out), contiguous :: A(:)
        integer,  intent(in)              :: ldA
        integer,  intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DLR2FULL'
    
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
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,r))) then
            info = -7
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldA, max(1,m))) then
            info = -9
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (r == 0) then
            call dlaset('a', m, n, ZERO, ZERO, A, lda)
        else
            call dgemm('n', 'n', m, n, r, ONE, U, ldu, VT, ldvt, ZERO, A, lda)
        end if
    end subroutine dlr2full

    !------------------------------------------------------------------------------------------------------------------------

    module function slrdotf &
    (m, n, r, U, ldU, VT, ldVT, k, A, ldA, BT, ldBT, work, lwork, info)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        ZERO => S_ZERO,            &
        ONE => S_ONE
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_la_core_mod,   only: &
        sroundup_lwork,            &
        sgedotf,                   &
        sgemm
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)              :: m
        integer,  intent(in)              :: n
        integer,  intent(in)              :: r
        real(WP), intent(in),  contiguous :: U(:)
        integer,  intent(in)              :: ldU
        real(WP), intent(in),  contiguous :: VT(:)
        integer,  intent(in)              :: ldVT
        integer,  intent(in)              :: k
        real(WP), intent(in),  contiguous :: A(:)
        integer,  intent(in)              :: ldA
        real(WP), intent(in),  contiguous :: BT(:)
        integer,  intent(in)              :: ldBT
        real(WP), intent(out), contiguous :: work(:)
        integer,  intent(in)              :: lwork
        integer,  intent(out)             :: info
        real(WP)                          :: slrdotf

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SLRDOTF'
        logical                 :: lw_query
        integer                 :: req_lw, mem_ATU, i_ATU, mem_BTV, i_BTV
    
    !-- Default values ---------------------------------------------------------
        slrdotf = ZERO

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
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,r))) then
            info = -7
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, k, 0)) then
            info = -8
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldA, max(1,m))) then
            info = -10
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldBT, max(1,k))) then
            info = -12
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. r == 0 .or. k == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) then
            mem_ATU = 0
            exit wrkwrk
        end if

    !-- Storage [s]: ATU, BTV --------------------------------------------------
        mem_ATU = k * r
        mem_BTV = k * r

        req_lw  = mem_ATU + mem_BTV

        lw_query = (lwork == -1)
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -14
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        !-- Slice workspace --------------------------------------------------------
        !-- |.....|.....|
        !--   ATU   BTV
        i_ATU = 1
        i_BTV = i_ATU + mem_ATU
    associate(ATU => work(i_ATU:), BTV => work(i_BTV:))
        call sgemm('t', 'n', k, r, m, ONE, A, ldA, U, ldU, ZERO, ATU, k)
        call sgemm('n', 't', k, r, n, ONE, BT, ldBT, VT, ldVT, ZERO, BTV, k)
        slrdotf = sgedotf('n', k, r, ATU, k, BTV, k, info)
    end associate
    end function slrdotf

    !------------------------------------------------------------------------------------------------------------------------

    module function dlrdotf &
    (m, n, r, U, ldU, VT, ldVT, k, A, ldA, BT, ldBT, work, lwork, info)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        ZERO => D_ZERO,            &
        ONE => D_ONE
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_la_core_mod,   only: &
        droundup_lwork,            &
        dgedotf,                   &
        dgemm
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)              :: m
        integer,  intent(in)              :: n
        integer,  intent(in)              :: r
        real(WP), intent(in),  contiguous :: U(:)
        integer,  intent(in)              :: ldU
        real(WP), intent(in),  contiguous :: VT(:)
        integer,  intent(in)              :: ldVT
        integer,  intent(in)              :: k
        real(WP), intent(in),  contiguous :: A(:)
        integer,  intent(in)              :: ldA
        real(WP), intent(in),  contiguous :: BT(:)
        integer,  intent(in)              :: ldBT
        real(WP), intent(out), contiguous :: work(:)
        integer,  intent(in)              :: lwork
        integer,  intent(out)             :: info
        real(WP)                          :: dlrdotf

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DLRDOTF'
        logical                 :: lw_query
        integer                 :: req_lw, mem_ATU, i_ATU, mem_BTV, i_BTV

    !-- Default values ---------------------------------------------------------
        dlrdotf = ZERO

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
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,r))) then
            info = -7
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, k, 0)) then
            info = -8
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldA, max(1,m))) then
            info = -10
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldBT, max(1,k))) then
            info = -12
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. r == 0 .or. k == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) then
            mem_ATU = 0
            exit wrkwrk
        end if

    !-- Storage [s]: ATU, BTV --------------------------------------------------
        mem_ATU = k * r
        mem_BTV = k * r

        req_lw  = mem_ATU + mem_BTV

        lw_query = (lwork == -1)
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -14
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        !-- Slice workspace --------------------------------------------------------
        !-- |.....|.....|
        !--   ATU   BTV
        i_ATU = 1
        i_BTV = i_ATU + mem_ATU
    associate(ATU => work(i_ATU:), BTV => work(i_BTV:))
        call dgemm('t', 'n', k, r, m, ONE, A, ldA, U, ldU, ZERO, ATU, k)
        call dgemm('n', 't', k, r, n, ONE, BT, ldBT, VT, ldVT, ZERO, BTV, k)
        dlrdotf = dgedotf('n', k, r, ATU, k, BTV, k, info)
    end associate
    end function dlrdotf

    !------------------------------------------------------------------------------------------------------------------------

    module function slrnrmf &
    (m, n, r, U, ldU, VT, ldVT, work, lwork, info)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        ZERO => S_ZERO,            &
        ONE => S_ONE
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_la_core_mod,   only: &
        sroundup_lwork
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)              :: m
        integer,  intent(in)              :: n
        integer,  intent(in)              :: r
        real(WP), intent(in),  contiguous :: U(:)
        integer,  intent(in)              :: ldU
        real(WP), intent(in),  contiguous :: VT(:)
        integer,  intent(in)              :: ldVT
        real(WP), intent(out), contiguous :: work(:)
        integer,  intent(in)              :: lwork
        integer,  intent(out)             :: info
        real(WP)                          :: slrnrmf

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SLRNRMF'
        logical                 :: lw_query
        integer                 :: req_lw, lw_slrdotf
        real(WP)                :: dot, foo(1)

    !-- Default values ---------------------------------------------------------
        slrnrmf = ZERO
        work(1) = ONE

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
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,r))) then
            info = -7
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. r == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) then
            exit wrkwrk
        end if

    !-- Query procedures -------------------------------------------------------
        dot = slrdotf(m, n, r, foo, ldU, foo, ldVT, r, foo, ldU, foo, ldVT, foo, -1, info)
        lw_slrdotf = int(foo(1))

        req_lw = lw_slrdotf

        lw_query = (lwork == -1)
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -9
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        slrnrmf = sqrt(slrdotf(m, n, r, U, ldU, VT, ldVT, r, U, ldU, VT, ldVT, work, lwork, info))
    end function slrnrmf

    !------------------------------------------------------------------------------------------------------------------------

    module function dlrnrmf &
    (m, n, r, U, ldU, VT, ldVT, work, lwork, info)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        ZERO => D_ZERO,            &
        ONE => D_ONE
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_la_core_mod,   only: &
        droundup_lwork
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)              :: m
        integer,  intent(in)              :: n
        integer,  intent(in)              :: r
        real(WP), intent(in),  contiguous :: U(:)
        integer,  intent(in)              :: ldU
        real(WP), intent(in),  contiguous :: VT(:)
        integer,  intent(in)              :: ldVT
        real(WP), intent(out), contiguous :: work(:)
        integer,  intent(in)              :: lwork
        integer,  intent(out)             :: info
        real(WP)                          :: dlrnrmf

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DLRNRMF'
        logical                 :: lw_query
        integer                 :: req_lw, lw_dlrdotf
        real(WP)                :: dot, foo(1)

    !-- Default values ---------------------------------------------------------
        dlrnrmf = ZERO
        work(1) = ONE

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
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,r))) then
            info = -7
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. r == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) then
            exit wrkwrk
        end if

    !-- Query procedures -------------------------------------------------------
        dot = dlrdotf(m, n, r, foo, ldU, foo, ldVT, r, foo, ldU, foo, ldVT, foo, -1, info)
        lw_dlrdotf = int(foo(1))

        req_lw = lw_dlrdotf

        lw_query = (lwork == -1)
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -9
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        dlrnrmf = sqrt(dlrdotf(m, n, r, U, ldU, VT, ldVT, r, U, ldU, VT, ldVT, work, lwork, info))
    end function dlrnrmf

    !------------------------------------------------------------------------------------------------------------------------

    module function slrnrmf_diff &
    (m, n, r, U, ldU, VT, ldVT, k, A, ldA, BT, ldBT, work, lwork, info)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        ZERO => S_ZERO,            &
        ONE => S_ONE
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_la_core_mod,   only: &
        sroundup_lwork
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)              :: m
        integer,  intent(in)              :: n
        integer,  intent(in)              :: r
        real(WP), intent(in),  contiguous :: U(:)
        integer,  intent(in)              :: ldU
        real(WP), intent(in),  contiguous :: VT(:)
        integer,  intent(in)              :: ldVT
        integer,  intent(in)              :: k
        real(WP), intent(in),  contiguous :: A(:)
        integer,  intent(in)              :: ldA
        real(WP), intent(in),  contiguous :: BT(:)
        integer,  intent(in)              :: ldBT
        real(WP), intent(out), contiguous :: work(:)
        integer,  intent(in)              :: lwork
        integer,  intent(out)             :: info
        real(WP)                          :: slrnrmf_diff

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SLRNRMF_DIFF'
        logical                 :: lw_query
        integer                 :: req_lw, lw_slrdotf
        real(WP)                :: nrm_UVT2, nrm_ABT2, dot, foo(1)

    !-- Default values ---------------------------------------------------------
        slrnrmf_diff = ZERO
        work(1) = ONE

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
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,r))) then
            info = -7
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, k, 0)) then
            info = -8
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldA, max(1,m))) then
            info = -10
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldBT, max(1,k))) then
            info = -12
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. (r == 0 .and. k == 0)) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) then
            exit wrkwrk
        end if

    !-- Query proceures ----------------------------------------------------
        dot = slrdotf(m, n, r, foo, ldU, foo, ldVT, k, foo, ldA, foo, ldBT, foo, -1, info)
        lw_slrdotf = int(foo(1))
        dot = slrdotf(m, n, r, foo, ldU, foo, ldVT, r, foo, ldU, foo, ldVT, foo, -1, info)
        lw_slrdotf = max(lw_slrdotf, int(foo(1)))
        dot = slrdotf(m, n, k, foo, ldA, foo, ldBT, k, foo, ldA, foo, ldBT, foo, -1, info)
        lw_slrdotf = max(lw_slrdotf, int(foo(1)))

        req_lw = lw_slrdotf

        lw_query = (lwork == -1)
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -14
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        nrm_UVT2 = slrdotf(m, n, r, U, ldU, VT, ldVT, r, U, ldU, VT, ldVT, work, lwork, info)
        nrm_ABT2 = slrdotf(m, n, k, A, ldA, BT, ldBT, k, A, ldA, BT, ldBT, work, lwork, info)
        dot = slrdotf(m, n, r, U, ldU, VT, ldVT, k, A, ldA, BT, ldBT, work, lwork, info)
        slrnrmf_diff = sqrt(nrm_UVT2 + nrm_ABT2 - 2 * dot)
    end function slrnrmf_diff

    !------------------------------------------------------------------------------------------------------------------------

    module function dlrnrmf_diff &
    (m, n, r, U, ldU, VT, ldVT, k, A, ldA, BT, ldBT, work, lwork, info)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        ZERO => D_ZERO,            &
        ONE => D_ONE
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_la_core_mod,   only: &
        droundup_lwork
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)              :: m
        integer,  intent(in)              :: n
        integer,  intent(in)              :: r
        real(WP), intent(in),  contiguous :: U(:)
        integer,  intent(in)              :: ldU
        real(WP), intent(in),  contiguous :: VT(:)
        integer,  intent(in)              :: ldVT
        integer,  intent(in)              :: k
        real(WP), intent(in),  contiguous :: A(:)
        integer,  intent(in)              :: ldA
        real(WP), intent(in),  contiguous :: BT(:)
        integer,  intent(in)              :: ldBT
        real(WP), intent(out), contiguous :: work(:)
        integer,  intent(in)              :: lwork
        integer,  intent(out)             :: info
        real(WP)                          :: dlrnrmf_diff

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DLRNRMF_DIFF'
        logical                 :: lw_query
        integer                 :: req_lw, lw_dlrdotf
        real(WP)                :: nrm_UVT2, nrm_ABT2, dot, foo(1)

    !-- Default values ---------------------------------------------------------
        dlrnrmf_diff = ZERO
        work(1) = ONE

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
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,r))) then
            info = -7
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, k, 0)) then
            info = -8
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldA, max(1,m))) then
            info = -10
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldBT, max(1,k))) then
            info = -12
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. (r == 0 .and. k == 0)) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) then
            exit wrkwrk
        end if

    !-- Query proceures ----------------------------------------------------
        dot = dlrdotf(m, n, r, foo, ldU, foo, ldVT, k, foo, ldA, foo, ldBT, foo, -1, info)
        lw_dlrdotf = int(foo(1))
        dot = dlrdotf(m, n, r, foo, ldU, foo, ldVT, r, foo, ldU, foo, ldVT, foo, -1, info)
        lw_dlrdotf = max(lw_dlrdotf, int(foo(1)))
        dot = dlrdotf(m, n, k, foo, ldA, foo, ldBT, k, foo, ldA, foo, ldBT, foo, -1, info)
        lw_dlrdotf = max(lw_dlrdotf, int(foo(1)))

        req_lw = lw_dlrdotf

        lw_query = (lwork == -1)
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -14
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        nrm_UVT2 = dlrdotf(m, n, r, U, ldU, VT, ldVT, r, U, ldU, VT, ldVT, work, lwork, info)
        nrm_ABT2 = dlrdotf(m, n, k, A, ldA, BT, ldBT, k, A, ldA, BT, ldBT, work, lwork, info)
        dot = dlrdotf(m, n, r, U, ldU, VT, ldVT, k, A, ldA, BT, ldBT, work, lwork, info)
        dlrnrmf_diff = sqrt(nrm_UVT2 + nrm_ABT2 - 2 * dot)
    end function dlrnrmf_diff

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine lrort_rank &
    (side, m, n, r, newr, info)
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        character(1), intent(in)  :: side
        integer,      intent(in)  :: m
        integer,      intent(in)  :: n
        integer,      intent(in)  :: r
        integer,      intent(out) :: newr
        integer,      intent(out) :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'LRORT_RANK'
        logical                 :: leftort

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        select case (side)
            case ('l', 'L')
                leftort = .true.
            case ('r', 'R')
                leftort = .false.
            case default
                leftort = .false.
                info = -1
                exit sanity
        end select
        if (arg_is_bad(BAD_IF_LESS, m, 0)) then
            info = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            info = -4
            exit sanity
        end if
    end block sanity

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (leftort) then
            newr = min(m, r)
        else
            newr = min(n, r)
        end if
    end subroutine lrort_rank

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine slrort &
    (side, m, n, r, U, ldu, VT, ldvt, tau, work, lwork, info)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        ONE => S_ONE
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_la_core_mod,   only: &
        sroundup_lwork,            &
        sgeqrf,                    &
        sgelqf,                    &
        strmm,                     &
        sgemm
    !-- Input/output arguments -------------------------------------------------
        character(1), intent(in)                :: side
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        integer,      intent(in)                :: r
        real(WP),     intent(inout), contiguous :: U(:)
        integer,      intent(in)                :: ldU
        real(WP),     intent(inout), contiguous :: VT(:)
        integer,      intent(in)                :: ldVT
        real(WP),     intent(out),   contiguous :: tau(:)
        real(WP),     intent(out),   contiguous :: work(:)
        integer,      intent(in)                :: lwork
        integer,      intent(out)               :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SLRORT'
        logical                 :: lw_query, leftort
        integer                 :: req_lw, lw_sgeqrf, lw_sgelqf
        real(WP)                :: foo(1) 
    
    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        select case (side)
            case ('l', 'L')
                leftort = .true.
            case ('r', 'R')
                leftort = .false.
            case default
                leftort = .false.
                info = -1
                exit sanity
        end select
        if (arg_is_bad(BAD_IF_LESS, m, 0)) then
            info = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -6
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,r))) then
            info = -8
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. r == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) then
            exit wrkwrk
        end if

        if (leftort) then
        !-- Query procedures ---------------------------------------------------
            call sgeqrf(m, r, foo, ldU, tau, foo, -1, info)
            lw_sgeqrf = int(foo(1))

            req_lw = lw_sgeqrf
        else
        !-- Query procedures ---------------------------------------------------
            call sgelqf(r, n, foo, ldVT, tau, foo, -1, info)
            lw_sgelqf = int(foo(1))

            req_lw = lw_sgelqf
        end if

        lw_query = (lwork == -1)
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -11
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (leftort) then
            call sgeqrf(m, r, U, ldU, tau, work, lwork, info)
            call strmm('l', 'u', 'n', 'n', min(m,r), n, ONE, U, ldU, VT, ldVT)
            if (m < r) then
            associate(subU => U(1 + m*ldU:), subVT => VT(m+1:))
                call sgemm('n', 'n', m, n, r-m, ONE, subU, ldU, subVT, ldVT, ONE, VT, ldVT)
            end associate
            end if
        else
            call sgelqf(r, n, VT, ldVT, tau, work, lwork, info)
            call strmm('r', 'l', 'n', 'n', m, min(n,r), ONE, VT, ldVT, U, ldU)
            if (n < r) then
            associate(subU => U(1 + n*ldU:), subVT => VT(n+1:))
                call sgemm('n', 'n', m, n, r-n, ONE, subU, ldU, subVT, ldVT, ONE, U, ldU)
            end associate
            end if
        end if
    end subroutine slrort

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dlrort &
    (side, m, n, r, U, ldu, VT, ldvt, tau, work, lwork, info)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        ONE => D_ONE
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_la_core_mod,   only: &
        droundup_lwork,            &
        dgeqrf,                    &
        dgelqf,                    &
        dtrmm,                     &
        dgemm
    !-- Input/output arguments -------------------------------------------------
        character(1), intent(in)                :: side
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        integer,      intent(in)                :: r
        real(WP),     intent(inout), contiguous :: U(:)
        integer,      intent(in)                :: ldU
        real(WP),     intent(inout), contiguous :: VT(:)
        integer,      intent(in)                :: ldVT
        real(WP),     intent(out),   contiguous :: tau(:)
        real(WP),     intent(out),   contiguous :: work(:)
        integer,      intent(in)                :: lwork
        integer,      intent(out)               :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DLRORT'
        logical                 :: lw_query, leftort
        integer                 :: req_lw, lw_dgeqrf, lw_dgelqf
        real(WP)                :: foo(1) 
    
    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        select case (side)
            case ('l', 'L')
                leftort = .true.
            case ('r', 'R')
                leftort = .false.
            case default
                leftort = .false.
                info = -1
                exit sanity
        end select
        if (arg_is_bad(BAD_IF_LESS, m, 0)) then
            info = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -6
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,r))) then
            info = -8
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. r == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) then
            exit wrkwrk
        end if

        if (leftort) then
        !-- Query procedures ---------------------------------------------------
            call dgeqrf(m, r, foo, ldU, tau, foo, -1, info)
            lw_dgeqrf = int(foo(1))

            req_lw = lw_dgeqrf
        else
        !-- Query procedures ---------------------------------------------------
            call dgelqf(r, n, foo, ldVT, tau, foo, -1, info)
            lw_dgelqf = int(foo(1))

            req_lw = lw_dgelqf
        end if

        lw_query = (lwork == -1)
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -11
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (leftort) then
            call dgeqrf(m, r, U, ldU, tau, work, lwork, info)
            call dtrmm('l', 'u', 'n', 'n', min(m,r), n, ONE, U, ldU, VT, ldVT)
            if (m < r) then
            associate(subU => U(1 + m*ldU:), subVT => VT(m+1:))
                call dgemm('n', 'n', m, n, r-m, ONE, subU, ldU, subVT, ldVT, ONE, VT, ldVT)
            end associate
            end if
        else
            call dgelqf(r, n, VT, ldVT, tau, work, lwork, info)
            call dtrmm('r', 'l', 'n', 'n', m, min(n,r), ONE, VT, ldVT, U, ldU)
            if (n < r) then
            associate(subU => U(1 + n*ldU:), subVT => VT(n+1:))
                call dgemm('n', 'n', m, n, r-n, ONE, subU, ldU, subVT, ldVT, ONE, U, ldU)
            end associate
            end if
        end if
    end subroutine dlrort

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine slrsvd_ort &
    (job, side, m, n, r, U, ldu, VT, ldvt, tau, S, Q, ldQ, PT, ldPT, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        ZERO => S_ZERO
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg,            &
        report_runtime_err,        &
        SVD_FAILED_ERR_CODE
    use maria_la_core_mod,   only: &
        sroundup_lwork,            &
        sgesdd_q,                  &
        sormqr,                    &
        sormlq,                    &
        slaset
    !-- Input/output arguments -------------------------------------------------
        character(1), intent(in)                :: job
        character(1), intent(in)                :: side
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        integer,      intent(in)                :: r
        real(WP),     intent(inout), contiguous :: U(:)
        integer,      intent(in)                :: ldU
        real(WP),     intent(inout), contiguous :: VT(:)
        integer,      intent(in)                :: ldVT
        real(WP),     intent(in),    contiguous :: tau(:)
        real(WP),     intent(out),   contiguous :: S(:)
        real(WP),     intent(out),   contiguous :: Q(:)
        integer,      intent(in)                :: ldQ
        real(WP),     intent(out),   contiguous :: PT(:)
        integer,      intent(in)                :: ldPT
        real(WP),     intent(out),   contiguous :: work(:)
        integer,      intent(in)                :: lwork
        integer,      intent(out),   contiguous :: iwork(:)
        integer,      intent(in)                :: liwork
        integer,      intent(out)               :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SLRSVD_ORT'
        logical                 :: lw_query, leftort, compute_uv
        integer                 :: req_lw, req_liw, lw_sgesdd, liw_sgesdd, &
                                   lw_sormqr, lw_sormlq, ifoo(1)
        real(WP)                :: foo(1) 
    
    !-- Default values ---------------------------------------------------------
        leftort = .true.

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        select case(job)
            case ('s', 'S')
                compute_uv = .true.
            case ('n', 'N')
                compute_uv = .false.
            case default
                compute_uv = .false.
                info = -1
                exit sanity
        end select
        select case (side)
            case ('l', 'L')
                leftort = .true.
            case ('r', 'R')
                leftort = .false.
            case default
                leftort = .false.
                info = -2
                exit sanity
        end select
        if (arg_is_bad(BAD_IF_LESS, m, 0)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            info = -5
            exit sanity
        end if
        if (leftort) then
            if (arg_is_bad(BAD_IF_LESS, m, r)) then
                info = -5
                exit sanity
            end if
        else
            if (arg_is_bad(BAD_IF_LESS, n, r)) then
                info = -5
                exit sanity
            end if
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -7
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,r))) then
            info = -9
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldQ, max(1,m))) then
            info = -13
            exit sanity
        end if
        if (leftort) then
            if (arg_is_bad(BAD_IF_LESS, ldPT, max(1,min(r,n)))) then
                info = -15
                exit sanity
            end if
        else
            if (arg_is_bad(BAD_IF_LESS, ldPT, max(1,min(m,r)))) then
                info = -15
                exit sanity
            end if
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. r == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) then
            exit wrkwrk
        end if

        if (leftort) then
        !-- Query procedures ---------------------------------------------------
            call sgesdd_q(job, r, n, foo, ldvt, foo, foo, ldq, foo, ldpt, foo, -1, ifoo, -1, info)
            lw_sgesdd = int(foo(1))
            liw_sgesdd = ifoo(1)
            if (compute_uv) then
                call sormqr(side, 'n', m, min(r,n), r, foo, ldu, foo, foo, ldq, foo, -1, info)
                lw_sormqr = int(foo(1))
            else
                lw_sormqr = 0
            end if

            req_lw = max(lw_sgesdd, lw_sormqr)
            req_liw = liw_sgesdd
        else
        !-- Query procedures ---------------------------------------------------
            call sgesdd_q(job, m, r, foo, ldu, foo, foo, ldq, foo, ldpt, foo, -1, ifoo, -1, info)
            lw_sgesdd = int(foo(1))
            liw_sgesdd = ifoo(1)
            if (compute_uv) then
                call sormlq(side, 'n', min(m,r), n, r, foo, ldvt, foo, foo, ldpt, foo, -1, info)
                lw_sormlq = int(foo(1))
            else
                lw_sormlq = 0
            end if

            req_lw = max(lw_sgesdd, lw_sormlq)
            req_liw = liw_sgesdd
        end if

        lw_query = (lwork == -1 .or. liwork == -1)
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -17
            exit wrkwrk
        end if
        if (arg_is_bad(BAD_IF_LESS, liwork, req_liw)) then
            info = -19
            exit wrkwrk
        end if

    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (leftort) then
            call slaset('a', m, min(r,n), ZERO, ZERO, Q, ldq)
            call sgesdd_q(job, r, n, VT, ldvt, S, Q, ldq, PT, ldpt, work, lwork, iwork, liwork, info)
            if (info > 0) then
                call report_runtime_err(SRNAME, SVD_FAILED_ERR_CODE)
                return
            end if
            if (compute_uv) then
                call sormqr(side, 'n', m, min(r,n), r, U, ldu, tau, Q, ldq, work, lwork, info)
            end if
        else
            call slaset('a', min(m,r), n, ZERO, ZERO, PT, ldpt)
            call sgesdd_q(job, m, r, U, ldu, S, Q, ldq, PT, ldpt, work, lwork, iwork, liwork, info)
            if (info > 0) then
                call report_runtime_err(SRNAME, SVD_FAILED_ERR_CODE)
                return
            end if
            if (compute_uv) then
                call sormlq(side, 'n', min(m,r), n, r, VT, ldvt, tau, PT, ldpt, work, lwork, info)
            end if
        end if
    end subroutine slrsvd_ort

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dlrsvd_ort &
    (job, side, m, n, r, U, ldu, VT, ldvt, tau, S, Q, ldQ, PT, ldPT, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        ZERO => D_ZERO
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg,            &
        report_runtime_err,        &
        SVD_FAILED_ERR_CODE
    use maria_la_core_mod,   only: &
        droundup_lwork,            &
        dgesdd_q,                  &
        dormqr,                    &
        dormlq,                    &
        dlaset
    !-- Input/output arguments -------------------------------------------------
        character(1), intent(in)                :: job
        character(1), intent(in)                :: side
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        integer,      intent(in)                :: r
        real(WP),     intent(inout), contiguous :: U(:)
        integer,      intent(in)                :: ldU
        real(WP),     intent(inout), contiguous :: VT(:)
        integer,      intent(in)                :: ldVT
        real(WP),     intent(in),    contiguous :: tau(:)
        real(WP),     intent(out),   contiguous :: S(:)
        real(WP),     intent(out),   contiguous :: Q(:)
        integer,      intent(in)                :: ldQ
        real(WP),     intent(out),   contiguous :: PT(:)
        integer,      intent(in)                :: ldPT
        real(WP),     intent(out),   contiguous :: work(:)
        integer,      intent(in)                :: lwork
        integer,      intent(out),   contiguous :: iwork(:)
        integer,      intent(in)                :: liwork
        integer,      intent(out)               :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DLRSVD_ORT'
        logical                 :: lw_query, leftort, compute_uv
        integer                 :: req_lw, req_liw, lw_dgesdd, liw_dgesdd, &
                                   lw_dormqr, lw_dormlq, ifoo(1)
        real(WP)                :: foo(1) 
    
    !-- Default values ---------------------------------------------------------
        leftort = .true.

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        select case(job)
            case ('s', 'S')
                compute_uv = .true.
            case ('n', 'N')
                compute_uv = .false.
            case default
                compute_uv = .false.
                info = -1
                exit sanity
        end select
        select case (side)
            case ('l', 'L')
                leftort = .true.
            case ('r', 'R')
                leftort = .false.
            case default
                leftort = .false.
                info = -2
                exit sanity
        end select
        if (arg_is_bad(BAD_IF_LESS, m, 0)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            info = -5
            exit sanity
        end if
        if (leftort) then
            if (arg_is_bad(BAD_IF_LESS, m, r)) then
                info = -5
                exit sanity
            end if
        else
            if (arg_is_bad(BAD_IF_LESS, n, r)) then
                info = -5
                exit sanity
            end if
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -7
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,r))) then
            info = -9
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldQ, max(1,m))) then
            info = -13
            exit sanity
        end if
        if (leftort) then
            if (arg_is_bad(BAD_IF_LESS, ldPT, max(1,min(r,n)))) then
                info = -15
                exit sanity
            end if
        else
            if (arg_is_bad(BAD_IF_LESS, ldPT, max(1,min(m,r)))) then
                info = -15
                exit sanity
            end if
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. r == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) then
            exit wrkwrk
        end if

        if (leftort) then
        !-- Query procedures ---------------------------------------------------
            call dgesdd_q(job, r, n, foo, ldvt, foo, foo, ldq, foo, ldpt, foo, -1, ifoo, -1, info)
            lw_dgesdd = int(foo(1))
            liw_dgesdd = ifoo(1)
            if (compute_uv) then
                call dormqr(side, 'n', m, min(r,n), r, foo, ldu, foo, foo, ldq, foo, -1, info)
                lw_dormqr = int(foo(1))
            else
                lw_dormqr = 0
            end if

            req_lw = max(lw_dgesdd, lw_dormqr)
            req_liw = liw_dgesdd
        else
        !-- Query procedures ---------------------------------------------------
            call dgesdd_q(job, m, r, foo, ldu, foo, foo, ldq, foo, ldpt, foo, -1, ifoo, -1, info)
            lw_dgesdd = int(foo(1))
            liw_dgesdd = ifoo(1)
            if (compute_uv) then
                call dormlq(side, 'n', min(m,r), n, r, foo, ldvt, foo, foo, ldpt, foo, -1, info)
                lw_dormlq = int(foo(1))
            else
                lw_dormlq = 0
            end if

            req_lw = max(lw_dgesdd, lw_dormlq)
            req_liw = liw_dgesdd
        end if

        lw_query = (lwork == -1 .or. liwork == -1)
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -17
            exit wrkwrk
        end if
        if (arg_is_bad(BAD_IF_LESS, liwork, req_liw)) then
            info = -19
            exit wrkwrk
        end if

    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (leftort) then
            call dlaset('a', m, min(r,n), ZERO, ZERO, Q, ldq)
            call dgesdd_q(job, r, n, VT, ldvt, S, Q, ldq, PT, ldpt, work, lwork, iwork, liwork, info)
            if (info > 0) then
                call report_runtime_err(SRNAME, SVD_FAILED_ERR_CODE)
                return
            end if
            if (compute_uv) then
                call dormqr(side, 'n', m, min(r,n), r, U, ldu, tau, Q, ldq, work, lwork, info)
            end if
        else
            call dlaset('a', min(m,r), n, ZERO, ZERO, PT, ldpt)
            call dgesdd_q(job, m, r, U, ldu, S, Q, ldq, PT, ldpt, work, lwork, iwork, liwork, info)
            if (info > 0) then
                call report_runtime_err(SRNAME, SVD_FAILED_ERR_CODE)
                return
            end if
            if (compute_uv) then
                call dormlq(side, 'n', min(m,r), n, r, VT, ldvt, tau, PT, ldpt, work, lwork, info)
            end if
        end if
    end subroutine dlrsvd_ort

    !------------------------------------------------------------------------------------------------------------------------

end submodule maria_lr_la_sub
