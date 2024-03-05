!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the implementation of the [[maria_la_core_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Implements the [[maria_la_core_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
submodule (maria_la_core_mod) maria_la_core_sub
implicit none (type, external)

contains
    !------------------------------------------------------------------------------------------------------------------------

    module subroutine sgesdd_q &
    (jobZ, m, n, A, ldA, S, U, ldU, VT, ldVT, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,    only: &
        WP => SP
    use maria_argcheck_mod, only: &
        arg_is_bad,               &
        BAD_IF_LESS
    use maria_reports_mod,  only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        character(1), intent(in)                :: jobZ
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        real(WP),     intent(inout), contiguous :: A(:)
        integer,      intent(in)                :: ldA
        real(WP),     intent(out),   contiguous :: S(:)
        real(WP),     intent(out),   contiguous :: U(:)
        integer,      intent(in)                :: ldU
        real(WP),     intent(out),   contiguous :: VT(:)
        integer,      intent(in)                :: ldVT
        real(WP),     intent(out),   contiguous :: work(:)
        integer,      intent(in)                :: lwork
        integer,      intent(out),   contiguous :: iwork(:)
        integer,      intent(in)                :: liwork
        integer,      intent(out)               :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SGESDD_Q'
        logical                 :: lw_query   
        integer                 :: req_lw, req_liw

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        call sgesdd(jobZ, m, n, A, ldA, S, U, ldU, VT, ldVT, work, -1, iwork, info)
    end block sanity

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) exit wrkwrk

        req_lw  = int(work(1)) 
        req_liw = 8 * min(m,n)

        lw_query = (lwork == -1 .or. liwork == -1)
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -12
            exit wrkwrk
        end if
        if (arg_is_bad(BAD_IF_LESS, liwork, req_liw)) then
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
        call sgesdd(jobZ, m, n, A, ldA, S, U, ldU, VT, ldVT, work, lwork, iwork, info) 
    end subroutine sgesdd_q

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dgesdd_q &
    (jobZ, m, n, A, ldA, S, U, ldU, VT, ldVT, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,      only: &
        WP => DP
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        character(1), intent(in)                :: jobZ
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        real(WP),     intent(inout), contiguous :: A(:)
        integer,      intent(in)                :: ldA
        real(WP),     intent(out),   contiguous :: S(:)
        real(WP),     intent(out),   contiguous :: U(:)
        integer,      intent(in)                :: ldU
        real(WP),     intent(out),   contiguous :: VT(:)
        integer,      intent(in)                :: ldVT
        real(WP),     intent(out),   contiguous :: work(:)
        integer,      intent(in)                :: lwork
        integer,      intent(out),   contiguous :: iwork(:)
        integer,      intent(in)                :: liwork
        integer,      intent(out)               :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DGESDD_Q'
        logical                 :: lw_query   
        integer                 :: req_lw, req_liw

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        call dgesdd(jobZ, m, n, A, ldA, S, U, ldU, VT, ldVT, work, -1, iwork, info)
    end block sanity

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) exit wrkwrk

        req_lw  = int(work(1)) 
        req_liw = 8 * min(m,n)

        lw_query = (lwork == -1 .or. liwork == -1)
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -12
            exit wrkwrk
        end if
        if (arg_is_bad(BAD_IF_LESS, liwork, req_liw)) then
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
        call dgesdd(jobZ, m, n, A, ldA, S, U, ldU, VT, ldVT, work, lwork, iwork, info) 
    end subroutine dgesdd_q

    !------------------------------------------------------------------------------------------------------------------------

    module function sroundup_lwork &
    (lwork)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        EPS => S_MACHTOL
    !-- Input/output arguments -------------------------------------------------
        integer, intent(in) :: lwork
        real(WP)            :: sroundup_lwork         

    !-- Executable section -----------------------------------------------------
        sroundup_lwork = real(lwork, WP)
        if (nint(sroundup_lwork) < lwork) then
            sroundup_lwork = sroundup_lwork * (1 + EPS)
        end if
    end function sroundup_lwork

    !------------------------------------------------------------------------------------------------------------------------

    module function droundup_lwork &
    (lwork)
    use maria_kinds_mod,      only: &
        WP => DP
    use maria_constants_mod,  only: &
        EPS => D_MACHTOL
    !-- Input/output arguments -------------------------------------------------
        integer, intent(in) :: lwork
        real(WP)            :: droundup_lwork         

    !-- Executable section -----------------------------------------------------       
        droundup_lwork = real(lwork, WP)
        if (nint(droundup_lwork) < lwork) then
            droundup_lwork = droundup_lwork * (1 + EPS)
        end if
    end function droundup_lwork

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine sgescal &
    (m, n, alpha, A, ldA, info)
    use maria_kinds_mod,      only: &
        WP => SP
    use maria_constants_mod,  only: &
        ONE => S_ONE
    use maria_comparison_mod, only: &
        safe_eq
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)                :: m
        integer,  intent(in)                :: n
        real(WP), intent(in)                :: alpha
        real(WP), intent(inout), contiguous :: A(:)
        integer,  intent(in)                :: ldA
        integer,  intent(out)               :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SGESCAL'
        integer                 :: j

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
        if (arg_is_bad(BAD_IF_LESS, ldA, max(1,m))) then
            info = -5
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0) return
        if (safe_eq(alpha, ONE)) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (lda == m) then
            call sscal(m*n, alpha, A, 1)
        else
            do j = 1, n
                call sscal(m, alpha, A(1 + (j-1)*ldA :), 1)
            end do
        end if
    end subroutine sgescal

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dgescal &
    (m, n, alpha, A, ldA, info)
    use maria_kinds_mod,      only: &
        WP => DP
    use maria_constants_mod,  only: &
        ONE => D_ONE
    use maria_comparison_mod, only: &
        safe_eq
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)                :: m
        integer,  intent(in)                :: n
        real(WP), intent(in)                :: alpha
        real(WP), intent(inout), contiguous :: A(:)
        integer,  intent(in)                :: ldA
        integer,  intent(out)               :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DGESCAL'
        integer                 :: j

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
        if (arg_is_bad(BAD_IF_LESS, ldA, max(1,m))) then
            info = -5
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0) return
        if (safe_eq(alpha, ONE)) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (lda == m) then
            call dscal(m*n, alpha, A, 1)
        else
            do j = 1, n
                call dscal(m, alpha, A(1 + (j-1)*ldA :), 1)
            end do
        end if
    end subroutine dgescal

    !------------------------------------------------------------------------------------------------------------------------

    module function sgedotf &
    (transA, m, n, A, ldA, B, ldB, info)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        ZERO => S_ZERO
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        character(1), intent(in)             :: transA
        integer,      intent(in)             :: m
        integer,      intent(in)             :: n
        real(WP),     intent(in), contiguous :: A(:)
        integer,      intent(in)             :: ldA
        real(WP),     intent(in), contiguous :: B(:)
        integer,      intent(in)             :: ldB
        integer,      intent(out)            :: info
        real(WP)                             :: sgedotf

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SGEDOTF'
        logical                 :: transposedA
        integer                 :: j

    !-- Default values ---------------------------------------------------------
    sgedotf = ZERO

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        select case (transA)
            case ('n', 'N')
                transposedA = .false.
            case ('t', 'T')
                transposedA = .true.
            case default
                transposedA = .false.
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
        if (transposedA) then
            if (arg_is_bad(BAD_IF_LESS, ldA, max(1,n))) then
                info = -5
                exit sanity
            end if
        else
            if (arg_is_bad(BAD_IF_LESS, ldA, max(1,m))) then
                info = -5
                exit sanity
            end if
        end if
        if (arg_is_bad(BAD_IF_LESS, ldB, max(1,m))) then
            info = -7
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
        if (transposedA) then
            do j = 1, n
                sgedotf = sgedotf + sdot(m, A(j:), ldA, B(1 + (j-1)*ldB:), 1)
            end do
        else
            if (ldA == m .and. ldA == m) then
                sgedotf = sdot(m*n, A, 1, B, 1)
            else
                do j = 1, n
                    sgedotf = sgedotf + sdot(m, A(1 + (j-1)*ldA:), 1, B(1 + (j-1)*ldB:), 1)
                end do
            end if
        end if
    end function sgedotf

    !------------------------------------------------------------------------------------------------------------------------

    module function dgedotf &
    (transA, m, n, A, ldA, B, ldB, info)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        ZERO => D_ZERO
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        character(1), intent(in)             :: transA
        integer,      intent(in)             :: m
        integer,      intent(in)             :: n
        real(WP),     intent(in), contiguous :: A(:)
        integer,      intent(in)             :: ldA
        real(WP),     intent(in), contiguous :: B(:)
        integer,      intent(in)             :: ldB
        integer,      intent(out)            :: info
        real(WP)                             :: dgedotf

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DGEDOTF'
        logical                 :: transposedA
        integer                 :: j

    !-- Default values ---------------------------------------------------------
    dgedotf = ZERO

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        select case (transA)
            case ('n', 'N')
                transposedA = .false.
            case ('t', 'T')
                transposedA = .true.
            case default
                transposedA = .false.
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
        if (transposedA) then
            if (arg_is_bad(BAD_IF_LESS, ldA, max(1,n))) then
                info = -5
                exit sanity
            end if
        else
            if (arg_is_bad(BAD_IF_LESS, ldA, max(1,m))) then
                info = -5
                exit sanity
            end if
        end if
        if (arg_is_bad(BAD_IF_LESS, ldB, max(1,m))) then
            info = -7
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
        if (transposedA) then
            do j = 1, n
                dgedotf = dgedotf + ddot(m, A(j:), ldA, B(1 + (j-1)*ldB:), 1)
            end do
        else
            if (ldA == m .and. ldA == m) then
                dgedotf = ddot(m*n, A, 1, B, 1)
            else
                do j = 1, n
                    dgedotf = dgedotf + ddot(m, A(1 + (j-1)*ldA:), 1, B(1 + (j-1)*ldB:), 1)
                end do
            end if
        end if
    end function dgedotf

    !------------------------------------------------------------------------------------------------------------------------

    module function sgenrmf &
    (m, n, A, ldA, info)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        ZERO => S_ZERO
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer,      intent(in)             :: m
        integer,      intent(in)             :: n
        real(WP),     intent(in), contiguous :: A(:)
        integer,      intent(in)             :: ldA
        integer,      intent(out)            :: info
        real(WP)                             :: sgenrmf

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SGENRMF'
        integer                 :: j

    !-- Default values ---------------------------------------------------------
    sgenrmf = ZERO

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
        if (arg_is_bad(BAD_IF_LESS, ldA, max(1,m))) then
            info = -4
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
        if (ldA == m .and. ldA == m) then
            sgenrmf = snrm2(m*n, A, 1)
        else
            do j = 1, n
                sgenrmf = sgenrmf + snrm2(m, A(1 + (j-1)*ldA:), 1)**2
            end do
            sgenrmf = sqrt(sgenrmf)
        end if
    end function sgenrmf

    !------------------------------------------------------------------------------------------------------------------------

    module function dgenrmf &
    (m, n, A, ldA, info)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        ZERO => D_ZERO
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer,      intent(in)             :: m
        integer,      intent(in)             :: n
        real(WP),     intent(in), contiguous :: A(:)
        integer,      intent(in)             :: ldA
        integer,      intent(out)            :: info
        real(WP)                             :: dgenrmf

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DGENRMF'
        integer                 :: j

    !-- Default values ---------------------------------------------------------
    dgenrmf = ZERO

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
        if (arg_is_bad(BAD_IF_LESS, ldA, max(1,m))) then
            info = -4
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
        if (ldA == m .and. ldA == m) then
            dgenrmf = dnrm2(m*n, A, 1)
        else
            do j = 1, n
                dgenrmf = dgenrmf + dnrm2(m, A(1 + (j-1)*ldA:), 1)**2
            end do
            dgenrmf = sqrt(dgenrmf)
        end if
    end function dgenrmf

    !------------------------------------------------------------------------------------------------------------------------

    module function sgenrmc &
    (m, n, A, ldA, info, pos)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        ZERO => S_ZERO
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer,      intent(in)             :: m
        integer,      intent(in)             :: n
        real(WP),     intent(in), contiguous :: A(:)
        integer,      intent(in)             :: ldA
        integer,      intent(out)            :: info
        integer,      intent(out), optional  :: pos(2)
        real(WP)                             :: sgenrmc

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SGENRMC'
        logical                 :: take
        integer                 :: i, j
        real(WP)                :: val

    !-- Default values ---------------------------------------------------------
        sgenrmc = ZERO
        if (present(pos)) then
            pos(1) = 1
            pos(2) = 1
        end if
 
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
        if (arg_is_bad(BAD_IF_LESS, ldA, max(1,m))) then
            info = -4
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
        do j = 1, n
            do i = 1, m
                val = abs(A(i + (j-1)*ldA))
                take = isnan(val)
                if (.not. take) take = (val > sgenrmc)
                if (take) then
                    sgenrmc = val
                    if (present(pos)) then
                        pos(1) = i
                        pos(2) = j
                    end if
                    if (isnan(val)) return
                end if
            end do
        end do
    end function sgenrmc

    !------------------------------------------------------------------------------------------------------------------------

    module function dgenrmc &
    (m, n, A, ldA, info, pos)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        ZERO => D_ZERO
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer,      intent(in)             :: m
        integer,      intent(in)             :: n
        real(WP),     intent(in), contiguous :: A(:)
        integer,      intent(in)             :: ldA
        integer,      intent(out)            :: info
        integer,      intent(out), optional  :: pos(2)
        real(WP)                             :: dgenrmc

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DGENRMC'
        logical                 :: take
        integer                 :: i, j
        real(WP)                :: val

    !-- Default values ---------------------------------------------------------
        dgenrmc = ZERO
        if (present(pos)) then
            pos(1) = 1
            pos(2) = 1
        end if

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
        if (arg_is_bad(BAD_IF_LESS, ldA, max(1,m))) then
            info = -4
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
        do j = 1, n
            do i = 1, m
                val = abs(A(i + (j-1)*ldA))
                take = isnan(val)
                if (.not. take) take = (val > dgenrmc)
                if (take) then
                    dgenrmc = val
                    if (present(pos)) then
                        pos(1) = i
                        pos(2) = j
                    end if
                    if (isnan(val)) return
                end if
            end do
        end do
    end function dgenrmc

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine sdgmm &
    (side, m, n, A, lda, D, incd, info)
    use maria_kinds_mod,    only: &
        WP => SP
    use maria_argcheck_mod, only: &
        arg_is_bad,               &
        BAD_IF_LESS
    use maria_reports_mod,  only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        character(1), intent(in)                :: side
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        real(WP),     intent(inout), contiguous :: A(:)
        integer,      intent(in)                :: ldA
        real(WP),     intent(in),    contiguous :: D(:)
        integer,      intent(in)                :: incD
        integer,      intent(out)               :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SDGMM'
        logical                 :: left
        integer                 :: i, posd

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        select case (side)
            case ('r', 'R')
                left= .false.
            case ('l', 'L')
                left = .true.
            case default
                left = .false.
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
        if (arg_is_bad(BAD_IF_LESS, ldA, max(1,m))) then
            info = -5
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
        if (incD == 0) then
            call sgescal(m, n, D(1), A, ldA, info)
            return
        end if

        if (left) then
            if (incd > 0) then
                posd = 1
            else
                posd = 1 - (m-1) * incd
            end if
            do i = 1, m
            associate(row => A(i:))
                call sscal(n, D(posd), row, lda)
                posd = posd + incd
            end associate
            end do
        else
            if (incd > 0) then
                posd = 1
            else
                posd = 1 - (n-1) * incd
            end if
            do i = 1, n
            associate(col => A(1 + (i-1)*lda:))
                call sscal(m, D(posd), col, 1)
                posd = posd + incd
            end associate
            end do
        end if
    end subroutine sdgmm

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine ddgmm &
    (side, m, n, A, lda, D, incd, info)
    use maria_kinds_mod,    only: &
        WP => DP
    use maria_argcheck_mod, only: &
        arg_is_bad,               &
        BAD_IF_LESS
    use maria_reports_mod,  only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        character(1), intent(in)                :: side
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        real(WP),     intent(inout), contiguous :: A(:)
        integer,      intent(in)                :: ldA
        real(WP),     intent(in),    contiguous :: D(:)
        integer,      intent(in)                :: incD
        integer,      intent(out)               :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DDGMM'
        logical                 :: left
        integer                 :: i, posd

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        select case (side)
            case ('r', 'R')
                left= .false.
            case ('l', 'L')
                left = .true.
            case default
                left = .false.
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
        if (arg_is_bad(BAD_IF_LESS, ldA, max(1,m))) then
            info = -5
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
        if (incD == 0) then
            call dgescal(m, n, D(1), A, ldA, info)
            return
        end if

        if (left) then
            if (incd > 0) then
                posd = 1
            else
                posd = 1 - (m-1) * incd
            end if
            do i = 1, m
            associate( x => A(i:) )
                call dscal(n, D(posd), x, lda)
            end associate
                posd = posd + incd
            end do
        else
            if (incd > 0) then
                posd = 1
            else
                posd = 1 - (n-1) * incd
            end if
            do i = 1, n
            associate( x => A(1 + (i-1)*lda:) )
                call dscal(m, D(posd), x, 1)
            end associate
                posd = posd + incd
            end do
        end if
    end subroutine ddgmm

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine sgepiv &
    (what, dir, m, n, A, lda, from, to, ipiv, info, order)
    use maria_kinds_mod,    only: &
        WP => SP
    use maria_utils_mod,    only: &
        swap_pair
    use maria_argcheck_mod, only: &
        arg_is_bad,               &
        BAD_IF_LESS,              &
        BAD_IF_MORE
    use maria_reports_mod,  only: &
        report_bad_arg,           &
        report_runtime_err,       &
        IND_OUT_OF_RANGE_ERR_CODE
    !-- Input/output arguments -------------------------------------------------
        character(1), intent(in)                          :: what
        character(1), intent(in)                          :: dir
        integer,      intent(in)                          :: m
        integer,      intent(in)                          :: n
        real(WP),     intent(inout), contiguous           :: A(:)
        integer,      intent(in)                          :: ldA
        integer,      intent(in)                          :: from
        integer,      intent(in)                          :: to
        integer,      intent(in),    contiguous           :: ipiv(:)
        integer,      intent(out)                         :: info
        integer,      intent(inout), contiguous, optional :: order(:)

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SGEPIV'
        logical                 :: rows, forward
        integer                 :: i, inc, begg, endd, ind

    !-- Sanity check -----------------------------------------------------------
        rows = .false.
        forward = .false.
    sanity: block
        info = 0
        select case (what)
            case ('c', 'C')
                rows = .false.
            case ('r', 'R')
                rows = .true.
            case default
                info = -1
                exit sanity
        end select
        select case (dir)
            case ('b', 'B')
                forward = .false.
            case ('f', 'F')
                forward = .true.
            case default
                info = -2
                exit sanity
        end select
        if (arg_is_bad(BAD_IF_LESS, m, 1)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, n, 1)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldA, max(1,m))) then
            info = -6
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, from, 1)) then
            info = -7
            exit sanity
        end if
        if (rows) then
            if (arg_is_bad(BAD_IF_MORE, from, m)) then
                info = -7
                exit sanity
            end if
        else
            if (arg_is_bad(BAD_IF_MORE, from, n)) then
                info = -7
                exit sanity
            end if
        end if
        if (arg_is_bad(BAD_IF_LESS, to, from)) then
            info = -8
            exit sanity
        end if
        if (rows) then
            if (arg_is_bad(BAD_IF_MORE, to, m)) then
                info = -8
                exit sanity
            end if
        else
            if (arg_is_bad(BAD_IF_MORE, to, n)) then
                info = -8
                exit sanity
            end if
        end if
    end block sanity

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (forward) then
            inc = 1
            begg = from
            endd = to
        else
            inc = -1
            begg = to
            endd = from
        end if

        if (rows) then
            do i = begg, endd, inc
                ind = ipiv(i - from + 1)
                if (ind < 1 .or. ind > m) then
                    info = i - from + 1
                    call report_runtime_err(SRNAME, IND_OUT_OF_RANGE_ERR_CODE)
                    return
                end if
                if (ind /= i) then
                    call sswap(n, A(i:), ldA, A(ind:), ldA)
                    if (present(order)) then
                        call swap_pair(order(i), order(ind))
                    end if
                end if
            end do
        else
            do i = begg, endd, inc
                ind = ipiv(i - from + 1)
                if (ind < 1 .or. ind > n) then
                    info = i - from + 1
                    call report_runtime_err(SRNAME, IND_OUT_OF_RANGE_ERR_CODE)
                    return
                end if
                if (ind /= i) then
                associate(col1 => A(1 + (i-1)*ldA:), col2 => A(1 + (ind-1)*ldA:))
                    call sswap(m, col1, 1, col2, 1)
                end associate
                    if (present(order)) then
                        call swap_pair(order(i), order(ind))
                    end if
                end if
            end do
        end if
    end subroutine sgepiv

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dgepiv &
    (what, dir, m, n, A, lda, from, to, ipiv, info, order)
    use maria_kinds_mod,    only: &
        WP => DP
    use maria_utils_mod,    only: &
        swap_pair
    use maria_argcheck_mod, only: &
        arg_is_bad,               &
        BAD_IF_LESS,              &
        BAD_IF_MORE
    use maria_reports_mod,  only: &
        report_bad_arg,           &
        report_runtime_err,       &
        IND_OUT_OF_RANGE_ERR_CODE
    !-- Input/output arguments -------------------------------------------------
        character(1), intent(in)                          :: what
        character(1), intent(in)                          :: dir
        integer,      intent(in)                          :: m
        integer,      intent(in)                          :: n
        real(WP),     intent(inout), contiguous           :: A(:)
        integer,      intent(in)                          :: ldA
        integer,      intent(in)                          :: from
        integer,      intent(in)                          :: to
        integer,      intent(in),    contiguous           :: ipiv(:)
        integer,      intent(out)                         :: info
        integer,      intent(inout), contiguous, optional :: order(:)

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DGEPIV'
        logical                 :: rows, forward
        integer                 :: i, inc, begg, endd, ind

    !-- Sanity check -----------------------------------------------------------
        rows = .false.
        forward = .false.
    sanity: block
        info = 0
        select case (what)
            case ('c', 'C')
                rows = .false.
            case ('r', 'R')
                rows = .true.
            case default
                info = -1
                exit sanity
        end select
        select case (dir)
            case ('b', 'B')
                forward = .false.
            case ('f', 'F')
                forward = .true.
            case default
                info = -2
                exit sanity
        end select
        if (arg_is_bad(BAD_IF_LESS, m, 1)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, n, 1)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldA, max(1,m))) then
            info = -6
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, from, 1)) then
            info = -7
            exit sanity
        end if
        if (rows) then
            if (arg_is_bad(BAD_IF_MORE, from, m)) then
                info = -7
                exit sanity
            end if
        else
            if (arg_is_bad(BAD_IF_MORE, from, n)) then
                info = -7
                exit sanity
            end if
        end if
        if (arg_is_bad(BAD_IF_LESS, to, from)) then
            info = -8
            exit sanity
        end if
        if (rows) then
            if (arg_is_bad(BAD_IF_MORE, to, m)) then
                info = -8
                exit sanity
            end if
        else
            if (arg_is_bad(BAD_IF_MORE, to, n)) then
                info = -8
                exit sanity
            end if
        end if
    end block sanity

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (forward) then
            inc = 1
            begg = from
            endd = to
        else
            inc = -1
            begg = to
            endd = from
        end if

        if (rows) then
            do i = begg, endd, inc
                ind = ipiv(i - from + 1)
                if (ind < 1 .or. ind > m) then
                    info = i - from + 1
                    call report_runtime_err(SRNAME, IND_OUT_OF_RANGE_ERR_CODE)
                    return
                end if
                if (ind /= i) then
                    call dswap(n, A(i:), ldA, A(ind:), ldA)
                    if (present(order)) then
                        call swap_pair(order(i), order(ind))
                    end if
                end if
            end do
        else
            do i = begg, endd, inc
                ind = ipiv(i - from + 1)
                if (ind < 1 .or. ind > n) then
                    info = i - from + 1
                    call report_runtime_err(SRNAME, IND_OUT_OF_RANGE_ERR_CODE)
                    return
                end if
                if (ind /= i) then
                associate (col1 => A(1 + (i-1)*ldA:), col2 => A(1 + (ind-1)*ldA:))
                    call dswap(m, col1, 1, col2, 1)
                end associate
                    if (present(order)) then
                        call swap_pair(order(i), order(ind))
                    end if
                end if
            end do
        end if
    end subroutine dgepiv

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine sorcangles &
    (what, m, n, U, ldU, V, ldV, cosines, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,    only: &
        WP => SP
    use maria_constants_mod, only: &
        ZERO => S_ZERO,            &
        ONE => S_ONE
    use maria_argcheck_mod, only: &
        arg_is_bad,               &
        BAD_IF_LESS,              &
        BAD_IF_MORE
    use maria_reports_mod,  only: &
        report_bad_arg,           &
        report_runtime_err,       &
        SVD_FAILED_ERR_CODE
    !-- Input/output arguments -------------------------------------------------
        character(1), intent(in)              :: what
        integer,      intent(in)              :: m
        integer,      intent(in)              :: n
        real(WP),     intent(in),  contiguous :: U(:)
        integer,      intent(in)              :: ldU
        real(WP),     intent(in),  contiguous :: V(:)
        integer,      intent(in)              :: ldV
        real(WP),     intent(out), contiguous :: cosines(:)
        real(WP),     intent(out), contiguous :: work(:)
        integer,      intent(in)              :: lwork
        integer,      intent(out), contiguous :: iwork(:)
        integer,      intent(in)              :: liwork
        integer,      intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SORCANGLES'
        logical                 :: colspan, lw_query
        integer                 :: req_lw, req_Liw, lw_sgesdd, liw_sgesdd, &
                                   mem_core, i_core, i_wrk, lwrk, ifoo(1)
        real(WP)                :: foo(1)

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        colspan = .false.
        info = 0
        select case (what)
            case ('r', 'R')
                colspan = .false.
            case ('c', 'C')
                colspan = .true.
            case default
                info = -1
                exit sanity
        end select
        if (arg_is_bad(BAD_IF_LESS, m, 0)) then
            info = -2
            exit sanity
        end if
        if (colspan) then
            if (arg_is_bad(BAD_IF_LESS, n, 0)) then
                info = -3
                exit sanity
            end if
            if (arg_is_bad(BAD_IF_MORE, n, m)) then
                info = -3
                exit sanity
            end if
        else
            if (arg_is_bad(BAD_IF_LESS, n, m)) then
                info = -3
                exit sanity
            end if
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldV, max(1,m))) then
            info = -7
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0) then
            work(1) = sroundup_lwork(1)
            iwork(1) = 1
            return
        end if
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) then
            mem_core = 0
            exit wrkwrk
        end if

        if (colspan) then
        !-- Memory: core [s] ---------------------------------------------------
            mem_core = n * n
        !-- Query procedures ---------------------------------------------------
            call sgesdd_q('n', n, n, foo, n, foo, foo, 1, foo, 1, foo, -1, ifoo, -1, info)
            lw_sgesdd = int(foo(1))
            liw_sgesdd = ifoo(1)
        else
        !-- Memory: core [s] ---------------------------------------------------
            mem_core = m * m
        !-- Query procedures ---------------------------------------------------
            call sgesdd_q('n', m, m, foo, m, foo, foo, 1, foo, 1, foo, -1, ifoo, -1, info)
            lw_sgesdd = int(foo(1))
            liw_sgesdd = ifoo(1)
        end if

        req_lw  = mem_core + lw_sgesdd
        req_liw = liw_sgesdd

        lw_query = (lwork == -1 .or. liwork == -1)
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -10
            exit wrkwrk
        end if
        if (arg_is_bad(BAD_IF_LESS, liwork, req_liw)) then
            info = -12
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (colspan) then
        !-- Slice workspace ----------------------------------------------------
        !-- |......|.....|
        !--   core   wrk
            i_core = 1
            i_wrk = i_core + mem_core
            lwrk = lwork - i_wrk + 1
        associate(core => work(i_core:), wrk => work(i_wrk:))
            call sgemm('t', 'n', n, n, m, ONE, U, ldU, V, ldV, ZERO, core, n)
            call sgesdd_q('n', n, n, core, n, cosines, foo, 1, foo, 1, wrk, lwrk, iwork, liwork, info)
        end associate
        else
        !-- Slice workspace ----------------------------------------------------
        !-- |......|.....|
        !--   core   wrk
            i_core = 1
            i_wrk = i_core + mem_core
            lwrk = lwork - i_wrk + 1
        associate(core => work(i_core:), wrk => work(i_wrk:))
            call sgemm('n', 't', m, m, n, ONE, U, ldU, V, ldV, ZERO, core, m)
            call sgesdd_q('n', m, m, core, m, cosines, foo, 1, foo, 1, wrk, lwrk, iwork, liwork, info)
        end associate
        end if
        if (info > 0) then
            call report_runtime_err(SRNAME, SVD_FAILED_ERR_CODE)
            return
        end if
    end subroutine sorcangles

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dorcangles &
    (what, m, n, U, ldU, V, ldV, cosines, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,    only: &
        WP => DP
    use maria_constants_mod, only: &
        ZERO => D_ZERO,            &
        ONE => D_ONE
    use maria_argcheck_mod, only: &
        arg_is_bad,               &
        BAD_IF_LESS,              &
        BAD_IF_MORE
    use maria_reports_mod,  only: &
        report_bad_arg,           &
        report_runtime_err,       &
        SVD_FAILED_ERR_CODE
    !-- Input/output arguments -------------------------------------------------
        character(1), intent(in)              :: what
        integer,      intent(in)              :: m
        integer,      intent(in)              :: n
        real(WP),     intent(in),  contiguous :: U(:)
        integer,      intent(in)              :: ldU
        real(WP),     intent(in),  contiguous :: V(:)
        integer,      intent(in)              :: ldV
        real(WP),     intent(out), contiguous :: cosines(:)
        real(WP),     intent(out), contiguous :: work(:)
        integer,      intent(in)              :: lwork
        integer,      intent(out), contiguous :: iwork(:)
        integer,      intent(in)              :: liwork
        integer,      intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DORCANGLES'
        logical                 :: colspan, lw_query
        integer                 :: req_lw, req_Liw, lw_dgesdd, liw_dgesdd, &
                                   mem_core, i_core, i_wrk, lwrk, ifoo(1)
        real(WP)                :: foo(1)

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        colspan = .false.
        info = 0
        select case (what)
            case ('r', 'R')
                colspan = .false.
            case ('c', 'C')
                colspan = .true.
            case default
                info = -1
                exit sanity
        end select
        if (arg_is_bad(BAD_IF_LESS, m, 0)) then
            info = -2
            exit sanity
        end if
        if (colspan) then
            if (arg_is_bad(BAD_IF_LESS, n, 0)) then
                info = -3
                exit sanity
            end if
            if (arg_is_bad(BAD_IF_MORE, n, m)) then
                info = -3
                exit sanity
            end if
        else
            if (arg_is_bad(BAD_IF_LESS, n, m)) then
                info = -3
                exit sanity
            end if
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldV, max(1,m))) then
            info = -7
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0) then
            work(1) = droundup_lwork(1)
            iwork(1) = 1
            return
        end if
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) then
            mem_core = 0
            exit wrkwrk
        end if

        if (colspan) then
        !-- Memory: core [d] ---------------------------------------------------
            mem_core = n * n
        !-- Query procedures ---------------------------------------------------
            call dgesdd_q('n', n, n, foo, n, foo, foo, 1, foo, 1, foo, -1, ifoo, -1, info)
            lw_dgesdd = int(foo(1))
            liw_dgesdd = ifoo(1)
        else
        !-- Memory: core [d] ---------------------------------------------------
            mem_core = m * m
        !-- Query procedures ---------------------------------------------------
            call dgesdd_q('n', m, m, foo, m, foo, foo, 1, foo, 1, foo, -1, ifoo, -1, info)
            lw_dgesdd = int(foo(1))
            liw_dgesdd = ifoo(1)
        end if

        req_lw  = mem_core + lw_dgesdd
        req_liw = liw_dgesdd

        lw_query = (lwork == -1 .or. liwork == -1)
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -10
            exit wrkwrk
        end if
        if (arg_is_bad(BAD_IF_LESS, liwork, req_liw)) then
            info = -12
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (colspan) then
        !-- Slice workspace ----------------------------------------------------
        !-- |......|.....|
        !--   core   wrk
            i_core = 1
            i_wrk = i_core + mem_core
            lwrk = lwork - i_wrk + 1
        associate(core => work(i_core:), wrk => work(i_wrk:))
            call dgemm('t', 'n', n, n, m, ONE, U, ldU, V, ldV, ZERO, core, n)
            call dgesdd_q('n', n, n, core, n, cosines, foo, 1, foo, 1, wrk, lwrk, iwork, liwork, info)
        end associate
        else
        !-- Slice workspace ----------------------------------------------------
        !-- |......|.....|
        !--   core   wrk
            i_core = 1
            i_wrk = i_core + mem_core
            lwrk = lwork - i_wrk + 1
        associate(core => work(i_core:), wrk => work(i_wrk:))
            call dgemm('n', 't', m, m, n, ONE, U, ldU, V, ldV, ZERO, core, m)
            call dgesdd_q('n', m, m, core, m, cosines, foo, 1, foo, 1, wrk, lwrk, iwork, liwork, info)
        end associate
        end if
        if (info > 0) then
            call report_runtime_err(SRNAME, SVD_FAILED_ERR_CODE)
            return
        end if
    end subroutine dorcangles

    !------------------------------------------------------------------------------------------------------------------------

end submodule maria_la_core_sub
