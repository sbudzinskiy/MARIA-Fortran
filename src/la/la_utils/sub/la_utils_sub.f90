!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the implementation of the [[maria_la_utils_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Implements the [[maria_la_utils_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
submodule (maria_la_utils_mod) maria_la_utils_sub
implicit none (type, external)

contains
    !------------------------------------------------------------------------------------------------------------------------

    module function sall_close &
    (n, x, incx, y, incy, info, atol, rtol)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        EPS => S_MACHTOL,          &
        ZERO => S_ZERO
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_utils_mod,     only: &
        optional_val,              &
        are_close
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)             :: n
        real(WP), intent(in), contiguous :: x(:)
        integer,  intent(in)             :: incx
        real(WP), intent(in), contiguous :: y(:)
        integer,  intent(in)             :: incy
        integer,  intent(out)            :: info
        real(WP), intent(in), optional   :: atol
        real(WP), intent(in), optional   :: rtol
        logical                          :: sall_close

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SALL_CLOSE'
        integer                 :: i, posx, posy
        real(WP)                :: atol_, rtol_

    !-- Process optional arguments and their default values --------------------
        atol_ = optional_val(EPS, atol)
        rtol_ = optional_val(ZERO, rtol)

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -1
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, atol_, ZERO)) then
            info = -7
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, rtol_, ZERO)) then
            info = -8
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        sall_close = .true.
        if (n == 0) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            sall_close = .false.
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (incx >= 0) then
            posx = 1
        else
            posx = 1 - (n-1) * incx
        end if

        if (incy >= 0) then
            posy = 1
        else
            posy = 1 - (n-1) * incy
        end if

        do i = 1, n
            sall_close = are_close(x(posx), y(posy), info, atol=atol_, rtol=rtol_)
            if (sall_close) then
                posx = posx + incx
                posy = posy + incy
            else
                return
            end if
        end do
    end function sall_close

    !------------------------------------------------------------------------------------------------------------------------

    module function dall_close &
    (n, x, incx, y, incy, info, atol, rtol)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        EPS => D_MACHTOL,          &
        ZERO => D_ZERO
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_utils_mod,     only: &
        optional_val,              &
        are_close
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)             :: n
        real(WP), intent(in), contiguous :: x(:)
        integer,  intent(in)             :: incx
        real(WP), intent(in), contiguous :: y(:)
        integer,  intent(in)             :: incy
        integer,  intent(out)            :: info
        real(WP), intent(in), optional   :: atol
        real(WP), intent(in), optional   :: rtol
        logical                          :: dall_close

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DALL_CLOSE'
        integer                 :: i, posx, posy
        real(WP)                :: atol_, rtol_

    !-- Process optional arguments and their default values --------------------
        atol_ = optional_val(EPS, atol)
        rtol_ = optional_val(ZERO, rtol)

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -1
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, atol_, ZERO)) then
            info = -7
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, rtol_, ZERO)) then
            info = -8
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        dall_close = .true.
        if (n == 0) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            dall_close = .false.
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (incx >= 0) then
            posx = 1
        else
            posx = 1 - (n-1) * incx
        end if

        if (incy >= 0) then
            posy = 1
        else
            posy = 1 - (n-1) * incy
        end if

        do i = 1, n
            dall_close = are_close(x(posx), y(posy), info, atol=atol_, rtol=rtol_)
            if (dall_close) then
                posx = posx + incx
                posy = posy + incy
            else
                return
            end if
        end do
    end function dall_close

    !------------------------------------------------------------------------------------------------------------------------

    module function sall_const &
    (n, x, incx, alpha, info, atol, rtol)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        EPS => S_MACHTOL,          &
        ZERO => S_ZERO
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_utils_mod,     only: &
        optional_val
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)             :: n
        real(WP), intent(in), contiguous :: x(:)
        integer,  intent(in)             :: incx
        real(WP), intent(in)             :: alpha
        integer,  intent(out)            :: info
        real(WP), intent(in), optional   :: atol
        real(WP), intent(in), optional   :: rtol
        logical                          :: sall_const

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SALL_CONST'
        real(WP)                :: atol_, rtol_, foo(1)

    !-- Process optional arguments and their default values --------------------
        atol_ = optional_val(EPS, atol)
        rtol_ = optional_val(ZERO, rtol)

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -1
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, atol_, ZERO)) then
            info = -6
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, rtol_, ZERO)) then
            info = -7
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        sall_const = .true.
        if (n == 0) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            sall_const = .false.
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------       
        foo(1) = alpha
        sall_const = all_close(n, x, incx, foo, 0, info, atol, rtol)
    end function sall_const

    !------------------------------------------------------------------------------------------------------------------------

    module function dall_const &
    (n, x, incx, alpha, info, atol, rtol)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        EPS => D_MACHTOL,          &
        ZERO => D_ZERO
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_utils_mod,     only: &
        optional_val
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)             :: n
        real(WP), intent(in), contiguous :: x(:)
        integer,  intent(in)             :: incx
        real(WP), intent(in)             :: alpha
        integer,  intent(out)            :: info
        real(WP), intent(in), optional   :: atol
        real(WP), intent(in), optional   :: rtol
        logical                          :: dall_const

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DALL_CONST'
        real(WP)                :: atol_, rtol_, foo(1)

    !-- Process optional arguments and their default values --------------------
        atol_ = optional_val(EPS, atol)
        rtol_ = optional_val(ZERO, rtol)

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -1
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, atol_, ZERO)) then
            info = -6
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, rtol_, ZERO)) then
            info = -7
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        dall_const = .true.
        if (n == 0) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            dall_const = .false.
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        foo(1) = alpha
        dall_const = all_close(n, x, incx, foo, 0, info, atol, rtol)
    end function dall_const

    !------------------------------------------------------------------------------------------------------------------------
    
    module function sgeall_close &
    (m, n, A, ldA, B, ldB, info, atol, rtol)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        EPS => S_MACHTOL,          &
        ZERO => S_ZERO
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_utils_mod,     only: &
        optional_val
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)             :: m
        integer,  intent(in)             :: n
        real(WP), intent(in), contiguous :: A(:)
        integer,  intent(in)             :: ldA
        real(WP), intent(in), contiguous :: B(:)
        integer,  intent(in)             :: ldB
        integer,  intent(out)            :: info
        real(WP), intent(in), optional   :: atol
        real(WP), intent(in), optional   :: rtol
        logical                          :: sgeall_close

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SGEALL_CLOSE'
        integer                 :: j
        real(WP)                :: atol_, rtol_

    !-- Process optional arguments and their default values --------------------
        atol_ = optional_val(EPS, atol)
        rtol_ = optional_val(ZERO, rtol)

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
        if (arg_is_bad(BAD_IF_LESS, ldB, max(1,m))) then
            info = -6
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, atol_, ZERO)) then
            info = -8
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, rtol_, ZERO)) then
            info = -9
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        sgeall_close = .true.
        if (m == 0 .or. n == 0) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            sgeall_close = .false.
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (ldA == m .and. ldB == n) then
            sgeall_close = all_close(m*n, A, 1, B, 1, info, atol=atol_, rtol=rtol_)
            return
        end if

        do j = 1, n
        associate(colA => A(1 + (j-1)*ldA:), colB => B(1 + (j-1)*ldB:))
            sgeall_close = all_close(m, colA, 1, colB, 1, info, atol=atol_, rtol=rtol_)
        end associate
            if (.not. sgeall_close) return
        end do
    end function sgeall_close

    !------------------------------------------------------------------------------------------------------------------------

    module function dgeall_close &
    (m, n, A, ldA, B, ldB, info, atol, rtol)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        EPS => D_MACHTOL,          &
        ZERO => D_ZERO
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_utils_mod,     only: &
        optional_val
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)             :: m
        integer,  intent(in)             :: n
        real(WP), intent(in), contiguous :: A(:)
        integer,  intent(in)             :: ldA
        real(WP), intent(in), contiguous :: B(:)
        integer,  intent(in)             :: ldB
        integer,  intent(out)            :: info
        real(WP), intent(in), optional   :: atol
        real(WP), intent(in), optional   :: rtol
        logical                          :: dgeall_close

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DGEALL_CLOSE'
        integer                 :: j
        real(WP)                :: atol_, rtol_

    !-- Process optional arguments and their default values --------------------
        atol_ = optional_val(EPS, atol)
        rtol_ = optional_val(ZERO, rtol)

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
        if (arg_is_bad(BAD_IF_LESS, ldB, max(1,m))) then
            info = -6
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, atol_, ZERO)) then
            info = -8
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, rtol_, ZERO)) then
            info = -9
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        dgeall_close = .true.
        if (m == 0 .or. n == 0) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            dgeall_close = .false.
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (ldA == m .and. ldB == n) then
            dgeall_close = all_close(m*n, A, 1, B, 1, info, atol=atol_, rtol=rtol_)
            return
        end if

        do j = 1, n
        associate(colA => A(1 + (j-1)*ldA:), colB => B(1 + (j-1)*ldB:))
            dgeall_close = all_close(m, colA, 1, colB, 1, info, atol=atol_, rtol=rtol_)
        end associate
            if (.not. dgeall_close) return
        end do
    end function dgeall_close

    !------------------------------------------------------------------------------------------------------------------------

    module function sgeall_const &
    (m, n, A, ldA, alpha, beta, info, atol, rtol)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        EPS => S_MACHTOL,          &
        ZERO => S_ZERO
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_utils_mod,     only: &
        optional_val,              &
        are_close
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)             :: m
        integer,  intent(in)             :: n
        real(WP), intent(in), contiguous :: A(:)
        integer,  intent(in)             :: ldA
        real(WP), intent(in)             :: alpha
        real(WP), intent(in)             :: beta
        integer,  intent(out)            :: info
        real(WP), intent(in), optional   :: atol
        real(WP), intent(in), optional   :: rtol
        logical                          :: sgeall_const

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SGEALL_CONST'
        integer                 :: i, j
        real(WP)                :: val, atol_, rtol_

    !-- Process optional arguments and their default values --------------------
        atol_ = optional_val(EPS, atol)
        rtol_ = optional_val(ZERO, rtol)

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
        if (arg_is_bad(BAD_IF_LESS, atol_, ZERO)) then
            info = -8
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, rtol_, ZERO)) then
            info = -9
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        sgeall_const = .true.
        if (m == 0 .or. n == 0) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            sgeall_const = .false.
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        do j = 1, n
            do i = 1, m
                if (i == j) then
                    val = beta
                else
                    val = alpha
                end if
                
                sgeall_const = are_close(A(i + (j-1)*ldA), val, info, atol=atol_, rtol=rtol_)
                if (.not. sgeall_const) return
            end do
        end do
    end function sgeall_const

    !------------------------------------------------------------------------------------------------------------------------

    module function dgeall_const &
    (m, n, A, ldA, alpha, beta, info, atol, rtol)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        EPS => D_MACHTOL,          &
        ZERO => D_ZERO
    use maria_argcheck_mod,  only: &
        arg_is_bad,                & 
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_utils_mod,     only: &
        optional_val,              &
        are_close
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)             :: m
        integer,  intent(in)             :: n
        real(WP), intent(in), contiguous :: A(:)
        integer,  intent(in)             :: ldA
        real(WP), intent(in)             :: alpha
        real(WP), intent(in)             :: beta
        integer,  intent(out)            :: info
        real(WP), intent(in), optional   :: atol
        real(WP), intent(in), optional   :: rtol
        logical                          :: dgeall_const

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DGEALL_CONST'
        integer                 :: i, j
        real(WP)                :: val, atol_, rtol_

    !-- Process optional arguments and their default values --------------------
        atol_ = optional_val(EPS, atol)
        rtol_ = optional_val(ZERO, rtol)

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
        if (arg_is_bad(BAD_IF_LESS, atol_, ZERO)) then
            info = -8
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, rtol_, ZERO)) then
            info = -9
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        dgeall_const = .true.
        if (m == 0 .or. n == 0) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            dgeall_const = .false.
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        do j = 1, n
            do i = 1, m
                if (i == j) then
                    val = beta
                else
                    val = alpha
                end if
                
                dgeall_const = are_close(A(i + (j-1)*ldA), val, info, atol=atol_, rtol=rtol_)
                if (.not. dgeall_const) return
            end do
        end do
    end function dgeall_const

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine srandort &
    (rng, m, n, Q, ldQ, work, lwork, info)
    use maria_prng_mod,      only: &
        prng
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        ZERO => S_ZERO,            &
        ONE => S_ONE
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS,               &
        BAD_IF_MORE
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_la_core_mod,   only: &
        sroundup_lwork,            &
        sgeqrf,                    &
        sorgqr,                    &
        sgelqf,                    &
        sorglq
    !-- Input/output arguments -------------------------------------------------
        class(prng), intent(in)              :: rng
        integer,     intent(in)              :: m
        integer,     intent(in)              :: n
        real(WP),    intent(out), contiguous :: Q(:)
        integer,     intent(in)              :: ldQ
        real(WP),    intent(out), contiguous :: work(:)
        integer,     intent(in)              :: lwork
        integer,     intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SRANDORT'
        logical                 :: lw_query
        integer                 :: req_lw, lw_sgeqrf, lw_sorgqr, &
                                   mem_tau, i_tau, i_wrk, lwrk, i
        real(WP)                :: foo(1)

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, m, 0)) then
            info = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldQ, max(1,m))) then
            info = -5
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0) then
            work(1) = sroundup_lwork(1)
            return
        end if
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) then
            mem_tau = 0
            exit wrkwrk
        end if
    !-- Memory: tau [s] --------------------------------------------------------
        mem_tau = min(m,n)
    !-- Query procedures -------------------------------------------------------
        if (m >= n) then
            call sgeqrf(m, n, foo, ldQ, foo, foo, -1, info)
            lw_sgeqrf = int(foo(1))
            call sorgqr(m, n, n, foo, ldQ, foo, foo, -1, info)
            lw_sorgqr = int(foo(1))
        else
            call sgelqf(m, n, foo, ldQ, foo, foo, -1, info)
            lw_sgeqrf = int(foo(1))
            call sorglq(m, n, m, foo, ldQ, foo, foo, -1, info)
            lw_sorgqr = int(foo(1))
        end if

        req_lw  = mem_tau + max(lw_sgeqrf, lw_sorgqr)

        lw_query = (lwork == -1)
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -7
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (ldQ == m) then
            call rng%snormal(m*n, Q, ZERO, ONE, info)
        else
            do i = 1, n
                call rng%snormal(m, Q(1 + (i-1)*ldQ:), ZERO, ONE, info)
            end do
        end if
    !-- Slice workspace --------------------------------------------------------
    !-- |.....|.....|
    !--   tau   wrk
        i_tau = 1
        i_wrk = i_tau + mem_tau
        lwrk = lwork - i_wrk + 1
    associate(tau => work(i_tau:), wrk => work(i_wrk:))
        if (m >= n) then
            call sgeqrf(m, n, Q, ldQ, tau, wrk, lwrk, info)
            call sorgqr(m, n, n, Q, ldQ, tau, wrk, lwrk, info)
        else
            call sgelqf(m, n, Q, ldQ, tau, wrk, lwrk, info)
            call sorglq(m, n, m, Q, ldQ, tau, wrk, lwrk, info)
        end if
    end associate
    end subroutine srandort

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine drandort &
    (rng, m, n, Q, ldQ, work, lwork, info)
    use maria_prng_mod,      only: &
        prng
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        ZERO => D_ZERO,            &
        ONE => D_ONE
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS,               &
        BAD_IF_MORE
    use maria_reports_mod,   only: &
        report_bad_arg
    use maria_la_core_mod,   only: &
        droundup_lwork,            &
        dgeqrf,                    &
        dorgqr,                    &
        dgelqf,                    &
        dorglq
    !-- Input/output arguments -------------------------------------------------
        class(prng), intent(in)              :: rng
        integer,     intent(in)              :: m
        integer,     intent(in)              :: n
        real(WP),    intent(out), contiguous :: Q(:)
        integer,     intent(in)              :: ldQ
        real(WP),    intent(out), contiguous :: work(:)
        integer,     intent(in)              :: lwork
        integer,     intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DRANDORT'
        logical                 :: lw_query
        integer                 :: req_lw, lw_dgeqrf, lw_dorgqr, &
                                   mem_tau, i_tau, i_wrk, lwrk, i
        real(WP)                :: foo(1)

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, m, 0)) then
            info = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldQ, max(1,m))) then
            info = -5
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0) then
            work(1) = real(1, WP)
            return
        end if
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) then
            mem_tau = 0 
            exit wrkwrk
        end if
    !-- Memory: tau [s] --------------------------------------------------------
        mem_tau = min(m,n)
    !-- Query procedures -------------------------------------------------------
        if (m >= n) then
            call dgeqrf(m, n, foo, ldQ, foo, foo, -1, info)
            lw_dgeqrf = int(foo(1))
            call dorgqr(m, n, n, foo, ldQ, foo, foo, -1, info)
            lw_dorgqr = int(foo(1))
        else
            call dgelqf(m, n, foo, ldQ, foo, foo, -1, info)
            lw_dgeqrf = int(foo(1))
            call dorglq(m, n, m, foo, ldQ, foo, foo, -1, info)
            lw_dorgqr = int(foo(1))
        end if

        req_lw  = mem_tau + max(lw_dgeqrf, lw_dorgqr)

        lw_query = (lwork == -1)
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -7
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (ldQ == m) then
            call rng%dnormal(m*n, Q, ZERO, ONE, info)
        else
            do i = 1, n
                call rng%dnormal(m, Q(1 + (i-1)*ldQ:), ZERO, ONE, info)
            end do
        end if
    !-- Slice workspace --------------------------------------------------------
    !-- |.....|.....|
    !--   tau   wrk
        i_tau = 1
        i_wrk = i_tau + mem_tau
        lwrk = lwork - i_wrk + 1
    associate(tau => work(i_tau:), wrk => work(i_wrk:))
        if (m >= n) then
            call dgeqrf(m, n, Q, ldQ, tau, wrk, lwrk, info)
            call dorgqr(m, n, n, Q, ldQ, tau, wrk, lwrk, info)
        else
            call dgelqf(m, n, Q, ldQ, tau, wrk, lwrk, info)
            call dorglq(m, n, m, Q, ldQ, tau, wrk, lwrk, info)
        end if
    end associate
    end subroutine drandort

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine srandsvd &
    (rng, m, n, A, ldA, S, work, lwork, info)
    use maria_prng_mod,      only: &
        prng
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
        sgeqrf,                    &
        sorgqr,                    &
        sdgmm,                     &
        sgemm
    !-- Input/output arguments -------------------------------------------------
        class(prng), intent(in)              :: rng
        integer,     intent(in)              :: m
        integer,     intent(in)              :: n
        real(WP),    intent(out), contiguous :: A(:)
        integer,     intent(in)              :: ldA
        real(WP),    intent(in),  contiguous :: S(:)
        real(WP),    intent(out), contiguous :: work(:)
        integer,     intent(in)              :: lwork
        integer,     intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SRANDSVD'
        logical                 :: lw_query
        integer                 :: req_lw, lw_srandort_mk, lw_srandort_nk, &
                                   mem_U, i_U, mem_V, i_V, i_wrk, lwrk, k
        real(WP)                :: foo(1)

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
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

        if (m == 0 .or. n == 0) then
            work(1) = real(1, WP)
            return
        end if
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) then
            mem_U = 0
            mem_V = 0
            exit wrkwrk
        end if

        k = min(m,n)
    !-- Memory [s]: U, V -------------------------------------------------------
        mem_U = m * k
        mem_V = n * k
    !-- Query procedures -------------------------------------------------------
        call srandort(rng, m, k, foo, m, foo, -1, info)
        lw_srandort_mk = int(foo(1))
        call srandort(rng, n, k, foo, n, foo, -1, info)
        lw_srandort_nk = int(foo(1))

        req_lw  = mem_U + mem_V + max(lw_srandort_mk, lw_srandort_nk)

        lw_query = (lwork == -1)
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -8
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
        !-- |...|...|.....|
        !--   U   V   wrk
        i_U = 1
        i_V = i_U + mem_U
        i_wrk = i_V + mem_V
        lwrk = lwork - i_wrk + 1
    associate(U => work(i_U:), V => work(i_V:), wrk => work(i_wrk:))
        call srandort(rng, m, k, U, m, wrk, lwrk, info)
        call srandort(rng, n, k, V, n, wrk, lwrk, info)
        call sdgmm('r', m, k, U, m, S, 1, info)
        call sgemm('n', 't', m, n, k, ONE, U, m, V, n, ZERO, A, ldA)
    end associate
    end subroutine srandsvd

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine drandsvd &
    (rng, m, n, A, ldA, S, work, lwork, info)
    use maria_prng_mod,      only: &
        prng
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
        dgeqrf,                    &
        dorgqr,                    &
        ddgmm,                     &
        dgemm
    !-- Input/output arguments -------------------------------------------------
        class(prng), intent(in)              :: rng
        integer,     intent(in)              :: m
        integer,     intent(in)              :: n
        real(WP),    intent(out), contiguous :: A(:)
        integer,     intent(in)              :: ldA
        real(WP),    intent(in),  contiguous :: S(:)
        real(WP),    intent(out), contiguous :: work(:)
        integer,     intent(in)              :: lwork
        integer,     intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DRANDSVD'
        logical                 :: lw_query
        integer                 :: req_lw, lw_drandort_mk, lw_drandort_nk, &
                                   mem_U, i_U, mem_V, i_V, i_wrk, lwrk, k
        real(WP)                :: foo(1)

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
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

        if (m == 0 .or. n == 0) then
            work(1) = real(1, WP)
            return
        end if
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) then
            mem_U = 0
            mem_V = 0
            exit wrkwrk
        end if

        k = min(m,n)
    !-- Memory [s]: U, V -------------------------------------------------------
        mem_U = m * k
        mem_V = n * k
    !-- Query procedures -------------------------------------------------------
        call drandort(rng, m, k, foo, m, foo, -1, info)
        lw_drandort_mk = int(foo(1))
        call drandort(rng, n, k, foo, n, foo, -1, info)
        lw_drandort_nk = int(foo(1))

        req_lw  = mem_U + mem_V + max(lw_drandort_mk, lw_drandort_nk)

        lw_query = (lwork == -1)
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -8
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
        !-- |...|...|.....|
        !--   U   V   wrk
        i_U = 1
        i_V = i_U + mem_U
        i_wrk = i_V + mem_V
        lwrk = lwork - i_wrk + 1
    associate(U => work(i_U:), V => work(i_V:), wrk => work(i_wrk:))
        call drandort(rng, m, k, U, m, wrk, lwrk, info)
        call drandort(rng, n, k, V, n, wrk, lwrk, info)
        call ddgmm('r', m, k, U, m, S, 1, info)
        call dgemm('n', 't', m, n, k, ONE, U, m, V, n, ZERO, A, ldA)
    end associate
    end subroutine drandsvd

    !------------------------------------------------------------------------------------------------------------------------

end submodule maria_la_utils_sub
