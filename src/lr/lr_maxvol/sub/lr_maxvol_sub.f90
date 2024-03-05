!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the implementation of the [[maria_lr_maxvol_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Implements the [[maria_lr_maxvol_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
submodule (maria_lr_maxvol_mod) maria_lr_maxvol_sub
implicit none (type, external)

contains
    !------------------------------------------------------------------------------------------------------------------------

    module function sgevolume &
    (m, n, r, A, lda, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,   only: &
        WP => SP
    use maria_constants_mod,  only: &
        ZERO => S_ZERO,             &
        ONE => S_ONE,               &
        EPS => S_MACHTOL
    use maria_la_core_mod,    only: &
        slacpy,                     &
        slaset,                     &
        sgetrf,                     &
        sgesdd_q,                   &
        sgemm,                      &
        strsm,                      &
        sdgmm,                      &
        sgepiv,                     &
        sroundup_lwork
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS,                &
        BAD_IF_MORE
    use maria_comparison_mod, only: &
        safe_eq,                    &
        safe_leq
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)                :: m
        integer,  intent(in)                :: n
        integer,  intent(in)                :: r
        real(WP), intent(inout), contiguous :: A(:)
        integer,  intent(in)                :: lda
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: info
        real(WP)                            :: sgevolume

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SGEVOLUME'
        logical                 :: lw_query, determinant
        integer                 :: req_lw, req_liw, lw_proc, liw_proc, i_S, mem_S, &
                                   lwrk, i_wrk, ifoo(1), i
        real(WP)                :: foo(1) 
    
    !-- Default values ---------------------------------------------------------
        sgevolume = ZERO

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
        if (arg_is_bad(BAD_IF_MORE, r, min(m,n))) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, lda, max(1,m))) then
            info = -5
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
        if (info /= 0) exit wrkwrk

        determinant = (m == r .and. n == r)

        if (determinant) then
            !-- Total
            req_lw = 0
            req_liw = r
        else
            !-- Storage: S
            mem_S = min(m, n)
            !-- Query procedures
            call sgesdd_q('n', m, n, foo, lda, foo, foo, m, foo, min(m,n), foo, -1, ifoo, -1, info)
            lw_proc = int(foo(1))
            liw_proc = ifoo(1)
            !-- Total
            req_lw = lw_proc + mem_S
            req_liw = liw_proc
        endif

        lw_query = (lwork == -1 .or. liwork == -1)
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -7
            exit wrkwrk
        end if
        if (arg_is_bad(BAD_IF_LESS, liwork, req_liw)) then
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
        sgevolume = ONE

        if (determinant) then
            call sgetrf(r, r, A, lda, iwork, info)
            do i = 1, r
                sgevolume = sgevolume * abs(A(i + (i-1)*lda))
            end do
        else
        !-- Workspace[s]: |..S..|..wrk..|
            i_S = 1
            i_wrk = i_S + mem_S
            lwrk = lwork - i_wrk + 1
        associate &
        (S => work(i_S:), wrk => work(i_wrk:))
            call sgesdd_q('n', m, n, A, lda, S, foo, m, foo, min(m,n), wrk, lwrk, iwork, liwork, info)
            do i = 1, r
                sgevolume = sgevolume * S(i)
            end do
        end associate
        end if
    end function sgevolume

    !------------------------------------------------------------------------------------------------------------------------

    module function dgevolume &
    (m, n, r, A, lda, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,   only: &
        WP => DP
    use maria_constants_mod,  only: &
        ZERO => D_ZERO,             &
        ONE => D_ONE,               &
        EPS => D_MACHTOL
    use maria_la_core_mod,    only: &
        dlacpy,                     &
        dlaset,                     &
        dgetrf,                     &
        dgesdd_q,                   &
        dgemm,                      &
        dtrsm,                      &
        ddgmm,                      &
        dgepiv,                     &
        droundup_lwork
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS,                &
        BAD_IF_MORE
    use maria_comparison_mod, only: &
        safe_eq,                    &
        safe_leq
    use maria_reports_mod,    only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer,  intent(in)                :: m
        integer,  intent(in)                :: n
        integer,  intent(in)                :: r
        real(WP), intent(inout), contiguous :: A(:)
        integer,  intent(in)                :: lda
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: info
        real(WP)                            :: dgevolume

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DGEVOLUME'
        logical                 :: lw_query, determinant
        integer                 :: req_lw, req_liw, lw_proc, liw_proc, i_S, mem_S, &
                                   lwrk, i_wrk, ifoo(1), i
        real(WP)                :: foo(1) 
    
    !-- Default values ---------------------------------------------------------
        dgevolume = ZERO

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
        if (arg_is_bad(BAD_IF_MORE, r, min(m,n))) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, lda, max(1,m))) then
            info = -5
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
        if (info /= 0) exit wrkwrk

        determinant = (m == r .and. n == r)

        if (determinant) then
            !-- Total
            req_lw = 0
            req_liw = r
        else
            !-- Storage: S
            mem_S = min(m, n)
            !-- Query procedures
            call dgesdd_q('n', m, n, foo, lda, foo, foo, m, foo, min(m,n), foo, -1, ifoo, -1, info)
            lw_proc = int(foo(1))
            liw_proc = ifoo(1)
            !-- Total
            req_lw = lw_proc + mem_S
            req_liw = liw_proc
        endif

        lw_query = (lwork == -1 .or. liwork == -1)
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -7
            exit wrkwrk
        end if
        if (arg_is_bad(BAD_IF_LESS, liwork, req_liw)) then
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
        dgevolume = ONE

        if (determinant) then
            call dgetrf(r, r, A, lda, iwork, info)
            do i = 1, r
                dgevolume = dgevolume * abs(A(i + (i-1)*lda))
            end do
        else
        !-- Workspace[s]: |..S..|..wrk..|
            i_S = 1
            i_wrk = i_S + mem_S
            lwrk = lwork - i_wrk + 1
        associate &
        (S => work(i_S:), wrk => work(i_wrk:))
            call dgesdd_q('n', m, n, A, lda, S, foo, m, foo, min(m,n), wrk, lwrk, iwork, liwork, info)
            do i = 1, r
                dgevolume = dgevolume * S(i)
            end do
        end associate
        end if
    end function dgevolume

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine sgemaxvol_swap_rows &
    (ort, m, n, A, lda, irow, thresh, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,   only: &
        WP => SP
    use maria_constants_mod,  only: &
        ZERO => S_ZERO,             &
        ONE => S_ONE
    use maria_utils_mod,      only: &
        isort,                      &
        arange,                     &
        swap_pair
    use maria_la_core_mod,    only: &
        scopy,                      &
        sger,                       &
        sgeqrf,                     &
        sorgqr,                     &
        sgetrf,                     &
        strsm,                      &
        sgepiv,                     &
        sgenrmc,                    &
        sroundup_lwork
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS,                &
        BAD_IF_MORE
    use maria_comparison_mod, only: &
        safe_leq
    use maria_reports_mod,    only: &
        report_bad_arg,             &
        report_runtime_err,         &
        SINGULAR_MATRIX_ERR_CODE
    !-- Input/output arguments -------------------------------------------------
        character(1), intent(in)            :: ort
        integer,  intent(in)                :: m
        integer,  intent(in)                :: n
        real(WP), intent(inout), contiguous :: A(:)
        integer,  intent(in)                :: lda
        integer,  intent(inout), contiguous :: irow(:)
        real(WP), intent(in)                :: thresh
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SGEMAXVOL_SWAP_ROWS'
        logical                 :: lw_query, orthogonalize
        integer                 :: req_lw, req_liw, lw_proc, i_growth, mem_growth, &
            i_ipiv, mem_ipiv, i_order, mem_order, i_tau, mem_tau, i_u, mem_u, i_v, mem_v, &
            lwrk, i_wrk, i, j, max_pos(2)
        real(WP)                :: foo(1), alpha, maxv 
    
    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        select case (ort)
            case ('n', 'N')
                orthogonalize = .false.
            case ('y', 'Y')
                orthogonalize = .true.
            case default
                orthogonalize = .false.
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
        if (arg_is_bad(BAD_IF_MORE, n, m)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, lda, max(1,m))) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, thresh, ONE)) then
            info = -7
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. m == n) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) exit wrkwrk

        if (orthogonalize) then
            !-- Storage: tau
            mem_tau = n
            !-- Query procedures
            call sgeqrf(m, n, foo, lda, foo, foo, -1, info)
            lw_proc = int(foo(1))
            call sorgqr(m, n, n, foo, lda, foo, foo, -1, info)
            lw_proc = max(lw_proc, int(foo(1)))
        end if

        !-- Storage: growth, u, v, order, ipiv
        mem_growth = 1
        mem_u = m - n
        mem_v = n
        mem_order = m
        mem_ipiv = n
        !-- Total
        req_lw = mem_u + mem_v + mem_growth
        req_liw = mem_order + mem_ipiv
        if (orthogonalize) then
           req_lw = max(req_lw, lw_proc + mem_tau)
        end if

        lw_query = (lwork == -1 .or. liwork == -1)
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -9
            exit wrkwrk
        end if
        if (arg_is_bad(BAD_IF_LESS, liwork, req_liw)) then
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
        if (orthogonalize) then
        !-- Workspace |..tau..|..wrk..|
            i_tau = 1
            i_wrk = i_tau + mem_tau
            lwrk = lwork - i_wrk + 1
        associate &
        (tau => work(i_tau:), wrk => work(i_wrk:))
            call sgeqrf(m, n, A, lda, tau, wrk, lwrk, info)
            call sorgqr(m, n, n, A, lda, tau, wrk, lwrk, info)
        end associate
        end if
        
    !-- Workspace |..growth..|..u..|..v..|..wrk..|  |..order..|..ipiv..|
        i_growth = 1
        i_u = i_growth + mem_growth
        i_v = i_u + mem_u
        i_wrk = i_v + mem_v
        lwrk = lwork - i_wrk + 1
        i_order = 1
        i_ipiv = i_order + mem_order
    associate &
    (growth => work(i_growth), u => work(i_u:), v => work(i_v:), wrk => work(i_wrk:), &
    order => iwork(i_order:), ipiv => iwork(i_ipiv:), B => A(n+1:))
        growth = ONE
        call isort('i', n, irow, info)
        call arange(m, order, 1, 1, info)
        call sgepiv('r', 'f', m, n, A, lda, 1, n, irow, info, order)
        
        !-- Invert the top and apply to bottom
        call sgetrf(n, n, A, lda, ipiv, info)
        if (info > 0) then
            call report_runtime_err(SRNAME, SINGULAR_MATRIX_ERR_CODE)
            return
        end if
        call strsm('r', 'u', 'n', 'n', m-n, n, ONE, A, lda, B, lda)
        call strsm('r', 'l', 'n', 'u', m-n, n, ONE, A, lda, B, lda)
        call sgepiv('c', 'b', m-n, n, B, lda, 1, n, ipiv, info)

        swap_rows: do
            maxv = sgenrmc(m-n, n, B, lda, info, max_pos)
            if (safe_leq(maxv, thresh)) exit swap_rows
            i = max_pos(1)
            j = max_pos(2)
            growth = growth * maxv

            call scopy(m-n, B(1 + (j-1)*lda:), 1, u, 1)
            u(i) = u(i) + 1
            call scopy(n, B(i:), lda, v, 1)
            v(j) = v(j) - 1
            alpha = ONE / B(i + (j-1)*lda)
            call sger(m-n, n, -alpha, u, 1, v, 1, B, lda)

            call swap_pair(order(i+n), order(j))
            ! info = info + 1
        end do swap_rows

        !-- Update the set of chosen rows
        irow(1:n) = order(1:n)
        call isort('i', n, irow, info)
    end associate 
    end subroutine sgemaxvol_swap_rows

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dgemaxvol_swap_rows &
    (ort, m, n, A, lda, irow, thresh, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,   only: &
        WP => DP
    use maria_constants_mod,  only: &
        ZERO => D_ZERO,             &
        ONE => D_ONE
    use maria_utils_mod,      only: &
        isort,                      &
        arange,                     &
        swap_pair
    use maria_la_core_mod,    only: &
        dcopy,                      &
        dger,                       &
        dgeqrf,                     &
        dorgqr,                     &
        dgetrf,                     &
        dtrsm,                      &
        dgepiv,                     &
        dgenrmc,                    &
        droundup_lwork
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS,                &
        BAD_IF_MORE
    use maria_comparison_mod, only: &
        safe_leq
    use maria_reports_mod,    only: &
        report_bad_arg,             &
        report_runtime_err,         &
        SINGULAR_MATRIX_ERR_CODE
    !-- Input/output arguments -------------------------------------------------
        character(1), intent(in)            :: ort
        integer,  intent(in)                :: m
        integer,  intent(in)                :: n
        real(WP), intent(inout), contiguous :: A(:)
        integer,  intent(in)                :: lda
        integer,  intent(inout), contiguous :: irow(:)
        real(WP), intent(in)                :: thresh
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DGEMAXVOL_SWAP_ROWS'
        logical                 :: lw_query, orthogonalize
        integer                 :: req_lw, req_liw, lw_proc, i_growth, mem_growth, &
            i_ipiv, mem_ipiv, i_order, mem_order, i_tau, mem_tau, i_u, mem_u, i_v, mem_v, &
            lwrk, i_wrk, i, j, max_pos(2)
        real(WP)                :: foo(1), alpha, maxv 
    
    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        select case (ort)
            case ('n', 'N')
                orthogonalize = .false.
            case ('y', 'Y')
                orthogonalize = .true.
            case default
                orthogonalize = .false.
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
        if (arg_is_bad(BAD_IF_MORE, n, m)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, lda, max(1,m))) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, thresh, ONE)) then
            info = -7
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. m == n) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) exit wrkwrk

        if (orthogonalize) then
            !-- Storage: tau
            mem_tau = n
            !-- Query procedures
            call dgeqrf(m, n, foo, lda, foo, foo, -1, info)
            lw_proc = int(foo(1))
            call dorgqr(m, n, n, foo, lda, foo, foo, -1, info)
            lw_proc = max(lw_proc, int(foo(1)))
        end if

        !-- Storage: growth, u, v, order, ipiv
        mem_growth = 1
        mem_u = m - n
        mem_v = n
        mem_order = m
        mem_ipiv = n
        !-- Total
        req_lw = mem_u + mem_v + mem_growth
        req_liw = mem_order + mem_ipiv
        if (orthogonalize) then
           req_lw = max(req_lw, lw_proc + mem_tau)
        end if

        lw_query = (lwork == -1 .or. liwork == -1)
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -9
            exit wrkwrk
        end if
        if (arg_is_bad(BAD_IF_LESS, liwork, req_liw)) then
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
        if (orthogonalize) then
        !-- Workspace |..tau..|..wrk..|
            i_tau = 1
            i_wrk = i_tau + mem_tau
            lwrk = lwork - i_wrk + 1
        associate &
        (tau => work(i_tau:), wrk => work(i_wrk:))
            call dgeqrf(m, n, A, lda, tau, wrk, lwrk, info)
            call dorgqr(m, n, n, A, lda, tau, wrk, lwrk, info)
        end associate
        end if
        
    !-- Workspace |..growth..|..u..|..v..|..wrk..|  |..order..|..ipiv..|
        i_growth = 1
        i_u = i_growth + mem_growth
        i_v = i_u + mem_u
        i_wrk = i_v + mem_v
        lwrk = lwork - i_wrk + 1
        i_order = 1
        i_ipiv = i_order + mem_order
    associate &
    (growth => work(i_growth), u => work(i_u:), v => work(i_v:), wrk => work(i_wrk:), &
    order => iwork(i_order:), ipiv => iwork(i_ipiv:), B => A(n+1:))
        growth = ONE
        call isort('i', n, irow, info)
        call arange(m, order, 1, 1, info)
        call dgepiv('r', 'f', m, n, A, lda, 1, n, irow, info, order)
        
        !-- Invert the top and apply to bottom
        call dgetrf(n, n, A, lda, ipiv, info)
        if (info > 0) then
            call report_runtime_err(SRNAME, SINGULAR_MATRIX_ERR_CODE)
            return
        end if
        call dtrsm('r', 'u', 'n', 'n', m-n, n, ONE, A, lda, B, lda)
        call dtrsm('r', 'l', 'n', 'u', m-n, n, ONE, A, lda, B, lda)
        call dgepiv('c', 'b', m-n, n, B, lda, 1, n, ipiv, info)

        swap_rows: do
            maxv = dgenrmc(m-n, n, B, lda, info, max_pos)
            if (safe_leq(maxv, thresh)) exit swap_rows
            i = max_pos(1)
            j = max_pos(2)
            growth = growth * maxv

            associate(x => B(1 + (j-1)*lda:))
                call dcopy(m-n, x, 1, u, 1)
            end associate
            u(i) = u(i) + 1
            associate(x => B(i:))
                call dcopy(n, x, lda, v, 1)
            end associate
            v(j) = v(j) - 1
            alpha = ONE / B(i + (j-1)*lda)
            call dger(m-n, n, -alpha, u, 1, v, 1, B, lda)

            call swap_pair(order(i+n), order(j))
            ! info = info + 1
        end do swap_rows

        !-- Update the set of chosen rows
        irow(1:n) = order(1:n)
        call isort('i', n, irow, info)
    end associate 
    end subroutine dgemaxvol_swap_rows

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine sgemaxvol &
    (m, n, matslc, r, irow, icol, niter, thresh, work, lwork, iwork, liwork, info)
    use maria_kinds_mod, only: &
        WP => SP
    use maria_access_matrix_mod, only: &
        MS => smatslc
    use maria_constants_mod,  only: &
        ZERO => S_ZERO,             &
        ONE => S_ONE
    use maria_utils_mod,      only: &
        isort,                      &
        arange,                     &
        swap_pair
    use maria_la_core_mod,    only: &
        scopy,                      &
        slacpy,                     &
        sger,                       &
        sgeqrf,                     &
        sorgqr,                     &
        sgetrf,                     &
        strsm,                      &
        sgepiv,                     &
        sgenrmc,                    &
        sroundup_lwork
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS,                &
        BAD_IF_MORE
    use maria_comparison_mod, only: &
        safe_leq
    use maria_reports_mod,    only: &
        report_bad_arg,             &
        report_runtime_err,         &
        SINGULAR_MATRIX_ERR_CODE
    !-- Input/output arguments -------------------------------------------------
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        procedure(MS), intent(in), pointer :: matslc
        integer,  intent(in)                :: r
        integer,  intent(inout), contiguous :: irow(:)
        integer,  intent(inout), contiguous :: icol(:)
        integer,  intent(in)                :: niter
        real(WP), intent(in)                :: thresh
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SGEMAXVOL'
        logical                 :: lw_query
        integer                 :: req_lw, req_liw, lw_proc, liw_proc, i_A, mem_A, &
            lwrk, i_wrk, i, j, ifoo(1)
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
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, r, min(m,n))) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, niter, 0)) then
            info = -7
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, thresh, ONE)) then
            info = -8
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) exit wrkwrk

        !-- Storage: A
        mem_A = max(m,n) * r
        !-- Query procedures
        call sgemaxvol_swap_rows('y', m, r, foo, m, ifoo, thresh, foo, -1, ifoo, -1, info)
        lw_proc = int(foo(1))
        liw_proc = ifoo(1)
        call sgemaxvol_swap_rows('y', n, r, foo, n, ifoo, thresh, foo, -1, ifoo, -1, info)
        lw_proc = max(lw_proc, int(foo(1)))
        liw_proc = max(liw_proc, ifoo(1))
        !-- Total
        req_lw = mem_A + lw_proc
        req_liw = liw_proc

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
        
        !-- Workspace |..A..|..wrk..|
        i_A = 1
        i_wrk = i_A + mem_A
        lwrk = lwork - i_wrk + 1
    associate &
    (A => work(i_A:), wrk => work(i_wrk:))
        do i = 1, niter
            !-- Compute columns and update the choice of rows
            do j = 1, r
                call matslc(m, n, 2, icol(j), A(1 + (j-1)*m:), 1, info)
            end do
            call sgemaxvol_swap_rows('y', m, r, A, m, irow, thresh, wrk, lwrk, iwork, liwork, info)
            if (info > 0) then
                ! report singular system
                return
            end if

            !-- Compute rows and update the choice of columns
            do j = 1, r
                call matslc(m, n, 1, irow(j), A(1 + (j-1)*n:), 1, info)
            end do
            call sgemaxvol_swap_rows('y', n, r, A, n, icol, thresh, wrk, lwrk, iwork, liwork, info)
        end do
    end associate 
    end subroutine sgemaxvol

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dgemaxvol &
    (m, n, matslc, r, irow, icol, niter, thresh, work, lwork, iwork, liwork, info)
    use maria_kinds_mod, only: &
        WP => DP
    use maria_access_matrix_mod, only: &
        MS => dmatslc
    use maria_constants_mod,  only: &
        ZERO => D_ZERO,             &
        ONE => D_ONE
    use maria_utils_mod,      only: &
        isort,                      &
        arange,                     &
        swap_pair
    use maria_la_core_mod,    only: &
        dcopy,                      &
        dlacpy,                     &
        dger,                       &
        dgeqrf,                     &
        dorgqr,                     &
        dgetrf,                     &
        dtrsm,                      &
        dgepiv,                     &
        dgenrmc,                    &
        droundup_lwork
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS,                &
        BAD_IF_MORE
    use maria_comparison_mod, only: &
        safe_leq
    use maria_reports_mod,    only: &
        report_bad_arg,             &
        report_runtime_err,         &
        SINGULAR_MATRIX_ERR_CODE
    !-- Input/output arguments -------------------------------------------------
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        procedure(MS), intent(in), pointer :: matslc
        integer,  intent(in)                :: r
        integer,  intent(inout), contiguous :: irow(:)
        integer,  intent(inout), contiguous :: icol(:)
        integer,  intent(in)                :: niter
        real(WP), intent(in)                :: thresh
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DGEMAXVOL'
        logical                 :: lw_query
        integer                 :: req_lw, req_liw, lw_proc, liw_proc, i_A, mem_A, &
            lwrk, i_wrk, i, j, ifoo(1)
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
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, r, min(m,n))) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, niter, 0)) then
            info = -7
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, thresh, ONE)) then
            info = -8
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) exit wrkwrk

        !-- Storage: A
        mem_A = max(m,n) * r
        !-- Query procedures
        call dgemaxvol_swap_rows('y', m, r, foo, m, ifoo, thresh, foo, -1, ifoo, -1, info)
        lw_proc = int(foo(1))
        liw_proc = ifoo(1)
        call dgemaxvol_swap_rows('y', n, r, foo, n, ifoo, thresh, foo, -1, ifoo, -1, info)
        lw_proc = max(lw_proc, int(foo(1)))
        liw_proc = max(liw_proc, ifoo(1))
        !-- Total
        req_lw = mem_A + lw_proc
        req_liw = liw_proc

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
        
        !-- Workspace |..A..|..wrk..|
        i_A = 1
        i_wrk = i_A + mem_A
        lwrk = lwork - i_wrk + 1
    associate &
    (A => work(i_A:), wrk => work(i_wrk:))
        do i = 1, niter
            !-- Compute columns and update the choice of rows
            do j = 1, r
            associate(x => A(1 + (j-1)*m:) )
                call matslc(m, n, 2, icol(j), x, 1, info)
            end associate
            end do
            call dgemaxvol_swap_rows('y', m, r, A, m, irow, thresh, wrk, lwrk, iwork, liwork, info)
            if (info > 0) then
                ! report singular system
                return
            end if

            !-- Compute rows and update the choice of columns
            do j = 1, r
            associate(x => A(1 + (j-1)*n:) )
                call matslc(m, n, 1, irow(j), x, 1, info)
            end associate
            end do
            call dgemaxvol_swap_rows('y', n, r, A, n, icol, thresh, wrk, lwrk, iwork, liwork, info)
        end do
    end associate 
    end subroutine dgemaxvol

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine sgemaxvol_rect_add_rows                                  &
    (ort, m, n, A, lda, k, irow, work, lwork, iwork, liwork, ierr)
    !-- Modules ----------------------------------------------------------------
    use maria_kinds_mod,                                                only:  &
        WP              => SP                                                               
    use maria_constants_mod,                                            only:  &
        ZERO            => S_ZERO,                                             &
        ONE             => S_ONE                                                           
    use maria_utils_mod,                                                only:  &
        isort,                                                                 &
        arange,                                                                &
        swap_pair                                                              
    use maria_la_core_mod,                                              only:  &
        scopy,                                                                 &
        sgemv,                                                                 &
        sswap,                                                                 &
        sger,                                                                  &
        sdot,                                                                  &
        slacpy,                                                                &
        sgeqrf,                                                                &
        sorgqr,                                                                &
        sgetrf,                                                                &
        strsm,                                                                 &
        sgepiv,                                                                &
        sgenrmc,                                                               &
        sroundup_lwork                                                         
    use maria_comparison_mod,                                           only:  &
        safe_leq
    use maria_reports_mod,                                              only:  &
        report_bad_arg,                                                        &
        report_runtime_err,                                                    &
        SINGULAR_MATRIX_ERR_CODE

    !-- Arguments --------------------------------------------------------------
        character(1),       intent(in   )               :: ort              !  1
        integer,            intent(in   )               :: m                !  2 
        integer,            intent(in   )               :: n                !  3
        real(WP),           intent(inout),  contiguous  :: A(:)             !  4
        integer,            intent(in   )               :: lda              !  5
        integer,            intent(in   )               :: k                !  6
        integer,            intent(inout),  contiguous  :: irow(:)          !  7
        real(WP),           intent(  out),  contiguous  :: work(:)          !  8
        integer,            intent(in   )               :: lwork            !  9
        integer,            intent(  out),  contiguous  :: iwork(:)         ! 10
        integer,            intent(in   )               :: liwork           ! 11
        integer,            intent(  out)               :: ierr             ! 12

    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'SGEMAXVOL_RECT_ADD_ROWS'

    !-- Variables --------------------------------------------------------------
        logical     :: lw_query, orthogonalize, bad_arg(12)
        integer     :: req_lw, req_liw, lw_proc, i_growth, mem_growth,         &
                i_ipiv, mem_ipiv, i_order, mem_order, i_tau, mem_tau, i_sqlen, &
                mem_sqlen, lwrk, i_wrk, i, j, mem_B, i_B, nn
        real(WP)    :: foo(1), alpha

    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(1) = (ort /= 'n') .and. (ort /= 'y')
        bad_arg(2) = (m < 0)
        bad_arg(3) = (n < 0) .or. (n > m)
        bad_arg(5) = (lda < max(1,m))
        bad_arg(6) = (k < n) .or. (k > m)

        orthogonalize = (ort == 'y')

    !-- Quick return if possible -----------------------------------------------
    quick: block
        if (any(bad_arg)) exit quick

        if (m == 0 .or. n == 0) return

        if (k == m) then
            call arange(m, irow, 1, 1, ierr)
            return
        end if

        if (k == n) then
            call isort('i', n, irow, ierr)
            return
        end if
    end block quick

    !-- Estimate workspace -----------------------------------------------------
    workspace: block
        if (any(bad_arg)) exit workspace

        if (orthogonalize) then
            !-- Storage: tau
            mem_tau = n
            !-- Query procedures
            call sgeqrf(m, n, foo, lda, foo, foo, -1, ierr)
            lw_proc = int(foo(1))
            call sorgqr(m, n, n, foo, lda, foo, foo, -1, ierr)
            lw_proc = max(lw_proc, int(foo(1)))
        end if

        !-- Storage: growth, sqlen, B, order, ipiv
        mem_growth = 1
        mem_sqlen = m - n
        mem_B = (m - n) * k
        mem_order = m
        mem_ipiv = n
        !-- Total
        req_lw = mem_sqlen + mem_B + mem_growth
        req_liw = mem_order + mem_ipiv
        if (orthogonalize) then
           req_lw = max(req_lw, lw_proc + mem_tau)
        end if

        lw_query = (lwork == -1) .or. (liwork == -1)
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        bad_arg(9) = (lwork < req_lw)
        bad_arg(11) = (liwork < req_liw)
    end block workspace

    !-- Report incorrect argumemts ---------------------------------------------
        if (any(bad_arg)) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------

    !-- Part 1. Replace A with the Q-factor of its QR decomposition
    !---------- for increased stability.

    part1: block
        if (orthogonalize) then
        !-- Slice work:     |..tau..|..wrk..|
            i_tau = 1
            i_wrk = i_tau + mem_tau
            lwrk = lwork - i_wrk + 1
        associate (tau => work(i_tau:), wrk => work(i_wrk:))
            call sgeqrf(m, n, A, lda, tau, wrk, lwrk, ierr)
            call sorgqr(m, n, n, A, lda, tau, wrk, lwrk, ierr)
        end associate
        end if
    end block part1
        
    !-- Part 2. Choose the rows to add that increase the volume the most.

    part2: block
    !-- Slice work:         |..growth..|..sqlen..|..B..|
        i_growth = 1
        i_sqlen = i_growth + mem_growth
        i_B = i_sqlen + mem_sqlen
    !-- Slice iwork:        |..order..|..ipiv..|
        i_order = 1
        i_ipiv = i_order + mem_order
    associate (growth => work(i_growth), sqlen => work(i_sqlen:),              &
    B => work(i_B:), order => iwork(i_order:), ipiv => iwork(i_ipiv:))
        growth = ONE
        call isort('i', n, irow, ierr)
        call arange(m, order, 1, 1, ierr)
        call sgepiv('r', 'f', m, n, A, lda, 1, n, irow, ierr, order)
        
        ! Invert the top and apply to bottom
        call sgetrf(n, n, A, lda, ipiv, ierr)
        if (ierr > 0) then
            call report_runtime_err(SRNAME, SINGULAR_MATRIX_ERR_CODE)
            return
        end if
        call strsm('r', 'u', 'n', 'n', m-n, n, ONE, A, lda, A(n+1:), lda)
        call strsm('r', 'l', 'n', 'u', m-n, n, ONE, A, lda, A(n+1:), lda)
        call sgepiv('c', 'b', m-n, n, A(n+1:), lda, 1, n, ipiv, ierr)
        call slacpy('a', m-n, n, A(n+1:), lda, B, m-n)

        do j = 1, m-n
            sqlen(j) = sdot(n, B(j:), m-n, B(j:), m-n)
        end do

        add_rows: do i = 1, k-n
            nn = n + i - 1
        associate (C => B(i:), new_col => B(i + 1 + nn*(m-n) :))
            j = maxloc(sqlen(i : m-n), dim=1)
            growth = growth * sqrt(1 + sqlen(i+j-1))
            if (j /= 1) then
                call swap_pair(order(n+i), order(n+i+j-1))
                call swap_pair(sqlen(i), sqlen(i+j-1))
                call sswap(nn, C(1:), m-n, C(j:), m-n)
            end if

            alpha = ONE / (1 + sqlen(i))
            call sgemv('n', m-nn-1, nn, alpha, C(2:), m-n, C(1:), m-n, ZERO, new_col, 1)
            call sger(m-nn-1, nn, -ONE, new_col, 1, C(1:), m-n, C(2:), m-n)
            do j = i+1, m-n
                sqlen(j) = sqlen(j) - (1 + sqlen(i)) * new_col(j-i)**2
            end do
        end associate
        end do add_rows

        !-- Update the set of chosen rows
        irow(1:k) = order(1:k)
        call isort('i', k, irow, ierr)
    end associate 
    end block part2
    end subroutine sgemaxvol_rect_add_rows

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dgemaxvol_rect_add_rows                                  &
    (ort, m, n, A, lda, k, irow, work, lwork, iwork, liwork, ierr)
    !-- Modules ----------------------------------------------------------------
    use maria_kinds_mod,                                                only:  &
        WP              => DP                                                               
    use maria_constants_mod,                                            only:  &
        ZERO            => D_ZERO,                                             &
        ONE             => D_ONE                                                           
    use maria_utils_mod,                                                only:  &
        isort,                                                                 &
        arange,                                                                &
        swap_pair                                                              
    use maria_la_core_mod,                                              only:  &
        dcopy,                                                                 &
        dgemv,                                                                 &
        dswap,                                                                 &
        dger,                                                                  &
        ddot,                                                                  &
        dlacpy,                                                                &
        dgeqrf,                                                                &
        dorgqr,                                                                &
        dgetrf,                                                                &
        dtrsm,                                                                 &
        dgepiv,                                                                &
        dgenrmc,                                                               &
        droundup_lwork                                                         
    use maria_comparison_mod,                                           only:  &
        safe_leq
    use maria_reports_mod,                                              only:  &
        report_bad_arg,                                                        &
        report_runtime_err,                                                    &
        SINGULAR_MATRIX_ERR_CODE

    !-- Arguments --------------------------------------------------------------
        character(1),       intent(in   )               :: ort              !  1
        integer,            intent(in   )               :: m                !  2 
        integer,            intent(in   )               :: n                !  3
        real(WP),           intent(inout),  contiguous  :: A(:)             !  4
        integer,            intent(in   )               :: lda              !  5
        integer,            intent(in   )               :: k                !  6
        integer,            intent(inout),  contiguous  :: irow(:)          !  7
        real(WP),           intent(  out),  contiguous  :: work(:)          !  8
        integer,            intent(in   )               :: lwork            !  9
        integer,            intent(  out),  contiguous  :: iwork(:)         ! 10
        integer,            intent(in   )               :: liwork           ! 11
        integer,            intent(  out)               :: ierr             ! 12

    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'DGEMAXVOL_RECT_ADD_ROWS'

    !-- Variables --------------------------------------------------------------
        logical     :: lw_query, orthogonalize, bad_arg(12)
        integer     :: req_lw, req_liw, lw_proc, i_growth, mem_growth,         &
                i_ipiv, mem_ipiv, i_order, mem_order, i_tau, mem_tau, i_sqlen, &
                mem_sqlen, lwrk, i_wrk, i, j, mem_B, i_B, nn
        real(WP)    :: foo(1), alpha

    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(1) = (ort /= 'n') .and. (ort /= 'y')
        bad_arg(2) = (m < 0)
        bad_arg(3) = (n < 0) .or. (n > m)
        bad_arg(5) = (lda < max(1,m))
        bad_arg(6) = (k < n) .or. (k > m)

        orthogonalize = (ort == 'y')

    !-- Quick return if possible -----------------------------------------------
    quick: block
        if (any(bad_arg)) exit quick

        if (m == 0 .or. n == 0) return

        if (k == m) then
            call arange(m, irow, 1, 1, ierr)
            return
        end if

        if (k == n) then
            call isort('i', n, irow, ierr)
            return
        end if
    end block quick

    !-- Estimate workspace -----------------------------------------------------
    workspace: block
        if (any(bad_arg)) exit workspace

        if (orthogonalize) then
            !-- Storage: tau
            mem_tau = n
            !-- Query procedures
            call dgeqrf(m, n, foo, lda, foo, foo, -1, ierr)
            lw_proc = int(foo(1))
            call dorgqr(m, n, n, foo, lda, foo, foo, -1, ierr)
            lw_proc = max(lw_proc, int(foo(1)))
        end if

        !-- Storage: growth, sqlen, B, order, ipiv
        mem_growth = 1
        mem_sqlen = m - n
        mem_B = (m - n) * k
        mem_order = m
        mem_ipiv = n
        !-- Total
        req_lw = mem_sqlen + mem_B + mem_growth
        req_liw = mem_order + mem_ipiv
        if (orthogonalize) then
           req_lw = max(req_lw, lw_proc + mem_tau)
        end if

        lw_query = (lwork == -1) .or. (liwork == -1)
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        bad_arg(9) = (lwork < req_lw)
        bad_arg(11) = (liwork < req_liw)
    end block workspace

    !-- Report incorrect argumemts ---------------------------------------------
        if (any(bad_arg)) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------

    !-- Part 1. Replace A with the Q-factor of its QR decomposition
    !---------- for increased stability.

    part1: block
        if (orthogonalize) then
        !-- Slice work:     |..tau..|..wrk..|
            i_tau = 1
            i_wrk = i_tau + mem_tau
            lwrk = lwork - i_wrk + 1
        associate (tau => work(i_tau:), wrk => work(i_wrk:))
            call dgeqrf(m, n, A, lda, tau, wrk, lwrk, ierr)
            call dorgqr(m, n, n, A, lda, tau, wrk, lwrk, ierr)
        end associate
        end if
    end block part1
        
    !-- Part 2. Choose the rows to add that increase the volume the most.

    part2: block
    !-- Slice work:         |..growth..|..sqlen..|..B..|
        i_growth = 1
        i_sqlen = i_growth + mem_growth
        i_B = i_sqlen + mem_sqlen
    !-- Slice iwork:        |..order..|..ipiv..|
        i_order = 1
        i_ipiv = i_order + mem_order
    associate (growth => work(i_growth), sqlen => work(i_sqlen:),              &
    B => work(i_B:), order => iwork(i_order:), ipiv => iwork(i_ipiv:))
        growth = ONE
        call isort('i', n, irow, ierr)
        call arange(m, order, 1, 1, ierr)
        call dgepiv('r', 'f', m, n, A, lda, 1, n, irow, ierr, order)
        
        ! Invert the top and apply to bottom
        call dgetrf(n, n, A, lda, ipiv, ierr)
        if (ierr > 0) then
            call report_runtime_err(SRNAME, SINGULAR_MATRIX_ERR_CODE)
            return
        end if
        call dtrsm('r', 'u', 'n', 'n', m-n, n, ONE, A, lda, A(n+1:), lda)
        call dtrsm('r', 'l', 'n', 'u', m-n, n, ONE, A, lda, A(n+1:), lda)
        call dgepiv('c', 'b', m-n, n, A(n+1:), lda, 1, n, ipiv, ierr)
        call dlacpy('a', m-n, n, A(n+1:), lda, B, m-n)

        do j = 1, m-n
            sqlen(j) = ddot(n, B(j:), m-n, B(j:), m-n)
        end do

        add_rows: do i = 1, k-n
            nn = n + i - 1
        associate (C => B(i:), new_col => B(i + 1 + nn*(m-n) :))
            j = maxloc(sqlen(i : m-n), dim=1)
            growth = growth * sqrt(1 + sqlen(i+j-1))
            if (j /= 1) then
                call swap_pair(order(n+i), order(n+i+j-1))
                call swap_pair(sqlen(i), sqlen(i+j-1))
                call dswap(nn, C(1:), m-n, C(j:), m-n)
            end if

            alpha = ONE / (1 + sqlen(i))
            call dgemv('n', m-nn-1, nn, alpha, C(2:), m-n, C(1:), m-n, ZERO, new_col, 1)
            call dger(m-nn-1, nn, -ONE, new_col, 1, C(1:), m-n, C(2:), m-n)
            do j = i+1, m-n
                sqlen(j) = sqlen(j) - (1 + sqlen(i)) * new_col(j-i)**2
            end do
        end associate
        end do add_rows

        !-- Update the set of chosen rows
        irow(1:k) = order(1:k)
        call isort('i', k, irow, ierr)
    end associate 
    end block part2
    end subroutine dgemaxvol_rect_add_rows

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine sgemaxvol_rect_swap_rows                                 &
    (ort, m, n, A, lda, k, irow, thresh, work, lwork, iwork, liwork, ierr)
    !-- Modules ----------------------------------------------------------------
    use maria_kinds_mod,                                                only:  &
        WP              => SP                                                               
    use maria_constants_mod,                                            only:  &
        ZERO            => S_ZERO,                                             &
        ONE             => S_ONE                                                           
    use maria_utils_mod,                                                only:  &
        isort,                                                                 &
        arange,                                                                &
        swap_pair                                                              
    use maria_la_core_mod,                                              only:  &
        scopy,                                                                 &
        sgemv,                                                                 &
        sswap,                                                                 &
        sger,                                                                  &
        sdot,                                                                  &
        slacpy,                                                                &
        sgeqrf,                                                                &
        sorgqr,                                                                &
        sormqr,                                                                &
        sgetrf,                                                                &
        strsm,                                                                 &
        sgepiv,                                                                &
        sgenrmc,                                                               &
        slaset,                                                                &
        sscal,                                                                 &
        sroundup_lwork                                                         
    use maria_comparison_mod,                                           only:  &
        safe_leq,                                                              &
        safe_less,                                                             &
        safe_eq
    use maria_reports_mod,                                              only:  &
        report_bad_arg,                                                        &
        report_runtime_err,                                                    &
        SINGULAR_MATRIX_ERR_CODE

    !-- Arguments --------------------------------------------------------------
        character(1),       intent(in   )               :: ort              !  1
        integer,            intent(in   )               :: m                !  2 
        integer,            intent(in   )               :: n                !  3
        real(WP),           intent(inout),  contiguous  :: A(:)             !  4
        integer,            intent(in   )               :: lda              !  5
        integer,            intent(in   )               :: k                !  6
        integer,            intent(inout),  contiguous  :: irow(:)          !  7
        real(WP),           intent(in   )               :: thresh           !  8
        real(WP),           intent(  out),  contiguous  :: work(:)          !  9
        integer,            intent(in   )               :: lwork            ! 10
        integer,            intent(  out),  contiguous  :: iwork(:)         ! 11
        integer,            intent(in   )               :: liwork           ! 12
        integer,            intent(  out)               :: ierr             ! 13

    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'SGEMAXVOL_RECT_SWAP_ROWS'

    !-- Variables --------------------------------------------------------------
        logical     :: lw_query, orthogonalize, bad_arg(13)
        integer     :: req_lw, req_liw, lw_proc, i_growth, mem_growth,         &
                i_ipiv, mem_ipiv, i_order, mem_order, i_tau, mem_tau, i_sqlen, &
                mem_sqlen, lwrk, i_wrk, i, j, mem_B, i_B, nn, i_new, i_old,    &
                i_col, mem_col, i_row, mem_row, i_scr, mem_scr, max_pos(2)
        real(WP)    :: foo(1), alpha, max_scr

    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(1) = (ort /= 'n') .and. (ort /= 'y')
        bad_arg(2) = (m < 0)
        bad_arg(3) = (n < 0) .or. (n > m)
        bad_arg(5) = (lda < max(1,m))
        bad_arg(6) = (k < n) .or. (k > m)
        bad_arg(8) = safe_less(thresh, ONE)

        orthogonalize = (ort == 'y')

    !-- Quick return if possible -----------------------------------------------
    quick: block
        if (any(bad_arg)) exit quick

        if (m == 0 .or. n == 0) return

        if (k == m) then
            call arange(m, irow, 1, 1, ierr)
            return
        end if
    end block quick

    !-- Estimate workspace -----------------------------------------------------
    workspace: block
        if (any(bad_arg)) exit workspace

        req_lw = 0
        req_liw = 0
        if (orthogonalize) then
            !-- Storage: tau
            mem_tau = n
            !-- Query procedures
            call sgeqrf(m, n, foo, lda, foo, foo, -1, ierr)
            lw_proc = int(foo(1))
            call sorgqr(m, n, n, foo, lda, foo, foo, -1, ierr)
            lw_proc = max(lw_proc, int(foo(1)))
            req_lw = mem_tau + lw_proc
        end if

        !-- Storage: growth, sqlen, B, tau, wrk, order
        mem_growth = 1
        mem_sqlen = m
        mem_B = m * k
        mem_tau = n 
        mem_order = m
        !-- Query procedures
        call sgeqrf(k, n, foo, lda, foo, foo, -1, ierr)
        lw_proc = int(foo(1))
        call sormqr('r', 't', m, k, n, foo, lda, foo, foo, m, foo, -1, ierr)
        lw_proc = max(lw_proc, int(foo(1)))
        !-- Total
        req_lw = max(req_lw, mem_sqlen + mem_B + mem_growth + mem_tau + lw_proc)
        req_liw = max(req_liw, mem_order)

        !-- Storage: growth, sqlen, B, scr, col, row
        mem_scr = (m-k) * k
        mem_col = m
        mem_row = k
        !-- Total
        req_lw = max(req_lw, mem_growth + mem_sqlen + mem_B + mem_scr + mem_col + mem_row)

        lw_query = (lwork == -1) .or. (liwork == -1)
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        bad_arg(9) = (lwork < req_lw)
        bad_arg(11) = (liwork < req_liw)
    end block workspace

    !-- Report incorrect argumemts ---------------------------------------------
        if (any(bad_arg)) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------

    !-- Part 1. Replace A with the Q-factor of its QR decomposition
    !---------- for increased stability.

    part1: block
        if (orthogonalize) then
        !-- Slice work:     |..tau..|..wrk..|
            i_tau = 1
            i_wrk = i_tau + mem_tau
            lwrk = lwork - i_wrk + 1
        associate (tau => work(i_tau:), wrk => work(i_wrk:))
            call sgeqrf(m, n, A, lda, tau, wrk, lwrk, ierr)
            call sorgqr(m, n, n, A, lda, tau, wrk, lwrk, ierr)
        end associate
        end if
    end block part1
        
    !-- Part 2. Choose the rows to add that increase the volume the most.-

    part2: block
    !-- Slice work:         |..growth..|..sqlen..|..B..|..tau..|..wrk..|
        i_growth = 1
        i_sqlen = i_growth + mem_growth
        i_B = i_sqlen + mem_sqlen
        i_tau = i_B + mem_B
        i_wrk = i_tau + mem_tau
        lwrk = lwork - i_wrk + 1
    !-- Slice iwork:        |..order..|
        i_order = 1
    associate (growth => work(i_growth), sqlen => work(i_sqlen:),              &
    B => work(i_B:), tau => work(i_tau:), wrk => work(i_wrk:),                 &
    order => iwork(i_order:))
        call isort('i', k, irow, ierr)
        call arange(m, order, 1, 1, ierr)
        call sgepiv('r', 'f', m, n, A, lda, 1, k, irow, ierr, order)
        
        ! Invert the top and apply to bottom
        call slaset('a', m, k, ZERO, ZERO, B, m)
        call slacpy('a', m, n, A, lda, B, m)
        call sgeqrf(k, n, A, lda, tau, wrk, lwrk, ierr)
        do i = 1, n
            if (safe_eq(A(i + (i-1)*lda), ZERO)) then
                ierr = 1
                call report_runtime_err(SRNAME, SINGULAR_MATRIX_ERR_CODE)
                return
            end if
        end do
        call strsm('r', 'u', 'n', 'n', m, n, ONE, A, lda, B, m)
        call sormqr('r', 't', m, k, n, A, lda, tau, B, m, wrk, lwrk, ierr)

        do j = 1, m
            sqlen(j) = sdot(k, B(j:), m, B(j:), m)
        end do
    end associate
    end block part2

    !-- Part 3.
    part3: block
    !-- Slice work:         |..growth..|..sqlen..|..B..|..scr..|..col..|..row..|
        i_growth = 1
        i_sqlen = i_growth + mem_growth
        i_B = i_sqlen + mem_sqlen
        i_scr = i_B + mem_B
        i_col = i_scr + mem_scr
        i_row = i_col + mem_col
    !-- Slice iwork:        |..order..|
        i_order = 1
    associate (growth => work(i_growth), sqlen => work(i_sqlen:),              &
    B => work(i_B:), scr => work(i_scr:), col => work(i_col:),                 &
    row => work(i_row:), order => iwork(i_order:))
        growth = ONE
        swap_rows: do 
            do j = 1, k
                do i = 1, m - k
                    scr(i + (j-1)*(m-k)) = B(k + i + (j-1)*m)**2 &
                            + (ONE + sqlen(k+i)) * (ONE - sqlen(j))
                end do
            end do
            max_scr = sqrt(sgenrmc(m-k, k, scr, m-k, ierr, max_pos))
            if (safe_leq(max_scr, thresh)) exit swap_rows
            i_new = max_pos(1) + k
            i_old = max_pos(2)
            growth = growth * max_scr

            ! Add row I_NEW as if it were the last
            alpha = sqlen(i_new) + ONE
            call scopy(k, B(i_new:), m, row, 1)
            call sgemv('n', m, k, ONE/alpha, B, m, row, 1, ZERO, col, 1)
            call sger(m, k, -ONE, col, 1, row, 1, B, m)
            do j = 1, m
                sqlen(j) = sqlen(j) - alpha * col(j)**2
            end do

            ! Swap I_NEW and I_OLD
            associate (new_row => B(i_new:), old_row => B(i_old:))
                call sswap(k, new_row, m, old_row, m)
            end associate
            call swap_pair(sqlen(i_new), sqlen(i_old))
            call swap_pair(order(i_new), order(i_old))
            call swap_pair(col(i_new), col(i_old))
            associate (old_col =>  B(1 + (i_old - 1)*m:))
                call sswap(m, col, 1, old_col, 1)
            end associate

            ! Remove row I_NEW as if it were the last
            alpha = ONE / (ONE - sqlen(i_new))
            call scopy(k, B(i_new:), m, row, 1)
            call sger(m, k, alpha, col, 1, row, 1, B, m)
            do j = 1, m
                sqlen(j) = sqlen(j) + alpha * col(j)**2
            end do
        end do swap_rows

        irow(1:k) = order(1:k)
        call isort('i', k, irow, ierr)
    end associate
    end block part3
    end subroutine sgemaxvol_rect_swap_rows

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dgemaxvol_rect_swap_rows                                 &
    (ort, m, n, A, lda, k, irow, thresh, work, lwork, iwork, liwork, ierr)
    !-- Modules ----------------------------------------------------------------
    use maria_kinds_mod,                                                only:  &
        WP              => DP                                                               
    use maria_constants_mod,                                            only:  &
        ZERO            => D_ZERO,                                             &
        ONE             => D_ONE                                                           
    use maria_utils_mod,                                                only:  &
        isort,                                                                 &
        arange,                                                                &
        swap_pair                                                              
    use maria_la_core_mod,                                              only:  &
        dcopy,                                                                 &
        dgemv,                                                                 &
        dswap,                                                                 &
        dger,                                                                  &
        ddot,                                                                  &
        dlacpy,                                                                &
        dgeqrf,                                                                &
        dorgqr,                                                                &
        dormqr,                                                                &
        dgetrf,                                                                &
        dtrsm,                                                                 &
        dgepiv,                                                                &
        dgenrmc,                                                               &
        dlaset,                                                                &
        dscal,                                                                 &
        droundup_lwork                                                         
    use maria_comparison_mod,                                           only:  &
        safe_leq,                                                              &
        safe_less,                                                             &
        safe_eq
    use maria_reports_mod,                                              only:  &
        report_bad_arg,                                                        &
        report_runtime_err,                                                    &
        SINGULAR_MATRIX_ERR_CODE

    !-- Arguments --------------------------------------------------------------
        character(1),       intent(in   )               :: ort              !  1
        integer,            intent(in   )               :: m                !  2 
        integer,            intent(in   )               :: n                !  3
        real(WP),           intent(inout),  contiguous  :: A(:)             !  4
        integer,            intent(in   )               :: lda              !  5
        integer,            intent(in   )               :: k                !  6
        integer,            intent(inout),  contiguous  :: irow(:)          !  7
        real(WP),           intent(in   )               :: thresh           !  8
        real(WP),           intent(  out),  contiguous  :: work(:)          !  9
        integer,            intent(in   )               :: lwork            ! 10
        integer,            intent(  out),  contiguous  :: iwork(:)         ! 11
        integer,            intent(in   )               :: liwork           ! 12
        integer,            intent(  out)               :: ierr             ! 13

    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'DGEMAXVOL_RECT_SWAP_ROWS'

    !-- Variables --------------------------------------------------------------
        logical     :: lw_query, orthogonalize, bad_arg(13)
        integer     :: req_lw, req_liw, lw_proc, i_growth, mem_growth,         &
                i_ipiv, mem_ipiv, i_order, mem_order, i_tau, mem_tau, i_sqlen, &
                mem_sqlen, lwrk, i_wrk, i, j, mem_B, i_B, nn, i_new, i_old,    &
                i_col, mem_col, i_row, mem_row, i_scr, mem_scr, max_pos(2)
        real(WP)    :: foo(1), alpha, max_scr

    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(1) = (ort /= 'n') .and. (ort /= 'y')
        bad_arg(2) = (m < 0)
        bad_arg(3) = (n < 0) .or. (n > m)
        bad_arg(5) = (lda < max(1,m))
        bad_arg(6) = (k < n) .or. (k > m)
        bad_arg(8) = safe_less(thresh, ONE)

        orthogonalize = (ort == 'y')

    !-- Quick return if possible -----------------------------------------------
    quick: block
        if (any(bad_arg)) exit quick

        if (m == 0 .or. n == 0) return

        if (k == m) then
            call arange(m, irow, 1, 1, ierr)
            return
        end if
    end block quick

    !-- Estimate workspace -----------------------------------------------------
    workspace: block
        if (any(bad_arg)) exit workspace

        req_lw = 0
        req_liw = 0
        if (orthogonalize) then
            !-- Storage: tau
            mem_tau = n
            !-- Query procedures
            call dgeqrf(m, n, foo, lda, foo, foo, -1, ierr)
            lw_proc = int(foo(1))
            call dorgqr(m, n, n, foo, lda, foo, foo, -1, ierr)
            lw_proc = max(lw_proc, int(foo(1)))
            req_lw = mem_tau + lw_proc
        end if

        !-- Storage: growth, sqlen, B, tau, wrk, order
        mem_growth = 1
        mem_sqlen = m
        mem_B = m * k
        mem_tau = n 
        mem_order = m
        !-- Query procedures
        call dgeqrf(k, n, foo, lda, foo, foo, -1, ierr)
        lw_proc = int(foo(1))
        call dormqr('r', 't', m, k, n, foo, lda, foo, foo, m, foo, -1, ierr)
        lw_proc = max(lw_proc, int(foo(1)))
        !-- Total
        req_lw = max(req_lw, mem_sqlen + mem_B + mem_growth + mem_tau + lw_proc)
        req_liw = max(req_liw, mem_order)

        !-- Storage: growth, sqlen, B, scr, col, row
        mem_scr = (m-k) * k
        mem_col = m
        mem_row = k
        !-- Total
        req_lw = max(req_lw, mem_growth + mem_sqlen + mem_B + mem_scr + mem_col + mem_row)

        lw_query = (lwork == -1) .or. (liwork == -1)
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        bad_arg(9) = (lwork < req_lw)
        bad_arg(11) = (liwork < req_liw)
    end block workspace

    !-- Report incorrect argumemts ---------------------------------------------
        if (any(bad_arg)) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------

    !-- Part 1. Replace A with the Q-factor of its QR decomposition
    !---------- for increased stability.

    part1: block
        if (orthogonalize) then
        !-- Slice work:     |..tau..|..wrk..|
            i_tau = 1
            i_wrk = i_tau + mem_tau
            lwrk = lwork - i_wrk + 1
        associate (tau => work(i_tau:), wrk => work(i_wrk:))
            call dgeqrf(m, n, A, lda, tau, wrk, lwrk, ierr)
            call dorgqr(m, n, n, A, lda, tau, wrk, lwrk, ierr)
        end associate
        end if
    end block part1
        
    !-- Part 2. Choose the rows to add that increase the volume the most.-

    part2: block
    !-- Slice work:         |..growth..|..sqlen..|..B..|..tau..|..wrk..|
        i_growth = 1
        i_sqlen = i_growth + mem_growth
        i_B = i_sqlen + mem_sqlen
        i_tau = i_B + mem_B
        i_wrk = i_tau + mem_tau
        lwrk = lwork - i_wrk + 1
    !-- Slice iwork:        |..order..|
        i_order = 1
    associate (growth => work(i_growth), sqlen => work(i_sqlen:),              &
    B => work(i_B:), tau => work(i_tau:), wrk => work(i_wrk:),                 &
    order => iwork(i_order:))
        call isort('i', k, irow, ierr)
        call arange(m, order, 1, 1, ierr)
        call dgepiv('r', 'f', m, n, A, lda, 1, k, irow, ierr, order)
        
        ! Invert the top and apply to bottom
        call dlaset('a', m, k, ZERO, ZERO, B, m)
        call dlacpy('a', m, n, A, lda, B, m)
        call dgeqrf(k, n, A, lda, tau, wrk, lwrk, ierr)
        do i = 1, n
            if (safe_eq(A(i + (i-1)*lda), ZERO)) then
                ierr = 1
                call report_runtime_err(SRNAME, SINGULAR_MATRIX_ERR_CODE)
                return
            end if
        end do
        call dtrsm('r', 'u', 'n', 'n', m, n, ONE, A, lda, B, m)
        call dormqr('r', 't', m, k, n, A, lda, tau, B, m, wrk, lwrk, ierr)

        do j = 1, m
            sqlen(j) = ddot(k, B(j:), m, B(j:), m)
        end do
    end associate
    end block part2

    !-- Part 3.
    part3: block
    !-- Slice work:         |..growth..|..sqlen..|..B..|..scr..|..col..|..row..|
        i_growth = 1
        i_sqlen = i_growth + mem_growth
        i_B = i_sqlen + mem_sqlen
        i_scr = i_B + mem_B
        i_col = i_scr + mem_scr
        i_row = i_col + mem_col
    !-- Slice iwork:        |..order..|
        i_order = 1
    associate (growth => work(i_growth), sqlen => work(i_sqlen:),              &
    B => work(i_B:), scr => work(i_scr:), col => work(i_col:),                 &
    row => work(i_row:), order => iwork(i_order:))
        growth = ONE
        swap_rows: do 
            do j = 1, k
                do i = 1, m - k
                    scr(i + (j-1)*(m-k)) = B(k + i + (j-1)*m)**2 &
                            + (ONE + sqlen(k+i)) * (ONE - sqlen(j))
                end do
            end do
            max_scr = sqrt(dgenrmc(m-k, k, scr, m-k, ierr, max_pos))
            if (safe_leq(max_scr, thresh)) exit swap_rows
            i_new = max_pos(1) + k
            i_old = max_pos(2)
            growth = growth * max_scr

            ! Add row I_NEW as if it were the last
            alpha = sqlen(i_new) + ONE
            call dcopy(k, B(i_new:), m, row, 1)
            call dgemv('n', m, k, ONE/alpha, B, m, row, 1, ZERO, col, 1)
            call dger(m, k, -ONE, col, 1, row, 1, B, m)
            do j = 1, m
                sqlen(j) = sqlen(j) - alpha * col(j)**2
            end do

            ! Swap I_NEW and I_OLD
            associate (new_row => B(i_new:), old_row => B(i_old:))
                call dswap(k, new_row, m, old_row, m)
            end associate
            call swap_pair(sqlen(i_new), sqlen(i_old))
            call swap_pair(order(i_new), order(i_old))
            call swap_pair(col(i_new), col(i_old))
            associate (old_col =>  B(1 + (i_old - 1)*m:))
                call dswap(m, col, 1, old_col, 1)
            end associate

            ! Remove row I_NEW as if it were the last
            alpha = ONE / (ONE - sqlen(i_new))
            call dcopy(k, B(i_new:), m, row, 1)
            call dger(m, k, alpha, col, 1, row, 1, B, m)
            do j = 1, m
                sqlen(j) = sqlen(j) + alpha * col(j)**2
            end do
        end do swap_rows

        irow(1:k) = order(1:k)
        call isort('i', k, irow, ierr)
    end associate
    end block part3
    end subroutine dgemaxvol_rect_swap_rows

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine sgemaxvol_rect_swap_cols                                 &
    (m, n, A, lda, k, icol, thresh, work, lwork, iwork, liwork, ierr)
    !-- Modules ----------------------------------------------------------------
    use maria_kinds_mod,                                                only:  &
        WP              => SP                                                               
    use maria_constants_mod,                                            only:  &
        ZERO            => S_ZERO,                                             &
        ONE             => S_ONE                                                           
    use maria_utils_mod,                                                only:  &
        isort,                                                                 &
        arange,                                                                &
        swap_pair                                                              
    use maria_la_core_mod,                                              only:  &
        scopy,                                                                 &
        sgemv,                                                                 &
        sswap,                                                                 &
        sger,                                                                  &
        sdot,                                                                  &
        slacpy,                                                                &
        sgeqrf,                                                                &
        sorgqr,                                                                &
        sormqr,                                                                &
        sgetrf,                                                                &
        strsm,                                                                 &
        sgepiv,                                                                &
        sgenrmc,                                                               &
        slaset,                                                                &
        sscal,                                                                 &
        strtri,                                                                &
        strmm,                                                                 &
        strsv,                                                                 &
        srotg,                                                                 &
        srot,                                                                  &
        slarfg,                                                                &
        slarf,                                                                 &
        saxpy,                                                                 &
        strmv,                                                                 &
        sroundup_lwork                                                         
    use maria_comparison_mod,                                           only:  &
        safe_leq,                                                              &
        safe_less,                                                             &
        safe_eq
    use maria_reports_mod,                                              only:  &
        report_bad_arg,                                                        &
        report_runtime_err,                                                    &
        SINGULAR_MATRIX_ERR_CODE

    !-- Arguments --------------------------------------------------------------
        integer,            intent(in   )               :: m                !  1 
        integer,            intent(in   )               :: n                !  2
        real(WP),           intent(inout),  contiguous  :: A(:)             !  3
        integer,            intent(in   )               :: lda              !  4
        integer,            intent(in   )               :: k                !  5
        integer,            intent(inout),  contiguous  :: icol(:)          !  6
        real(WP),           intent(in   )               :: thresh           !  7
        real(WP),           intent(  out),  contiguous  :: work(:)          !  8
        integer,            intent(in   )               :: lwork            !  9
        integer,            intent(  out),  contiguous  :: iwork(:)         ! 10
        integer,            intent(in   )               :: liwork           ! 11
        integer,            intent(  out)               :: ierr             ! 12

    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'SGEMAXVOL_RECT_SWAP_COLS'

    !-- Variables --------------------------------------------------------------
        logical     :: lw_query, orthogonalize, bad_arg(12)
        integer     :: req_lw, req_liw, lw_proc, i_growth, mem_growth,         &
                i_ipiv, mem_ipiv, i_order, mem_order, i_tau, mem_tau, i_nrm2C, &
                mem_nrm2C, lwrk, i_wrk, i, j, mem_nrm2UU, i_nrm2UU, i_new,     &
                i_old, mem_UU, i_UU, i_C, i_scr, mem_scr, max_pos(2),   &
                i_U, i_UUB, i_UUx
        real(WP)    :: foo(1), alph, max_scr, beta, mu, rotcos, rotsin, tauC,  &
                    zeta

    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(1) = (m < 0)
        bad_arg(2) = (n < 0)
        bad_arg(4) = (lda < max(1,m))
        bad_arg(5) = (k < 0) .or. (k > min(m,n))
        bad_arg(7) = safe_less(thresh, ONE)

    !-- Quick return if possible -----------------------------------------------
    quick: block
        if (any(bad_arg)) exit quick

        if (m == 0 .or. n == 0 .or. k == 0) return

        if (k == n) then
            call arange(n, icol, 1, 1, ierr)
            return
        end if
    end block quick

    !-- Estimate workspace -----------------------------------------------------
    workspace: block
        if (any(bad_arg)) exit workspace

        req_lw = 0
        req_liw = 0

        !-- Part 1
        !-- Storage: growth, nrm2UU, nrm2C, UU, tau, order
        mem_growth = 1
        mem_nrm2UU = k
        mem_nrm2C = n - k
        mem_UU = k * k
        mem_tau = k
        mem_order = n
        !-- Query procedures
        call sgeqrf(m, k, foo, lda, foo, foo, -1, ierr)
        lw_proc = int(foo(1))
        call sormqr('l', 't', m, n-k, k, foo, lda, foo, foo, lda, foo, -1, ierr)
        lw_proc = max(lw_proc, int(foo(1)))
        !-- Total
        req_lw = mem_growth + mem_nrm2UU + mem_nrm2C + mem_UU + mem_tau + lw_proc
        req_liw = mem_order

        !-- Part 2
        !-- Storage: growth, nrm2UU, nrm2C, scr, order, ipiv
        mem_scr = k * (n-k)
        mem_ipiv = k - 1
        !-- Total
        req_lw = max(req_lw, mem_growth + mem_nrm2UU + mem_nrm2C + mem_scr)
        req_liw = max(req_liw, mem_order + mem_ipiv)

        lw_query = (lwork == -1) .or. (liwork == -1)
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        bad_arg(9) = (lwork < req_lw)
        bad_arg(11) = (liwork < req_liw)
    end block workspace

    !-- Report incorrect argumemts ---------------------------------------------
        if (any(bad_arg)) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------

    !-- Part 1. Choose the rows to add that increase the volume the most.-

    part2: block
    !-- Slice work:         |..growth..|..nrm2UU..|..nrm2C..|..UU..|..tau..|..wrk..|
        i_growth = 1
        i_nrm2UU = i_growth + mem_growth
        i_nrm2C = i_nrm2UU + mem_nrm2UU
        i_UU = i_nrm2C + mem_nrm2C
        i_tau = i_UU + mem_UU
        i_wrk = i_tau + mem_tau
        lwrk = lwork - i_wrk + 1
    !-- Slice iwork:        |..order..|
        i_order = 1
    !-- Slice A:            _____________
    !----------           |   |         |
    !----------           | U |   UUB   | k
    !----------      A =  |___|_________|
    !----------           |___|____C____| m-k
    !----------             k     n-k
        i_U = 1
        i_UUB = 1 + k * lda
        i_C = k+1 + k * lda
    associate (growth => work(i_growth), nrm2UU => work(i_nrm2UU:),            &
    nrm2C => work(i_nrm2C:), UU => work(i_UU:), tau => work(i_tau:),           &
    wrk => work(i_wrk:), order => iwork(i_order:), U => A(i_U:),               &
    UUB => A(i_UUB:), C => A(i_C:))
        ! Move seelcted columns to the left
        call isort('i', k, icol, ierr)
        call arange(n, order, 1, 1, ierr)
        call sgepiv('c', 'f', m, n, A, lda, 1, k, icol, ierr, order)
        
        ! Find QR of the left m x k block and apply Q^T to the right block
        call sgeqrf(m, k, A, lda, tau, wrk, lwrk, ierr)
        call sormqr('l' ,'t', m, n-k, k, A, lda, tau, UUB, lda, wrk, lwrk, ierr)

        ! Explicitly form the R factor
        call slaset('l', k-1, k-1, ZERO, ZERO, U(2:), lda)

        ! Compute column norms of C
        do i = 1, n-k
        associate (col => C(1 + (i-1)*lda:))
            nrm2C(i) = sdot(m-k, col, 1, col, 1)
        end associate
        end do

        ! Compute row norms of inverted U
        call slacpy('a', k, k, U, lda, UU, k)
        call strtri('u', 'n', k, UU, k, ierr)
        if (ierr > 0) then
            call report_runtime_err(SRNAME, SINGULAR_MATRIX_ERR_CODE)
            return
        end if
        do i = 1, k
        associate (row => UU(i + (i-1)*k:))
            nrm2UU(i) = sdot(k-i+1, row, k, row, k)
        end associate
        end do

        ! Apply inverted U to the top right block
        call strmm('l', 'u', 'n', 'n', k, n-k, ONE, UU, k, UUB, lda)
    end associate
    end block part2

    !-- Part 3.
    part3: block
    !-- Slice work:         |..growth..|..nrm2UU..|..nrm2C..|..scr..|
    !----------                                                wrk
    !----------                                                UUx
        i_growth = 1
        i_nrm2UU = i_growth + mem_growth
        i_nrm2C = i_nrm2UU + mem_nrm2UU
        i_scr = i_nrm2C + mem_nrm2C
        i_UUx = i_scr
        i_wrk = i_scr
        lwrk = lwork - i_wrk + 1
    !-- Slice iwork:        |..order..|..ipiv..|
        i_order = 1
        i_ipiv = i_order + mem_order
    !-- Slice A:            _____________
    !----------           |   |         |
    !----------           | U |   UUB   | k
    !----------      A =  |___|_________|
    !----------           |___|____C____| m-k
    !----------             k     n-k
        i_U = 1
        i_UUB = 1 + k * lda
        i_C = k+1 + k * lda
    associate (growth => work(i_growth), nrm2UU => work(i_nrm2UU:),           &
    nrm2C => work(i_nrm2C:), order => iwork(i_order:), ipiv => iwork(i_ipiv:), &
    U => A(i_U:), UUB => A(i_UUB:), C => A(i_C:))
        growth = ONE
        swap_cols: do 
        !-- Pick column to swap
        associate (scr => work(i_scr:)) 
            do j = 1, n - k
                do i = 1, k
                    scr(i + (j-1)*k) = UUB(i + (j-1)*lda)**2 &
                            + nrm2C(j) * nrm2UU(i)
                end do
            end do
            max_scr = sqrt(sgenrmc(k, n-k, scr, k, ierr, max_pos))
            if (safe_leq(max_scr, thresh)) exit swap_cols
            i_old = max_pos(1)
            i_new = max_pos(2)
            growth = growth * max_scr
        end associate 

            ! Move K+I_NEW to K+1
            if (i_new /= 1) then
                call swap_pair(order(k+1), order(k + i_new))
                associate (old_col => UUB, new_col => UUB(1 + (i_new-1)*lda:))
                    call sswap(m, old_col, 1, new_col, 1)
                end associate
                call swap_pair(nrm2C(1), nrm2C(i_new))
            end if

            ! Move I_OLD to K
            if (i_old /= k) then
                ! Cyclic shift
                do i = 1, k - i_old
                    ipiv(i) = i_old + i
                end do
                call sgepiv('c', 'f', k, k, U, lda, i_old, k-1, ipiv, ierr, order)
                call sgepiv('c', 'f', 1, k, nrm2UU, 1, i_old, k-1, ipiv, ierr)
                call sgepiv('r', 'f', k, n-k, UUB, lda, i_old, k-1, ipiv, ierr)

                ! Hessenberg reduction of the shifted U
                do i = i_old, k-1
                    call srotg(U(i + (i-1)*lda), U(i+1 + (i-1)*lda), rotcos, rotsin)
                    U(i+1 + (i-1)*lda) = ZERO
                    associate (row => U(i + i*lda:), next_row => U(i+1 + i*lda:))
                        call srot(k-i, row, lda, next_row, lda, rotcos, rotsin)
                    end associate
                end do
            end if


            ! Householder reflection of the first column of C
            if (m > k) then
                call slarfg(m-k, C(1), C(2:), 1, tauC)
                alph = C(1)
                C(1) = ONE
                associate (wrk => work(i_wrk:))
                    call slarf('l', m-k, n-k-1, C, 1, tauC, C(1 + lda:), lda, wrk)
                end associate
                C(1) = alph
                call slaset('a', m-k-1, 1, ZERO, ZERO, C(2:), m-k-1)
            end if

            ! Compute angle to rotate the last row of B and the first row of C
            beta = U(k + (k-1)*lda)
            zeta = beta * UUB(k)
            mu = zeta
            if (m > k) then
                call srotg(mu, C(1), rotcos, rotsin)
            end if

            associate (UUX => work(i_UUx:))
            ! Find 'last column' of UU (its scaled top part)
            call scopy(k-1, U(1 + (k-1)*lda:), 1, UUx, 1)
            call strsv('u', 'n', 'n', k-1, U, lda, UUx, 1)

            ! Compute first column of B [y^T zeta]^T via UUy and update U
            call saxpy(k-1, zeta/beta, UUx, 1, UUB, 1) ! stores UUy
            call scopy(k-1, UUB, 1, U(1 + (k-1)*lda:), 1)
            call strmv('u', 'n', 'n', k-1, U, lda, U(1 + (k-1)*lda:), 1)
            U(k + (k-1)*lda) = mu

            ! Update row norms of UU
            do i = 1, k-1
                nrm2UU(i) = nrm2UU(i) - (UUx(i)/beta)**2 + (UUB(i)/mu)**2
            end do
            nrm2UU(k) = ONE / mu**2

            ! Remove the impcate of the previous last column of U from UUB
            call sger(k-1, n-k-1, ONE, UUx, 1, UUB(k+lda:), lda, UUB(1+lda:), lda)

            ! Remove the impact of the first row of C on its column norms 
            ! starting from 2
            if (m > k) then
                do i = 2, n-k
                    nrm2C(i) = nrm2C(i) - C(1 + (i-1)*lda)**2
                end do
            end if

            ! Rotate the last row of B and the first row of C
            call sscal(n-k-1, beta, UUB(k+lda:), lda)
            if (m > k) then
                call srot(n-k-1, UUB(k+lda:), lda, C(1 + lda:), lda, rotcos, rotsin)
                UUB(k) = rotcos * beta / mu
                C(1) = -rotsin * beta
            else
                UUB(k) = beta / mu
            end if
            
            ! Update column norms of C
            if (m > k) then
                nrm2C(1) = C(1)**2
                do i = 2, n-k
                    nrm2C(i) = nrm2C(i) + C(1 + (i-1)*lda)**2
                end do
            end if

            ! Update UUB
            call sscal(n-k-1, ONE/mu, UUB(k+lda:), lda)
            call sger(k-1, n-k-1, -ONE, UUB, 1, UUB(k+lda:), lda, UUB(1+lda:), lda)
            call sscal(k-1, -UUB(k), UUB, 1)
            call saxpy(k-1, ONE, UUx, 1, UUB, 1)
            end associate

            call swap_pair(order(k), order(k+1))
        end do swap_cols

        icol(1:k) = order(1:k)
        call isort('i', k, icol, ierr)
    end associate
    end block part3
    end subroutine sgemaxvol_rect_swap_cols

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dgemaxvol_rect_swap_cols                                 &
    (m, n, A, lda, k, icol, thresh, work, lwork, iwork, liwork, ierr)
    !-- Modules ----------------------------------------------------------------
    use maria_kinds_mod,                                                only:  &
        WP              => DP                                                               
    use maria_constants_mod,                                            only:  &
        ZERO            => D_ZERO,                                             &
        ONE             => D_ONE                                                           
    use maria_utils_mod,                                                only:  &
        isort,                                                                 &
        arange,                                                                &
        swap_pair                                                              
    use maria_la_core_mod,                                              only:  &
        dcopy,                                                                 &
        dgemv,                                                                 &
        dswap,                                                                 &
        dger,                                                                  &
        ddot,                                                                  &
        dlacpy,                                                                &
        dgeqrf,                                                                &
        dorgqr,                                                                &
        dormqr,                                                                &
        dgetrf,                                                                &
        dtrsm,                                                                 &
        dgepiv,                                                                &
        dgenrmc,                                                               &
        dlaset,                                                                &
        dscal,                                                                 &
        dtrtri,                                                                &
        dtrmm,                                                                 &
        dtrsv,                                                                 &
        drotg,                                                                 &
        drot,                                                                  &
        dlarfg,                                                                &
        dlarf,                                                                 &
        daxpy,                                                                 &
        dtrmv,                                                                 &
        droundup_lwork                                                         
    use maria_comparison_mod,                                           only:  &
        safe_leq,                                                              &
        safe_less,                                                             &
        safe_eq
    use maria_reports_mod,                                              only:  &
        report_bad_arg,                                                        &
        report_runtime_err,                                                    &
        SINGULAR_MATRIX_ERR_CODE

    !-- Arguments --------------------------------------------------------------
        integer,            intent(in   )               :: m                !  1 
        integer,            intent(in   )               :: n                !  2
        real(WP),           intent(inout),  contiguous  :: A(:)             !  3
        integer,            intent(in   )               :: lda              !  4
        integer,            intent(in   )               :: k                !  5
        integer,            intent(inout),  contiguous  :: icol(:)          !  6
        real(WP),           intent(in   )               :: thresh           !  7
        real(WP),           intent(  out),  contiguous  :: work(:)          !  8
        integer,            intent(in   )               :: lwork            !  9
        integer,            intent(  out),  contiguous  :: iwork(:)         ! 10
        integer,            intent(in   )               :: liwork           ! 11
        integer,            intent(  out)               :: ierr             ! 12

    !-- Parameters -------------------------------------------------------------
        character(*), parameter :: SRNAME = 'DGEMAXVOL_RECT_SWAP_COLS'

    !-- Variables --------------------------------------------------------------
        logical     :: lw_query, orthogonalize, bad_arg(12)
        integer     :: req_lw, req_liw, lw_proc, i_growth, mem_growth,         &
                i_ipiv, mem_ipiv, i_order, mem_order, i_tau, mem_tau, i_nrm2C, &
                mem_nrm2C, lwrk, i_wrk, i, j, mem_nrm2UU, i_nrm2UU, i_new,     &
                i_old, mem_UU, i_UU, i_C, i_scr, mem_scr, max_pos(2),   &
                i_U, i_UUB, i_UUx
        real(WP)    :: foo(1), alph, max_scr, beta, mu, rotcos, rotsin, tauC,  &
                    zeta

    !-- Sanity check -----------------------------------------------------------
        ierr = 0
        bad_arg = .false.

        bad_arg(1) = (m < 0)
        bad_arg(2) = (n < 0)
        bad_arg(4) = (lda < max(1,m))
        bad_arg(5) = (k < 0) .or. (k > min(m,n))
        bad_arg(7) = safe_less(thresh, ONE)

    !-- Quick return if possible -----------------------------------------------
    quick: block
        if (any(bad_arg)) exit quick

        if (m == 0 .or. n == 0 .or. k == 0) return

        if (k == n) then
            call arange(n, icol, 1, 1, ierr)
            return
        end if
    end block quick

    !-- Estimate workspace -----------------------------------------------------
    workspace: block
        if (any(bad_arg)) exit workspace

        req_lw = 0
        req_liw = 0

        !-- Part 1
        !-- Storage: growth, nrm2UU, nrm2C, UU, tau, order
        mem_growth = 1
        mem_nrm2UU = k
        mem_nrm2C = n - k
        mem_UU = k * k
        mem_tau = k
        mem_order = n
        !-- Query procedures
        call dgeqrf(m, k, foo, lda, foo, foo, -1, ierr)
        lw_proc = int(foo(1))
        call dormqr('l', 't', m, n-k, k, foo, lda, foo, foo, lda, foo, -1, ierr)
        lw_proc = max(lw_proc, int(foo(1)))
        !-- Total
        req_lw = mem_growth + mem_nrm2UU + mem_nrm2C + mem_UU + mem_tau + lw_proc
        req_liw = mem_order

        !-- Part 2
        !-- Storage: growth, nrm2UU, nrm2C, scr, order, ipiv
        mem_scr = k * (n-k)
        mem_ipiv = k - 1
        !-- Total
        req_lw = max(req_lw, mem_growth + mem_nrm2UU + mem_nrm2C + mem_scr)
        req_liw = max(req_liw, mem_order + mem_ipiv)

        lw_query = (lwork == -1) .or. (liwork == -1)
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        bad_arg(9) = (lwork < req_lw)
        bad_arg(11) = (liwork < req_liw)
    end block workspace

    !-- Report incorrect argumemts ---------------------------------------------
        if (any(bad_arg)) then
            ierr = -findloc(bad_arg, .true., dim=1)
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------

    !-- Part 1. Choose the rows to add that increase the volume the most.-

    part2: block
    !-- Slice work:         |..growth..|..nrm2UU..|..nrm2C..|..UU..|..tau..|..wrk..|
        i_growth = 1
        i_nrm2UU = i_growth + mem_growth
        i_nrm2C = i_nrm2UU + mem_nrm2UU
        i_UU = i_nrm2C + mem_nrm2C
        i_tau = i_UU + mem_UU
        i_wrk = i_tau + mem_tau
        lwrk = lwork - i_wrk + 1
    !-- Slice iwork:        |..order..|
        i_order = 1
    !-- Slice A:            _____________
    !----------           |   |         |
    !----------           | U |   UUB   | k
    !----------      A =  |___|_________|
    !----------           |___|____C____| m-k
    !----------             k     n-k
        i_U = 1
        i_UUB = 1 + k * lda
        i_C = k+1 + k * lda
    associate (growth => work(i_growth), nrm2UU => work(i_nrm2UU:),            &
    nrm2C => work(i_nrm2C:), UU => work(i_UU:), tau => work(i_tau:),           &
    wrk => work(i_wrk:), order => iwork(i_order:), U => A(i_U:),               &
    UUB => A(i_UUB:), C => A(i_C:))
        ! Move seelcted columns to the left
        call isort('i', k, icol, ierr)
        call arange(n, order, 1, 1, ierr)
        call dgepiv('c', 'f', m, n, A, lda, 1, k, icol, ierr, order)
        
        ! Find QR of the left m x k block and apply Q^T to the right block
        call dgeqrf(m, k, A, lda, tau, wrk, lwrk, ierr)
        call dormqr('l' ,'t', m, n-k, k, A, lda, tau, UUB, lda, wrk, lwrk, ierr)

        ! Explicitly form the R factor
        associate(x => U(2:))
            call dlaset('l', k-1, k-1, ZERO, ZERO, x, lda)
        end associate

        ! Compute column norms of C
        do i = 1, n-k
        associate (col => C(1 + (i-1)*lda:))
            nrm2C(i) = ddot(m-k, col, 1, col, 1)
        end associate
        end do

        ! Compute row norms of inverted U
        call dlacpy('a', k, k, U, lda, UU, k)
        call dtrtri('u', 'n', k, UU, k, ierr)
        if (ierr > 0) then
            call report_runtime_err(SRNAME, SINGULAR_MATRIX_ERR_CODE)
            return
        end if
        do i = 1, k
        associate (row => UU(i + (i-1)*k:))
            nrm2UU(i) = ddot(k-i+1, row, k, row, k)
        end associate
        end do

        ! Apply inverted U to the top right block
        call dtrmm('l', 'u', 'n', 'n', k, n-k, ONE, UU, k, UUB, lda)
    end associate
    end block part2

    !-- Part 3.
    part3: block
    !-- Slice work:         |..growth..|..nrm2UU..|..nrm2C..|..scr..|
    !----------                                                wrk
    !----------                                                UUx
        i_growth = 1
        i_nrm2UU = i_growth + mem_growth
        i_nrm2C = i_nrm2UU + mem_nrm2UU
        i_scr = i_nrm2C + mem_nrm2C
        i_UUx = i_scr
        i_wrk = i_scr
        lwrk = lwork - i_wrk + 1
    !-- Slice iwork:        |..order..|..ipiv..|
        i_order = 1
        i_ipiv = i_order + mem_order
    !-- Slice A:            _____________
    !----------           |   |         |
    !----------           | U |   UUB   | k
    !----------      A =  |___|_________|
    !----------           |___|____C____| m-k
    !----------             k     n-k
        i_U = 1
        i_UUB = 1 + k * lda
        i_C = k+1 + k * lda
    associate (growth => work(i_growth), nrm2UU => work(i_nrm2UU:),           &
    nrm2C => work(i_nrm2C:), order => iwork(i_order:), ipiv => iwork(i_ipiv:), &
    U => A(i_U:), UUB => A(i_UUB:), C => A(i_C:))
        growth = ONE
        swap_cols: do 
        !-- Pick column to swap
        associate (scr => work(i_scr:)) 
            do j = 1, n - k
                do i = 1, k
                    scr(i + (j-1)*k) = UUB(i + (j-1)*lda)**2 &
                            + nrm2C(j) * nrm2UU(i)
                end do
            end do
            max_scr = sqrt(dgenrmc(k, n-k, scr, k, ierr, max_pos))
            if (safe_leq(max_scr, thresh)) exit swap_cols
            i_old = max_pos(1)
            i_new = max_pos(2)
            growth = growth * max_scr
        end associate 

            ! Move K+I_NEW to K+1
            if (i_new /= 1) then
                call swap_pair(order(k+1), order(k + i_new))
                associate (old_col => UUB, new_col => UUB(1 + (i_new-1)*lda:))
                    call dswap(m, old_col, 1, new_col, 1)
                end associate
                call swap_pair(nrm2C(1), nrm2C(i_new))
            end if

            ! Move I_OLD to K
            if (i_old /= k) then
                ! Cyclic shift
                do i = 1, k - i_old
                    ipiv(i) = i_old + i
                end do
                call dgepiv('c', 'f', k, k, U, lda, i_old, k-1, ipiv, ierr, order)
                call dgepiv('c', 'f', 1, k, nrm2UU, 1, i_old, k-1, ipiv, ierr)
                call dgepiv('r', 'f', k, n-k, UUB, lda, i_old, k-1, ipiv, ierr)

                ! Hessenberg reduction of the shifted U
                do i = i_old, k-1
                    call drotg(U(i + (i-1)*lda), U(i+1 + (i-1)*lda), rotcos, rotsin)
                    U(i+1 + (i-1)*lda) = ZERO
                    associate (row => U(i + i*lda:), next_row => U(i+1 + i*lda:))
                        call drot(k-i, row, lda, next_row, lda, rotcos, rotsin)
                    end associate
                end do
            end if

            ! Householder reflection of the first column of C
            if (m > k) then
                associate(x => C(2:))
                    call dlarfg(m-k, C(1), x, 1, tauC)
                end associate
                alph = C(1)
                C(1) = ONE
                associate (wrk => work(i_wrk:), x => C(1 + lda:) )
                    call dlarf('l', m-k, n-k-1, C, 1, tauC, x, lda, wrk)
                end associate
                C(1) = alph
                associate(x => C(2:))
                    call dlaset('a', m-k-1, 1, ZERO, ZERO, x, m-k-1)
                end associate
            end if

            ! Compute angle to rotate the last row of B and the first row of C
            beta = U(k + (k-1)*lda)
            zeta = beta * UUB(k)
            mu = zeta
            if (m > k) then
                call drotg(mu, C(1), rotcos, rotsin)
            end if

            associate (UUX => work(i_UUx:))
            associate (lastU => U(1 + (k-1)*lda:))
            ! Find 'last column' of UU (its scaled top part)
            call dcopy(k-1, lastU, 1, UUx, 1)
            call dtrsv('u', 'n', 'n', k-1, U, lda, UUx, 1)

            ! Compute first column of B [y^T zeta]^T via UUy and update U
            call daxpy(k-1, zeta/beta, UUx, 1, UUB, 1) ! stores UUy
            call dcopy(k-1, UUB, 1, lastU, 1)
            call dtrmv('u', 'n', 'n', k-1, U, lda, lastU, 1)
            U(k + (k-1)*lda) = mu
            end associate

            ! Update row norms of UU
            do i = 1, k-1
                nrm2UU(i) = nrm2UU(i) - (UUx(i)/beta)**2 + (UUB(i)/mu)**2
            end do
            nrm2UU(k) = ONE / mu**2

            ! Remove the impact of the previous last column of U from UUB
            associate(x => UUB(k+lda:), y => UUB(1+lda:))
                call dger(k-1, n-k-1, ONE, UUx, 1, x, lda, y, lda)
            end associate

            ! Remove the impact of the first row of C on its column norms 
            ! starting from 2
            if (m > k) then
                do i = 2, n-k
                    nrm2C(i) = nrm2C(i) - C(1 + (i-1)*lda)**2
                end do
            end if

            ! Rotate the last row of B and the first row of C
            associate(x => UUB(k+lda:))
                call dscal(n-k-1, beta, x, lda)
            end associate
            if (m > k) then
                associate(x => UUB(k+lda:), y => C(1+lda:))
                    call drot(n-k-1, x, lda, y, lda, rotcos, rotsin)
                end associate
                UUB(k) = rotcos * beta / mu
                C(1) = -rotsin * beta
            else
                UUB(k) = beta / mu
            end if
           
            ! Update column norms of C
            if (m > k) then
                nrm2C(1) = C(1)**2
                do i = 2, n-k
                    nrm2C(i) = nrm2C(i) + C(1 + (i-1)*lda)**2
                end do
            end if

            ! Update UUB
            associate(x => UUB(k+lda:), y => UUB(1+lda:))
                call dscal(n-k-1, ONE/mu, x, lda)
                call dger(k-1, n-k-1, -ONE, UUB, 1, x, lda, y, lda)
            end associate
            call dscal(k-1, -UUB(k), UUB, 1)
            call daxpy(k-1, ONE, UUx, 1, UUB, 1)
            end associate

            call swap_pair(order(k), order(k+1))
        end do swap_cols

        icol(1:k) = order(1:k)
        call isort('i', k, icol, ierr)
    end associate
    end block part3
    end subroutine dgemaxvol_rect_swap_cols

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine sgemaxvol_proj &
    (m, n, matslc, r, kr, irow, icol_short, kc, icol, irow_short, niter, thresh, work, lwork, iwork, liwork, ierr)
    use maria_kinds_mod, only: &
        WP => SP
    use maria_access_matrix_mod, only: &
        MS => smatslc
    use maria_constants_mod,  only: &
        ZERO => S_ZERO,             &
        ONE => S_ONE
    use maria_utils_mod,      only: &
        isort,                      &
        arange,                     &
        swap_pair
    use maria_la_core_mod,    only: &
        scopy,                      &
        slacpy,                     &
        sger,                       &
        sgeqrf,                     &
        sorgqr,                     &
        sgetrf,                     &
        strsm,                      &
        sgepiv,                     &
        sgenrmc,                    &
        sroundup_lwork
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS,                &
        BAD_IF_MORE
    use maria_comparison_mod, only: &
        safe_leq
    use maria_reports_mod,    only: &
        report_bad_arg,             &
        report_runtime_err,         &
        SINGULAR_MATRIX_ERR_CODE
    !-- Input/output arguments -------------------------------------------------
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        procedure(MS), intent(in), pointer :: matslc
        integer,  intent(in)                :: r
        integer,  intent(in)        :: kr
        integer,  intent(inout), contiguous :: irow(:)
        integer, intent(inout), contiguous :: icol_short(:)
        integer, intent(in)         :: kc
        integer,  intent(inout), contiguous :: icol(:)
        integer, intent(inout), contiguous :: irow_short(:)
        integer,  intent(in)                :: niter
        real(WP), intent(in)                :: thresh
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: ierr

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SGEMAXVOL_PROJ'
        logical                 :: lw_query
        integer                 :: req_lw, req_liw, lw_proc, liw_proc, i_A, mem_A, &
            lwrk, i_wrk, i, j, ifoo(1)
        real(WP)                :: foo(1)
    
    !-- Sanity check -----------------------------------------------------------
    sanity: block
        ierr = 0
        if (arg_is_bad(BAD_IF_LESS, m, 0)) then
            ierr = -1
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            ierr = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            ierr = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, r, min(m,n))) then
            ierr = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, kr, r)) then
            ierr = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, kr, m)) then
            ierr = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, kc, r)) then
            ierr = -8
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, kc, n)) then
            ierr = -8
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, niter, 0)) then
            ierr = -11
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, thresh, ONE)) then
            ierr = -12
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (ierr /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. r == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (ierr /= 0) exit wrkwrk

        !-- Storage: A
        mem_A = max(m * kc, n * kr)
        !-- Query procedures
        call sgemaxvol_rect_swap_cols(kr, n, foo, kr, r, ifoo, thresh, foo, -1, ifoo, -1, ierr)
        lw_proc = int(foo(1))
        liw_proc = ifoo(1)
        call sgemaxvol_rect_swap_rows('y', m, r, foo, m, kr, ifoo, thresh, foo, -1, ifoo, -1, ierr)
        lw_proc = max(lw_proc, int(foo(1)))
        liw_proc = max(liw_proc, ifoo(1))
        call sgemaxvol_rect_swap_cols(kc, m, foo, kc, r, ifoo, thresh, foo, -1, ifoo, -1, ierr)
        lw_proc = max(lw_proc, int(foo(1)))
        liw_proc = max(liw_proc, ifoo(1))
        call sgemaxvol_rect_swap_rows('y', n, r, foo, n, kc, ifoo, thresh, foo, -1, ifoo, -1, ierr)
        lw_proc = max(lw_proc, int(foo(1)))
        liw_proc = max(liw_proc, ifoo(1))
        !-- Total
        req_lw = mem_A + lw_proc
        req_liw = liw_proc

        lw_query = (lwork == -1 .or. liwork == -1)
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            ierr = -14
            exit wrkwrk
        end if
        if (arg_is_bad(BAD_IF_LESS, liwork, req_liw)) then
            ierr = -16
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (ierr /= 0) then
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------
        
        !-- Workspace |..A..|..wrk..|
        i_A = 1
        i_wrk = i_A + mem_A
        lwrk = lwork - i_wrk + 1
    associate &
    (A => work(i_A:), wrk => work(i_wrk:))
        ! Find kr rows that contain a dominant kr x r submatrix
        do i = 1, niter
            !-- Compute r columns and update the choice of kr rows
            do j = 1, r
                call matslc(m, n, 2, icol_short(j), A(1 + (j-1)*m:), 1, ierr)
            end do
            call sgemaxvol_rect_swap_rows('y', m, r, A, m, kr, irow, thresh, wrk, lwrk, iwork, liwork, ierr)
            if (ierr > 0) then
                ! report singular system
            end if

            ! Compute kr rows and update the choice of r columns
            do j = 1, kr
                call matslc(m, n, 1, irow(j), A(j:), kr, ierr)
            end do
            call sgemaxvol_rect_swap_cols(kr, n, A, kr, r, icol_short, thresh, wrk, lwrk, iwork, liwork, ierr)
        end do

        do i = 1, niter
            !-- Compute r rows and update the choice of kc columns
            do j = 1, r
                call matslc(m, n, 1, irow_short(j), A(1 + (j-1)*n:), 1, ierr)
            end do
            call sgemaxvol_rect_swap_rows('y', n, r, A, n, kc, icol, thresh, wrk, lwrk, iwork, liwork, ierr)
            if (ierr > 0) then
                ! report singular system
            end if

            ! Compute kc columns and update the choice of r rows
            do j = 1, kc
                call matslc(m, n, 2, icol(j), A(j:), kc, ierr)
            end do
            call sgemaxvol_rect_swap_cols(kc, m, A, kc, r, irow_short, thresh, wrk, lwrk, iwork, liwork, ierr)
        end do
    end associate 
    end subroutine sgemaxvol_proj

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dgemaxvol_proj &
    (m, n, matslc, r, kr, irow, icol_short, kc, icol, irow_short, niter, thresh, work, lwork, iwork, liwork, ierr)
    use maria_kinds_mod, only: &
        WP => DP
    use maria_access_matrix_mod, only: &
        MS => dmatslc
    use maria_constants_mod,  only: &
        ZERO => D_ZERO,             &
        ONE => D_ONE
    use maria_utils_mod,      only: &
        isort,                      &
        arange,                     &
        swap_pair
    use maria_la_core_mod,    only: &
        dcopy,                      &
        dlacpy,                     &
        dger,                       &
        dgeqrf,                     &
        dorgqr,                     &
        dgetrf,                     &
        dtrsm,                      &
        dgepiv,                     &
        dgenrmc,                    &
        droundup_lwork
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS,                &
        BAD_IF_MORE
    use maria_comparison_mod, only: &
        safe_leq
    use maria_reports_mod,    only: &
        report_bad_arg,             &
        report_runtime_err,         &
        SINGULAR_MATRIX_ERR_CODE
    !-- Input/output arguments -------------------------------------------------
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        procedure(MS), intent(in), pointer :: matslc
        integer,  intent(in)                :: r
        integer,  intent(in)        :: kr
        integer,  intent(inout), contiguous :: irow(:)
        integer, intent(inout), contiguous :: icol_short(:)
        integer, intent(in)         :: kc
        integer,  intent(inout), contiguous :: icol(:)
        integer, intent(inout), contiguous :: irow_short(:)
        integer,  intent(in)                :: niter
        real(WP), intent(in)                :: thresh
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: ierr

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DGEMAXVOL_PROJ'
        logical                 :: lw_query
        integer                 :: req_lw, req_liw, lw_proc, liw_proc, i_A, mem_A, &
            lwrk, i_wrk, i, j, ifoo(1)
        real(WP)                :: foo(1)
    
    !-- Sanity check -----------------------------------------------------------
    sanity: block
        ierr = 0
        if (arg_is_bad(BAD_IF_LESS, m, 0)) then
            ierr = -1
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            ierr = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            ierr = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, r, min(m,n))) then
            ierr = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, kr, r)) then
            ierr = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, kr, m)) then
            ierr = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, kc, r)) then
            ierr = -8
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, kc, n)) then
            ierr = -8
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, niter, 0)) then
            ierr = -11
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, thresh, ONE)) then
            ierr = -12
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (ierr /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. r == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (ierr /= 0) exit wrkwrk

        !-- Storage: A
        mem_A = max(m * kc, n * kr)
        !-- Query procedures
        call dgemaxvol_rect_swap_cols(kr, n, foo, kr, r, ifoo, thresh, foo, -1, ifoo, -1, ierr)
        lw_proc = int(foo(1))
        liw_proc = ifoo(1)
        call dgemaxvol_rect_swap_rows('y', m, r, foo, m, kr, ifoo, thresh, foo, -1, ifoo, -1, ierr)
        lw_proc = max(lw_proc, int(foo(1)))
        liw_proc = max(liw_proc, ifoo(1))
        call dgemaxvol_rect_swap_cols(kc, m, foo, kc, r, ifoo, thresh, foo, -1, ifoo, -1, ierr)
        lw_proc = max(lw_proc, int(foo(1)))
        liw_proc = max(liw_proc, ifoo(1))
        call dgemaxvol_rect_swap_rows('y', n, r, foo, n, kc, ifoo, thresh, foo, -1, ifoo, -1, ierr)
        lw_proc = max(lw_proc, int(foo(1)))
        liw_proc = max(liw_proc, ifoo(1))
        !-- Total
        req_lw = mem_A + lw_proc
        req_liw = liw_proc

        lw_query = (lwork == -1 .or. liwork == -1)
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            ierr = -14
            exit wrkwrk
        end if
        if (arg_is_bad(BAD_IF_LESS, liwork, req_liw)) then
            ierr = -16
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (ierr /= 0) then
            call report_bad_arg(SRNAME, -ierr)
            return
        end if

    !-- Executable section -----------------------------------------------------
        
        !-- Workspace |..A..|..wrk..|
        i_A = 1
        i_wrk = i_A + mem_A
        lwrk = lwork - i_wrk + 1
    associate &
    (A => work(i_A:), wrk => work(i_wrk:))
        ! Find kr rows that contain a dominant kr x r submatrix
        do i = 1, niter
            !-- Compute r columns and update the choice of kr rows
            do j = 1, r
            associate(x => A(1 + (j-1)*m:))
                call matslc(m, n, 2, icol_short(j), x, 1, ierr)
            end associate
            end do
            call dgemaxvol_rect_swap_rows('y', m, r, A, m, kr, irow, thresh, wrk, lwrk, iwork, liwork, ierr)
            if (ierr > 0) then
                ! report singular system
                exit
            end if

            ! Compute kr rows and update the choice of r columns
            do j = 1, kr
            associate(x => A(j:))
                call matslc(m, n, 1, irow(j), x, kr, ierr)
            end associate
            end do
            call dgemaxvol_rect_swap_cols(kr, n, A, kr, r, icol_short, thresh, wrk, lwrk, iwork, liwork, ierr)
        end do

        do i = 1, niter
            !-- Compute r rows and update the choice of kc columns
            do j = 1, r
            associate(x => A(1 + (j-1)*n:))
                call matslc(m, n, 1, irow_short(j), x, 1, ierr)
            end associate
            end do
            call dgemaxvol_rect_swap_rows('y', n, r, A, n, kc, icol, thresh, wrk, lwrk, iwork, liwork, ierr)
            if (ierr > 0) then
                ! report singular system
                exit
            end if

            ! Compute kc columns and update the choice of r rows
            do j = 1, kc
            associate(x => A(j:))
                call matslc(m, n, 2, icol(j), x, kc, ierr)
            end associate
            end do
            call dgemaxvol_rect_swap_cols(kc, m, A, kc, r, irow_short, thresh, wrk, lwrk, iwork, liwork, ierr)
        end do
    end associate 
    end subroutine dgemaxvol_proj

    !------------------------------------------------------------------------------------------------------------------------
end submodule maria_lr_maxvol_sub
