!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the implementation of the [[maria_lr_cross_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Implements the [[maria_lr_cross_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
submodule (maria_lr_cross_mod) maria_lr_cross_sub
implicit none (type, external)

contains
    !------------------------------------------------------------------------------------------------------------------------

    module subroutine smatcross_top &
    (m, n, nc, cols, ldc, nr, rows, ldr, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
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
        report_bad_arg,             &
        report_runtime_err,         &
        SINGULAR_MATRIX_ERR_CODE
    !-- Input/output arguments -------------------------------------------------
        integer,       intent(in)                :: m
        integer,       intent(in)                :: n
        integer,       intent(in)                :: nc
        real(WP),      intent(inout), contiguous :: cols(:)
        integer,       intent(in)                :: ldc
        integer,       intent(in)                :: nr
        real(WP),      intent(inout), contiguous :: rows(:)
        integer,       intent(in)                :: ldr
        integer,       intent(in)                :: r
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
        character(*), parameter :: SRNAME = 'SMATCROSS_TOP'
        logical                 :: lw_query, square
        integer                 :: req_lw, req_liw, lw_proc, liw_proc, i_A, mem_A, &
                                   mem_S, i_S, mem_tmpU, i_tmpU, mem_tmpVT, i_tmpVT, &
                                   mem_ipiv, lwrk, i_wrk, ifoo(1)
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
        if (arg_is_bad(BAD_IF_LESS, nc, 0)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, nc, n)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldc, max(1,m))) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, nr, 0)) then
            info = -6
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, nr, m)) then
            info = -6
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldr, max(1,nr))) then
            info = -8
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            info = -9
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, r, min(nc,nr))) then
            info = -9
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -11
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,r))) then
            info = -13
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. nc == 0 .or. nr == 0 .or. r == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) exit wrkwrk

        square = (nc == r .and. nr == r)

        if (square) then
            !-- Storage: ipiv
            mem_ipiv = r
            !-- Total
            req_lw = 0
            req_liw = mem_ipiv
        else
            !-- Storage: A, S, tmpU, tmpVT
            mem_A = nr * nc
            mem_S = min(nr, nc)
            mem_tmpU = nr * min(nr, nc)
            mem_tmpVT = min(nr, nc) * nc
            !-- Query procedures
            call sgesdd_q('s', nr, nc, foo, nr, foo, foo, nr, foo, min(nr,nc), foo, -1, ifoo, -1, info)
            lw_proc = int(foo(1))
            liw_proc = ifoo(1)
            !-- Total
            req_lw = lw_proc + mem_A + mem_S + mem_tmpU + mem_tmpVT
            req_liw = liw_proc
        endif

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

    if (square) then
    !-- Workspace[i] |..ipiv..|
    associate &
    (ipiv => iwork, bot => U(r+1:))
        call slacpy('a', m, r, cols, ldc, U, ldu)
        call slacpy('a', r, n, rows, ldr, VT, ldvt)
        call sgetrf(r, r, U, ldu, ipiv, info)
        if (info > 0) then
            call report_runtime_err(SRNAME, SINGULAR_MATRIX_ERR_CODE)
            return
        end if
        call strsm('r', 'u', 'n', 'n', m-r, r, ONE, U, ldu, bot, ldu)
        call strsm('r', 'l', 'n', 'u', m-r, r, ONE, U, ldu, bot, ldu)
        call sgepiv('c', 'b', m-r, r, bot, ldu, 1, r, ipiv, info)
        call slaset('a', r, r, ZERO, ONE, U, ldu)
    end associate
    else
    !-- Workspace[s]: |..A..|..S..|..tmpU..|..tmpVT..|..wrk..|
        i_A = 1
        i_S = i_A + mem_A
        i_tmpU = i_S + mem_S
        i_tmpVT = i_tmpU + mem_tmpU
        i_wrk = i_tmpVT + mem_tmpVT
        lwrk = lwork - i_wrk + 1
    associate &
    (A => work(i_A:), S => work(i_S:), tmpU => work(i_tmpU:), tmpVT => work(i_tmpVT:), wrk => work(i_wrk:))
        call slacpy('a', nr, nc, cols, ldc, A, nr)
        call sgesdd_q('s', nr, nc, A, nr, S, tmpU, nr, tmpVT, min(nr,nc), wrk, lwrk, iwork, liwork, info)
        if (safe_eq(S(r), ZERO) .or. safe_leq(S(r), S(1) * EPS)) then
            info = 1
            call report_runtime_err(SRNAME, SINGULAR_MATRIX_ERR_CODE)
            return
        end if
        S(1:r) = ONE / S(1:r)
        call sdgmm('r', nr, r, tmpU, nr, S, 1, info)
        call sgemm('n', 't', m, r, nc, ONE, cols, ldc, tmpVT, min(nr,nc), ZERO, U, ldu)
        call sgemm('t', 'n', r, n, nr, ONE, tmpU, nr, rows, ldr, ZERO, VT, ldvt)
    end associate
    end if
    end subroutine smatcross_top

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dmatcross_top &
    (m, n, nc, cols, ldc, nr, rows, ldr, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
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
        report_bad_arg,             &
        report_runtime_err,         &
        SINGULAR_MATRIX_ERR_CODE
    !-- Input/output arguments -------------------------------------------------
        integer,       intent(in)                :: m
        integer,       intent(in)                :: n
        integer,       intent(in)                :: nc
        real(WP),      intent(inout), contiguous :: cols(:)
        integer,       intent(in)                :: ldc
        integer,       intent(in)                :: nr
        real(WP),      intent(inout), contiguous :: rows(:)
        integer,       intent(in)                :: ldr
        integer,       intent(in)                :: r
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
        character(*), parameter :: SRNAME = 'DMATCROSS_TOP'
        logical                 :: lw_query, square
        integer                 :: req_lw, req_liw, lw_proc, liw_proc, i_A, mem_A, &
                                   mem_S, i_S, mem_tmpU, i_tmpU, mem_tmpVT, i_tmpVT, &
                                   mem_ipiv, lwrk, i_wrk, ifoo(1)
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
        if (arg_is_bad(BAD_IF_LESS, nc, 0)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, nc, n)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldc, max(1,m))) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, nr, 0)) then
            info = -6
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, nr, m)) then
            info = -6
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldr, max(1,nr))) then
            info = -8
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            info = -9
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, r, min(nc,nr))) then
            info = -9
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -11
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,r))) then
            info = -13
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. nc == 0 .or. nr == 0 .or. r == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) exit wrkwrk

        square = (nc == r .and. nr == r)

        if (square) then
            !-- Storage: ipiv
            mem_ipiv = r
            !-- Total
            req_lw = 0
            req_liw = mem_ipiv
        else
            !-- Storage: A, S, tmpU, tmpVT
            mem_A = nr * nc
            mem_S = min(nr, nc)
            mem_tmpU = nr * min(nr, nc)
            mem_tmpVT = min(nr, nc) * nc
            !-- Query procedures
            call dgesdd_q('s', nr, nc, foo, nr, foo, foo, nr, foo, min(nr,nc), foo, -1, ifoo, -1, info)
            lw_proc = int(foo(1))
            liw_proc = ifoo(1)
            !-- Total
            req_lw = lw_proc + mem_A + mem_S + mem_tmpU + mem_tmpVT
            req_liw = liw_proc
        endif

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

    if (square) then
    !-- Workspace[i] |..ipiv..|
    associate &
    (ipiv => iwork, bot => U(r+1:))
        call dlacpy('a', m, r, cols, ldc, U, ldu)
        call dlacpy('a', r, n, rows, ldr, VT, ldvt)
        call dgetrf(r, r, U, ldu, ipiv, info)
        if (info > 0) then
            call report_runtime_err(SRNAME, SINGULAR_MATRIX_ERR_CODE)
            return
        end if
        call dtrsm('r', 'u', 'n', 'n', m-r, r, ONE, U, ldu, bot, ldu)
        call dtrsm('r', 'l', 'n', 'u', m-r, r, ONE, U, ldu, bot, ldu)
        call dgepiv('c', 'b', m-r, r, bot, ldu, 1, r, ipiv, info)
        call dlaset('a', r, r, ZERO, ONE, U, ldu)
    end associate
    else
    !-- Workspace[s]: |..A..|..S..|..tmpU..|..tmpVT..|..wrk..|
        i_A = 1
        i_S = i_A + mem_A
        i_tmpU = i_S + mem_S
        i_tmpVT = i_tmpU + mem_tmpU
        i_wrk = i_tmpVT + mem_tmpVT
        lwrk = lwork - i_wrk + 1
    associate &
    (A => work(i_A:), S => work(i_S:), tmpU => work(i_tmpU:), tmpVT => work(i_tmpVT:), wrk => work(i_wrk:))
        call dlacpy('a', nr, nc, cols, ldc, A, nr)
        call dgesdd_q('s', nr, nc, A, nr, S, tmpU, nr, tmpVT, min(nr,nc), wrk, lwrk, iwork, liwork, info)
        if (safe_eq(S(r), ZERO) .or. safe_leq(S(r), S(1) * EPS)) then
            info = 1
            call report_runtime_err(SRNAME, SINGULAR_MATRIX_ERR_CODE)
            return
        end if
        S(1:r) = ONE / S(1:r)
        call ddgmm('r', nr, r, tmpU, nr, S, 1, info)
        call dgemm('n', 't', m, r, nc, ONE, cols, ldc, tmpVT, min(nr,nc), ZERO, U, ldu)
        call dgemm('t', 'n', r, n, nr, ONE, tmpU, nr, rows, ldr, ZERO, VT, ldvt)
    end associate
    end if
    end subroutine dmatcross_top

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine smatcross &
    (m, n, fun, nc, icol, nr, irow, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,   only: &
        WP => SP
    use maria_access_matrix_mod, only: &
        MS => smatslc
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
    use maria_utils_mod,      only: &
        sort
    use maria_reports_mod,    only: &
        report_bad_arg,             &
        report_runtime_err,         &
        SINGULAR_MATRIX_ERR_CODE
    !-- Input/output arguments -------------------------------------------------
        integer,       intent(in)                :: m
        integer,       intent(in)                :: n
        procedure(MS), intent(in),    pointer    :: fun
        integer,       intent(in)                :: nc
        integer,       intent(inout), contiguous :: icol(:)
        integer,       intent(in)                :: nr
        integer,       intent(inout), contiguous :: irow(:)
        integer,       intent(in)                :: r
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
        character(*), parameter :: SRNAME = 'SMATCROSS'
        logical                 :: lw_query
        integer                 :: req_lw, req_liw, lw_proc, liw_proc, i_cols, mem_cols, &
                                   mem_rows, i_rows, lwrk, i_wrk, ifoo(1), i
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
        if (arg_is_bad(BAD_IF_LESS, nc, 0)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, nc, n)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, nr, 0)) then
            info = -6
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, nr, m)) then
            info = -6
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            info = -8
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, r, min(nc,nr))) then
            info = -8
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -10
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,r))) then
            info = -12
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. nc == 0 .or. nr == 0 .or. r == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) exit wrkwrk
        
        !-- Storage: cols, rows
        mem_cols = m * nc
        mem_rows = nr * n
        !-- Query procedures
        call smatcross_top(m, n, nc, foo, m, nr, foo, nr, r, foo, ldu, foo, ldvt, foo, -1, ifoo, -1, info)
        lw_proc = int(foo(1))
        liw_proc = ifoo(1)

        req_lw = lw_proc + mem_cols + mem_rows
        req_liw = liw_proc

        lw_query = (lwork == -1 .or. liwork == -1)
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -14
            exit wrkwrk
        end if
        if (arg_is_bad(BAD_IF_LESS, liwork, req_liw)) then
            info = -16
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        
    !-- Workspace: |..cols..|..rows..|..wrk..|
        i_cols = 1
        i_rows = i_cols + mem_cols
        i_wrk = i_rows + mem_rows
        lwrk = lwork - i_wrk + 1
    associate &
    (cols => work(i_cols:), rows => work(i_rows:), wrk => work(i_wrk:))
        call sort('i', nr, irow, info)
        call sort('i', nc, icol, info)
        do i = 1, nr
            call fun(m, n, 1, irow(i), rows(i:), nr, info)
        end do
        do i = 1, nc
            call fun(m, n, 2, icol(i), cols(1 + (i-1)*m:), 1, info)
        end do
        call sgepiv('r', 'f', m, nc, cols, m, 1, nr, irow, info)
        call smatcross_top(m, n, nc, cols, m, nr, rows, nr, r, U, ldu, VT, ldvt, wrk, lwrk, iwork, liwork, info)
        if (info > 0) then
            call report_runtime_err(SRNAME, SINGULAR_MATRIX_ERR_CODE)
            return
        end if
        call sgepiv('r', 'b', m, r, U, ldu, 1, nr, irow, info)
    end associate
    end subroutine smatcross

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dmatcross &
    (m, n, fun, nc, icol, nr, irow, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,   only: &
        WP => DP
    use maria_access_matrix_mod, only: &
        MS => dmatslc
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
    use maria_utils_mod,      only: &
        sort
    use maria_reports_mod,    only: &
        report_bad_arg,             &
        report_runtime_err,         &
        SINGULAR_MATRIX_ERR_CODE
    !-- Input/output arguments -------------------------------------------------
        integer,       intent(in)                :: m
        integer,       intent(in)                :: n
        procedure(MS), intent(in),    pointer    :: fun
        integer,       intent(in)                :: nc
        integer,       intent(inout), contiguous :: icol(:)
        integer,       intent(in)                :: nr
        integer,       intent(inout), contiguous :: irow(:)
        integer,       intent(in)                :: r
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
        character(*), parameter :: SRNAME = 'DMATCROSS'
        logical                 :: lw_query
        integer                 :: req_lw, req_liw, lw_proc, liw_proc, i_cols, mem_cols, &
                                   mem_rows, i_rows, lwrk, i_wrk, ifoo(1), i
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
        if (arg_is_bad(BAD_IF_LESS, nc, 0)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, nc, n)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, nr, 0)) then
            info = -6
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, nr, m)) then
            info = -6
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, r, 0)) then
            info = -8
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, r, min(nc,nr))) then
            info = -8
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info = -10
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,r))) then
            info = -12
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. nc == 0 .or. nr == 0 .or. r == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info /= 0) exit wrkwrk
        
        !-- Storage: cols, rows
        mem_cols = m * nc
        mem_rows = nr * n
        !-- Query procedures
        call dmatcross_top(m, n, nc, foo, m, nr, foo, nr, r, foo, ldu, foo, ldvt, foo, -1, ifoo, -1, info)
        lw_proc = int(foo(1))
        liw_proc = ifoo(1)

        req_lw = lw_proc + mem_cols + mem_rows
        req_liw = liw_proc

        lw_query = (lwork == -1 .or. liwork == -1)
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            iwork(1) = req_liw
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info = -14
            exit wrkwrk
        end if
        if (arg_is_bad(BAD_IF_LESS, liwork, req_liw)) then
            info = -16
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        
    !-- Workspace: |..cols..|..rows..|..wrk..|
        i_cols = 1
        i_rows = i_cols + mem_cols
        i_wrk = i_rows + mem_rows
        lwrk = lwork - i_wrk + 1
    associate &
    (cols => work(i_cols:), rows => work(i_rows:), wrk => work(i_wrk:))
        call sort('i', nr, irow, info)
        call sort('i', nc, icol, info)
        do i = 1, nr
        associate(x => rows(i:))
            call fun(m, n, 1, irow(i), x, nr, info)
        end associate
        end do
        do i = 1, nc
        associate(x => cols(1 + (i-1)*m:))
            call fun(m, n, 2, icol(i), x, 1, info)
        end associate
        end do
        call dgepiv('r', 'f', m, nc, cols, m, 1, nr, irow, info)
        call dmatcross_top(m, n, nc, cols, m, nr, rows, nr, r, U, ldu, VT, ldvt, wrk, lwrk, iwork, liwork, info)
        if (info > 0) then
            call report_runtime_err(SRNAME, SINGULAR_MATRIX_ERR_CODE)
            return
        end if
        call dgepiv('r', 'b', m, r, U, ldu, 1, nr, irow, info)
    end associate
    end subroutine dmatcross

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine smatcross_aca &
    (m, n, matval, col_order, maxr, niter_rook, r, U, ldu, VT, ldvt, irow, icol, work, lwork, info, rtol)
    use maria_kinds_mod,   only: &
        WP => SP
    use maria_constants_mod,  only: &
        ZERO => S_ZERO,             &
        ONE => S_ONE
    use maria_comparison_mod, only: &
        safe_less,                  &
        safe_eq
    use maria_access_matrix_mod, only: &
        MV => smatval
    use maria_la_core_mod,    only: &
        sgenrmc,                    &
        snrm2,                      &
        sdot,                       &
        sgemv,                      &
        sroundup_lwork
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS,                &
        BAD_IF_MORE
    use maria_lr_la_mod,      only: &
        slrval
    use maria_reports_mod,    only: &
        report_bad_arg,             &
        report_runtime_err,         &
        IND_OUT_OF_RANGE_ERR_CODE,  &
        SINGULAR_MATRIX_ERR_CODE
    !-- Input/output arguments -------------------------------------------------
        integer,       intent(in)                :: m
        integer,       intent(in)                :: n
        procedure(MV), intent(in),    pointer    :: matval
        integer,       intent(in),    contiguous :: col_order(:)
        integer,       intent(in)                :: maxr
        integer,       intent(in)                :: niter_rook
        integer,       intent(out)               :: r
        real(WP),      intent(out),   contiguous :: U(:)
        integer,       intent(in)                :: ldu
        real(WP),      intent(out),   contiguous :: VT(:)
        integer,       intent(in)                :: ldvt
        integer,       intent(out),   contiguous :: irow(:)
        integer,       intent(out),   contiguous :: icol(:)
        real(WP),      intent(out),   contiguous :: work(:)
        integer,       intent(in)                :: lwork
        integer,       intent(out)               :: info(2)
        real(WP),      intent(in),    optional   :: rtol

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'SMATCROSS_ACA'
        logical                 :: lw_query, end_loop
        integer                 :: req_lw, cur_icol, cur_irow, cur_icol_ind, i, j, &
                                   i_uu,  i_vv, pos(2)
        real(WP)                :: total_norm2, new_norm, appr_val, pivot, &
                                   maxv, nrmu, nrmv

    !-- Default values ---------------------------------------------------------
        r = 0

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, m, 0)) then
            info(1) = -1
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info(1) = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, maxr, 0)) then
            info(1) = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, maxr, min(m,n))) then
            info(1) = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, niter_rook, 0)) then
            info(1) = -6
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info(1) = -9
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,maxr))) then
            info(1) = -11
            exit sanity
        end if
        if (present(rtol)) then
            if (arg_is_bad(BAD_IF_LESS, rtol, ZERO)) then
                info(1) = -17
                exit sanity
            end if
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info(1) /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. maxr == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info(1) /= 0) exit wrkwrk
        
        if (present(rtol)) then
            req_lw = 2 * maxr
        else
            req_lw = 0
        end if

        lw_query = (lwork == -1)
        if (lw_query) then
            work(1) = sroundup_lwork(req_lw)
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info(1) = -15
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info(1) /= 0) then
            call report_bad_arg(SRNAME, -info(1))
            return
        end if

    !-- Executable section -----------------------------------------------------
        cur_icol_ind = 0
        total_norm2 = ZERO

        if (present(rtol)) then
            i_uu = 1
            i_vv = 1 + maxr
        end if

        increase_rank: do
            next_col: do
                cur_icol_ind = mod(cur_icol_ind, n) + 1
                if (all(icol(1:r) /= col_order(cur_icol_ind))) exit next_col
            end do next_col
            cur_icol = col_order(cur_icol_ind)
            if (cur_icol < 1 .or. cur_icol > n) then
                info = 1
                call report_runtime_err(SRNAME, IND_OUT_OF_RANGE_ERR_CODE)
                return
            end if

        associate &
        (col => U(1 + r*ldu:), row => VT(1 + r:))
            better_col: do i = 1, niter_rook
                do j = 1, m
                    appr_val = slrval(m, n, j, cur_icol, r, U, ldu, VT, ldvt, info(1))
                    col(j) = matval(m, n, j, cur_icol, info(1)) - appr_val
                end do
                maxv = sgenrmc(m, 1, col, m, info(1), pos)
                cur_irow = pos(1)

                do j = 1, n
                    appr_val = slrval(m, n, cur_irow, j, r, U, ldu, VT, ldvt, info(1))
                    row(1 + (j-1)*ldvt) = matval(m, n, cur_irow, j, info(1)) - appr_val
                end do
                maxv = sgenrmc(1, n, row, ldvt, info(1), pos)
                if (pos(2) == cur_icol) then
                    exit better_col
                else
                    cur_icol = pos(2)
                    info(2) = info(2) + 1
                end if
            end do better_col

            find_pivot: do j = 1, m
                appr_val = slrval(m, n, j, cur_icol, r, U, ldu, VT, ldvt, info(1))
                col(j) = matval(m, n, j, cur_icol, info(1)) - appr_val
            end do find_pivot
            maxv = sgenrmc(m, 1, col, m, info(1), pos)
            cur_irow = pos(1)

        !-- Rank-1 update
            appr_val = slrval(m, n, cur_irow, cur_icol, r, U, ldu, VT, ldvt, info(1))
            pivot = matval(m, n, cur_irow, cur_icol, info(1)) - appr_val
            if (safe_eq(pivot, ZERO)) then
                info = 2
                call report_runtime_err(SRNAME, SINGULAR_MATRIX_ERR_CODE)
                return
            end if
            do j = 1, n
                appr_val = slrval(m, n, cur_irow, j, r, U, ldu, VT, ldvt, info(1))
                row(1 + (j-1)*ldvt) = matval(m, n, cur_irow, j, info(1)) - appr_val
                row(1 + (j-1)*ldvt) = row(1 + (j-1)*ldvt) / pivot
            end do

            irow(r+1) = cur_irow
            icol(r+1) = cur_icol

        !-- Stopping criteria
            end_loop = (r+1 == maxr)
            if (present(rtol)) then
                nrmu = snrm2(m, col, 1)
                nrmv = snrm2(n, row, ldvt)
                new_norm = nrmu * nrmv
                total_norm2 = total_norm2 + new_norm**2
                if (r > 0) then
                associate &
                (uu => work(i_uu:), vv => work(i_vv:))
                    call sgemv('t', m, r, ONE, U, ldu, col, 1, ZERO, uu, 1)
                    call sgemv('n', r, n, ONE, VT, ldvt, row, ldvt, ZERO, vv, 1)
                    total_norm2 = total_norm2 + 2 * sdot(r, uu, 1, vv, 1)
                end associate
                end if
                end_loop = end_loop .or. safe_less(new_norm, rtol * sqrt(total_norm2))
            end if
        end associate

            r = r + 1
            if (end_loop) return
        end do increase_rank
    end subroutine smatcross_aca

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dmatcross_aca &
    (m, n, matval, col_order, maxr, niter_rook, r, U, ldu, VT, ldvt, irow, icol, work, lwork, info, rtol)
    use maria_kinds_mod,   only: &
        WP => DP
    use maria_constants_mod,  only: &
        ZERO => D_ZERO,             &
        ONE => D_ONE
    use maria_comparison_mod, only: &
        safe_less,                  &
        safe_eq
    use maria_access_matrix_mod, only: &
        MV => dmatval
    use maria_la_core_mod,    only: &
        dgenrmc,                    &
        dnrm2,                      &
        ddot,                       &
        dgemv,                      &
        droundup_lwork
    use maria_argcheck_mod,   only: &
        arg_is_bad,                 &
        BAD_IF_LESS,                &
        BAD_IF_MORE
    use maria_lr_la_mod,      only: &
        dlrval
    use maria_reports_mod,    only: &
        report_bad_arg,             &
        report_runtime_err,         &
        IND_OUT_OF_RANGE_ERR_CODE,  &
        SINGULAR_MATRIX_ERR_CODE
    !-- Input/output arguments -------------------------------------------------
        integer,       intent(in)                :: m
        integer,       intent(in)                :: n
        procedure(MV), intent(in),    pointer    :: matval
        integer,       intent(in),    contiguous :: col_order(:)
        integer,       intent(in)                :: maxr
        integer,       intent(in)                :: niter_rook
        integer,       intent(out)               :: r
        real(WP),      intent(out),   contiguous :: U(:)
        integer,       intent(in)                :: ldu
        real(WP),      intent(out),   contiguous :: VT(:)
        integer,       intent(in)                :: ldvt
        integer,       intent(out),   contiguous :: irow(:)
        integer,       intent(out),   contiguous :: icol(:)
        real(WP),      intent(out),   contiguous :: work(:)
        integer,       intent(in)                :: lwork
        integer,       intent(out)               :: info(2)
        real(WP),      intent(in),    optional   :: rtol

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = 'DMATCROSS_ACA'
        logical                 :: lw_query, end_loop
        integer                 :: req_lw, cur_icol, cur_irow, cur_icol_ind, i, j, &
                                   i_uu,  i_vv, pos(2)
        real(WP)                :: total_norm2, new_norm, appr_val, pivot, &
                                   maxv, nrmu, nrmv

    !-- Default values ---------------------------------------------------------
        r = 0

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, m, 0)) then
            info(1) = -1
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info(1) = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, maxr, 0)) then
            info(1) = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, maxr, min(m,n))) then
            info(1) = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, niter_rook, 0)) then
            info(1) = -6
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldU, max(1,m))) then
            info(1) = -9
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ldVT, max(1,maxr))) then
            info(1) = -11
            exit sanity
        end if
        if (present(rtol)) then
            if (arg_is_bad(BAD_IF_LESS, rtol, ZERO)) then
                info(1) = -17
                exit sanity
            end if
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info(1) /= 0) exit quickr

        if (m == 0 .or. n == 0 .or. maxr == 0) return
    end block quickr

    !-- Compute required workspace, return/check -------------------------------
    wrkwrk: block
        if (info(1) /= 0) exit wrkwrk
        
        if (present(rtol)) then
            req_lw = 2 * maxr
        else
            req_lw = 0
        end if

        lw_query = (lwork == -1)
        if (lw_query) then
            work(1) = droundup_lwork(req_lw)
            return
        end if

        if (arg_is_bad(BAD_IF_LESS, lwork, req_lw)) then
            info(1) = -15
            exit wrkwrk
        end if
    end block wrkwrk

    !-- Report bad input -------------------------------------------------------
        if (info(1) /= 0) then
            call report_bad_arg(SRNAME, -info(1))
            return
        end if

    !-- Executable section -----------------------------------------------------
        cur_icol_ind = 0
        total_norm2 = ZERO

        if (present(rtol)) then
            i_uu = 1
            i_vv = 1 + maxr
        end if

        increase_rank: do
            next_col: do
                cur_icol_ind = mod(cur_icol_ind, n) + 1
                if (all(icol(1:r) /= col_order(cur_icol_ind))) exit next_col
            end do next_col
            cur_icol = col_order(cur_icol_ind)
            if (cur_icol < 1 .or. cur_icol > n) then
                info = 1
                call report_runtime_err(SRNAME, IND_OUT_OF_RANGE_ERR_CODE)
                return
            end if

        associate &
        (col => U(1 + r*ldu:), row => VT(1 + r:))
            better_col: do i = 1, niter_rook
                do j = 1, m
                    appr_val = dlrval(m, n, j, cur_icol, r, U, ldu, VT, ldvt, info(1))
                    col(j) = matval(m, n, j, cur_icol, info(1)) - appr_val
                end do
                maxv = dgenrmc(m, 1, col, m, info(1), pos)
                cur_irow = pos(1)

                do j = 1, n
                    appr_val = dlrval(m, n, cur_irow, j, r, U, ldu, VT, ldvt, info(1))
                    row(1 + (j-1)*ldvt) = matval(m, n, cur_irow, j, info(1)) - appr_val
                end do
                maxv = dgenrmc(1, n, row, ldvt, info(1), pos)
                if (pos(2) == cur_icol) then
                    exit better_col
                else
                    cur_icol = pos(2)
                    info(2) = info(2) + 1
                end if
            end do better_col

            find_pivot: do j = 1, m
                appr_val = dlrval(m, n, j, cur_icol, r, U, ldu, VT, ldvt, info(1))
                col(j) = matval(m, n, j, cur_icol, info(1)) - appr_val
            end do find_pivot
            maxv = dgenrmc(m, 1, col, m, info(1), pos)
            cur_irow = pos(1)

        !-- Rank-1 update
            appr_val = dlrval(m, n, cur_irow, cur_icol, r, U, ldu, VT, ldvt, info(1))
            pivot = matval(m, n, cur_irow, cur_icol, info(1)) - appr_val
            if (safe_eq(pivot, ZERO)) then
                info = 2
                call report_runtime_err(SRNAME, SINGULAR_MATRIX_ERR_CODE)
                return
            end if
            do j = 1, n
                appr_val = dlrval(m, n, cur_irow, j, r, U, ldu, VT, ldvt, info(1))
                row(1 + (j-1)*ldvt) = matval(m, n, cur_irow, j, info(1)) - appr_val
                row(1 + (j-1)*ldvt) = row(1 + (j-1)*ldvt) / pivot
            end do

            irow(r+1) = cur_irow
            icol(r+1) = cur_icol

        !-- Stopping criteria
            end_loop = (r+1 == maxr)
            if (present(rtol)) then
                nrmu = dnrm2(m, col, 1)
                nrmv = dnrm2(n, row, ldvt)
                new_norm = nrmu * nrmv
                total_norm2 = total_norm2 + new_norm**2
                if (r > 0) then
                associate &
                (uu => work(i_uu:), vv => work(i_vv:))
                    call dgemv('t', m, r, ONE, U, ldu, col, 1, ZERO, uu, 1)
                    call dgemv('n', r, n, ONE, VT, ldvt, row, ldvt, ZERO, vv, 1)
                    total_norm2 = total_norm2 + 2 * ddot(r, uu, 1, vv, 1)
                end associate
                end if
                end_loop = end_loop .or. safe_less(new_norm, rtol * sqrt(total_norm2))
            end if
        end associate

            r = r + 1
            if (end_loop) return
        end do increase_rank
    end subroutine dmatcross_aca

    !------------------------------------------------------------------------------------------------------------------------
end submodule maria_lr_cross_sub
