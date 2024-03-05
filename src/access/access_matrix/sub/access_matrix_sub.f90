!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the implementation of the [[maria_access_matrix_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Implements the [[maria_access_matrix_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------

submodule (maria_access_matrix_mod) maria_access_matrix_sub
implicit none (type, external)

contains
    !------------------------------------------------------------------------------------------------------------------------

    module subroutine smatval2slc &
    (matval, m, n, mode, ind, x, incx, info)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS,               &
        BAD_IF_MORE
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        procedure(smatval), intent(in),  pointer    :: matval
        integer,            intent(in)              :: m
        integer,            intent(in)              :: n
        integer,            intent(in)              :: mode
        integer,            intent(in)              :: ind
        real(WP),           intent(out), contiguous :: x(:)
        integer,            intent(in)              :: incx
        integer,            intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "SMATVAL2SLC"
        integer                 :: i, posx

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, m, 1)) then
            info = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, n, 1)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, mode, 1)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, mode, 2)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ind, 1)) then
            info = -5
            exit sanity
        end if
        if (mode == 1) then
            if (arg_is_bad(BAD_IF_MORE, ind, m)) then
                info = -5
                exit sanity
            end if
        else
            if (arg_is_bad(BAD_IF_MORE, ind, n)) then
                info = -5
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
        if (mode == 1) then
            if (incx >= 0) then
                posx = 1
             else
                posx = 1 - (n-1) * incx
            end if
            do i = 1, n
                x(posx) = matval(m, n, ind, i, info)
                if (info /= 0) return
                posx = posx + incx
            end do
        else
            if (incx >= 0) then
                posx = 1
             else
                posx = 1 - (m-1) * incx
            end if
            do i = 1, m
                x(posx) = matval(m, n, i, ind, info)
                if (info /= 0) return
                posx = posx + incx
            end do
        end if
    end subroutine smatval2slc

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dmatval2slc &
    (matval, m, n, mode, ind, x, incx, info)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS,               &
        BAD_IF_MORE
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        procedure(dmatval), intent(in),  pointer    :: matval
        integer,            intent(in)              :: m
        integer,            intent(in)              :: n
        integer,            intent(in)              :: mode
        integer,            intent(in)              :: ind
        real(WP),           intent(out), contiguous :: x(:)
        integer,            intent(in)              :: incx
        integer,            intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "DMATVAL2SLC"
        integer                 :: i, posx

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, m, 1)) then
            info = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, n, 1)) then
            info = -3
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, mode, 1)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_MORE, mode, 2)) then
            info = -4
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, ind, 1)) then
            info = -5
            exit sanity
        end if
        if (mode == 1) then
            if (arg_is_bad(BAD_IF_MORE, ind, m)) then
                info = -5
                exit sanity
            end if
        else
            if (arg_is_bad(BAD_IF_MORE, ind, n)) then
                info = -5
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
        if (mode == 1) then
            if (incx >= 0) then
                posx = 1
             else
                posx = 1 - (n-1) * incx
            end if
            do i = 1, n
                x(posx) = matval(m, n, ind, i, info)
                if (info /= 0) return
                posx = posx + incx
            end do
        else
            if (incx >= 0) then
                posx = 1
             else
                posx = 1 - (m-1) * incx
            end if
            do i = 1, m
                x(posx) = matval(m, n, i, ind, info)
                if (info /= 0) return
                posx = posx + incx
            end do
        end if
    end subroutine dmatval2slc

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine smatval2full &
    (matval, trans, m, n, A, lda, info)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        character(1),       intent(in)              :: trans
        procedure(smatval), intent(in),  pointer    :: matval
        integer,            intent(in)              :: m
        integer,            intent(in)              :: n
        real(WP),           intent(out), contiguous :: A(:)
        integer,            intent(in)              :: lda
        integer,            intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "SMATVAL2FULL"
        logical                 :: transposed
        integer                 :: i, j

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        select case (trans)
            case ('n', 'N')
                transposed = .false.
            case ('t', 'T')
                transposed = .true.
            case default
                transposed = .false.
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
        if (transposed) then
            if (arg_is_bad(BAD_IF_LESS, lda, max(1,n))) then
                info = -6
                exit sanity
            end if
        else
            if (arg_is_bad(BAD_IF_LESS, lda, max(1,m))) then
                info = -6
                exit sanity
            end if
        endif
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
        if (transposed) then
            do i = 1, m
                do j = 1, n
                    A(j + (i-1)*lda) = matval(m, n, i, j, info)
                end do
            end do
        else
            do j = 1, n
                do i = 1, m
                    A(i + (j-1)*lda) = matval(m, n, i, j, info)
                end do
            end do
        end if
    end subroutine smatval2full

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dmatval2full &
    (matval, trans, m, n, A, lda, info)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        character(1),       intent(in)              :: trans
        procedure(dmatval), intent(in),  pointer    :: matval
        integer,            intent(in)              :: m
        integer,            intent(in)              :: n
        real(WP),           intent(out), contiguous :: A(:)
        integer,            intent(in)              :: lda
        integer,            intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "DMATVAL2FULL"
        logical                 :: transposed
        integer                 :: i, j

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        select case (trans)
            case ('n', 'N')
                transposed = .false.
            case ('t', 'T')
                transposed = .true.
            case default
                transposed = .false.
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
        if (transposed) then
            if (arg_is_bad(BAD_IF_LESS, lda, max(1,n))) then
                info = -6
                exit sanity
            end if
        else
            if (arg_is_bad(BAD_IF_LESS, lda, max(1,m))) then
                info = -6
                exit sanity
            end if
        endif
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
        if (transposed) then
            do i = 1, m
                do j = 1, n
                    A(j + (i-1)*lda) = matval(m, n, i, j, info)
                end do
            end do
        else
            do j = 1, n
                do i = 1, m
                    A(i + (j-1)*lda) = matval(m, n, i, j, info)
                end do
            end do
        end if
    end subroutine dmatval2full

    !------------------------------------------------------------------------------------------------------------------------

    module function smatval_hilbert &
    (m, n, i, j, info)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        ONE => S_ONE
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS,               &
        BAD_IF_MORE
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer, intent(in)  :: m
        integer, intent(in)  :: n
        integer, intent(in)  :: i
        integer, intent(in)  :: j
        integer, intent(out) :: info
        real(WP)             :: smatval_hilbert

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "SMATVAL_HILBERT"

    !-- Default value ----------------------------------------------------------
    smatval_hilbert = ONE

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
    end block sanity

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        smatval_hilbert = ONE / (i + j - ONE)
    end function smatval_hilbert

    !------------------------------------------------------------------------------------------------------------------------

    module function dmatval_hilbert &
    (m, n, i, j, info)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        ONE => D_ONE
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS,               &
        BAD_IF_MORE
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer, intent(in)  :: m
        integer, intent(in)  :: n
        integer, intent(in)  :: i
        integer, intent(in)  :: j
        integer, intent(out) :: info
        real(WP)             :: dmatval_hilbert

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "DMATVAL_HILBERT"

    !-- Default value ----------------------------------------------------------
    dmatval_hilbert = ONE

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
    end block sanity

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        dmatval_hilbert = ONE / (i + j - ONE)
    end function dmatval_hilbert

    !------------------------------------------------------------------------------------------------------------------------

    module function smatval_lotkin &
    (m, n, i, j, info)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        ONE => S_ONE
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS,               &
        BAD_IF_MORE
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer, intent(in)  :: m
        integer, intent(in)  :: n
        integer, intent(in)  :: i
        integer, intent(in)  :: j
        integer, intent(out) :: info
        real(WP)             :: smatval_lotkin

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "SMATVAL_LOTKIN"

    !-- Default value ----------------------------------------------------------
    smatval_lotkin = ONE

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
    end block sanity

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (i == 1) then
            smatval_lotkin = ONE
        else
            smatval_lotkin = ONE / (i + j - ONE)
        end if
    end function smatval_lotkin

    !------------------------------------------------------------------------------------------------------------------------

    module function dmatval_lotkin &
    (m, n, i, j, info)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        ONE => D_ONE
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS,               &
        BAD_IF_MORE
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer, intent(in)  :: m
        integer, intent(in)  :: n
        integer, intent(in)  :: i
        integer, intent(in)  :: j
        integer, intent(out) :: info
        real(WP)             :: dmatval_lotkin

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "DMATVAL_LOTKIN"

    !-- Default value ----------------------------------------------------------
    dmatval_lotkin = ONE

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
    end block sanity

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        if (i == 1) then
            dmatval_lotkin = ONE
        else
            dmatval_lotkin = ONE / (i + j - ONE)
        end if
    end function dmatval_lotkin

    !------------------------------------------------------------------------------------------------------------------------

    module function smatval_ballistic &
    (m, n, i, j, info)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        ONE => S_ONE
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS,               &
        BAD_IF_MORE
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer, intent(in)  :: m
        integer, intent(in)  :: n
        integer, intent(in)  :: i
        integer, intent(in)  :: j
        integer, intent(out) :: info
        real(WP)             :: smatval_ballistic

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "SMATVAL_BALLISTIC"

    !-- Default value ----------------------------------------------------------
    smatval_ballistic = ONE

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
    end block sanity

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        smatval_ballistic = sqrt(ONE / i + ONE / j) * (i**(ONE/3) + j**(ONE/3))
    end function smatval_ballistic

    !------------------------------------------------------------------------------------------------------------------------

    module function dmatval_ballistic &
    (m, n, i, j, info)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        ONE => D_ONE
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS,               &
        BAD_IF_MORE
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        integer, intent(in)  :: m
        integer, intent(in)  :: n
        integer, intent(in)  :: i
        integer, intent(in)  :: j
        integer, intent(out) :: info
        real(WP)             :: dmatval_ballistic

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "DMATVAL_BALLISTIC"

    !-- Default value ----------------------------------------------------------
    dmatval_ballistic = ONE

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
    end block sanity

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        dmatval_ballistic = sqrt(ONE / i + ONE / j) * (i**(ONE/3) + j**(ONE/3))
    end function dmatval_ballistic

    !------------------------------------------------------------------------------------------------------------------------
end submodule maria_access_matrix_sub
