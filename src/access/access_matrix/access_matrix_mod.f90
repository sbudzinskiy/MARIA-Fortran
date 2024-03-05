!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_access_matrix_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Interfaces for accessing matrices via their elements, columns, and rows.
!----------------------------------------------------------------------------------------------------------------------------
module maria_access_matrix_mod
implicit none (type, external)

! Interfaces
public :: smatval, &
          dmatval, &
          smatslc, &
          dmatslc

! Procedures
public :: smatval2slc,       &
          dmatval2slc,       &
          smatval2full,      &
          dmatval2full,      &
          smatval_hilbert,   &
          dmatval_hilbert,   &
          smatval_lotkin,    &
          dmatval_lotkin,    &
          smatval_ballistic, &
          dmatval_ballistic

abstract interface
    !------------------------------------------------------------------------------------------------------------------------

    function smatval &
    (m, n, i, j, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer, intent(in)  :: m
        integer, intent(in)  :: n
        integer, intent(in)  :: i
        integer, intent(in)  :: j
        integer, intent(out) :: info
        real(WP)             :: smatval
    end function smatval

    !------------------------------------------------------------------------------------------------------------------------

    function dmatval &
    (m, n, i, j, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer, intent(in)  :: m
        integer, intent(in)  :: n
        integer, intent(in)  :: i
        integer, intent(in)  :: j
        integer, intent(out) :: info
        real(WP)             :: dmatval
    end function dmatval

    !------------------------------------------------------------------------------------------------------------------------

    subroutine smatslc &
    (m, n, mode, ind, x, incx, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,  intent(in)              :: m
        integer,  intent(in)              :: n
        integer,  intent(in)              :: mode
        integer,  intent(in)              :: ind
        real(WP), intent(out), contiguous :: x(:)
        integer,  intent(in)              :: incx
        integer,  intent(out)             :: info
    end subroutine smatslc

    !------------------------------------------------------------------------------------------------------------------------

    subroutine dmatslc &
    (m, n, mode, ind, x, incx, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,  intent(in)              :: m
        integer,  intent(in)              :: n
        integer,  intent(in)              :: mode
        integer,  intent(in)              :: ind
        real(WP), intent(out), contiguous :: x(:)
        integer,  intent(in)              :: incx
        integer,  intent(out)             :: info
    end subroutine dmatslc

    !------------------------------------------------------------------------------------------------------------------------
end interface

interface
    !------------------------------------------------------------------------------------------------------------------------

    !> Computes columns and rows of a matrix element by element.
    module subroutine smatval2slc &
    (matval, m, n, mode, ind, x, incx, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        procedure(smatval), intent(in),  pointer    :: matval
        integer,            intent(in)              :: m
        integer,            intent(in)              :: n
        integer,            intent(in)              :: mode
        integer,            intent(in)              :: ind
        real(WP),           intent(out), contiguous :: x(:)
        integer,            intent(in)              :: incx
        integer,            intent(out)             :: info
    end subroutine smatval2slc

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes columns and rows of a matrix element by element.
    module subroutine dmatval2slc &
    (matval, m, n, mode, ind, x, incx, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        procedure(dmatval), intent(in),  pointer    :: matval
        integer,            intent(in)              :: m
        integer,            intent(in)              :: n
        integer,            intent(in)              :: mode
        integer,            intent(in)              :: ind
        real(WP),           intent(out), contiguous :: x(:)
        integer,            intent(in)              :: incx
        integer,            intent(out)             :: info
    end subroutine dmatval2slc

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes a matrix element by element.
    module subroutine smatval2full &
    (matval, trans, m, n, A, lda, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        character(1),       intent(in)              :: trans
        procedure(smatval), intent(in),  pointer    :: matval
        integer,            intent(in)              :: m
        integer,            intent(in)              :: n
        real(WP),           intent(out), contiguous :: A(:)
        integer,            intent(in)              :: lda
        integer,            intent(out)             :: info
    end subroutine smatval2full

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes a matrix element by element.
    module subroutine dmatval2full &
    (matval, trans, m, n, A, lda, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        character(1),       intent(in)              :: trans
        procedure(dmatval), intent(in),  pointer    :: matval
        integer,            intent(in)              :: m
        integer,            intent(in)              :: n
        real(WP),           intent(out), contiguous :: A(:)
        integer,            intent(in)              :: lda
        integer,            intent(out)             :: info
    end subroutine dmatval2full

    !------------------------------------------------------------------------------------------------------------------------

    module function smatval_hilbert &
    (m, n, i, j, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer, intent(in)  :: m
        integer, intent(in)  :: n
        integer, intent(in)  :: i
        integer, intent(in)  :: j
        integer, intent(out) :: info
        real(WP)             :: smatval_hilbert
    end function smatval_hilbert

    !------------------------------------------------------------------------------------------------------------------------

    module function dmatval_hilbert &
    (m, n, i, j, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer, intent(in)  :: m
        integer, intent(in)  :: n
        integer, intent(in)  :: i
        integer, intent(in)  :: j
        integer, intent(out) :: info
        real(WP)             :: dmatval_hilbert
    end function dmatval_hilbert

    !------------------------------------------------------------------------------------------------------------------------

    module function smatval_lotkin &
    (m, n, i, j, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer, intent(in)  :: m
        integer, intent(in)  :: n
        integer, intent(in)  :: i
        integer, intent(in)  :: j
        integer, intent(out) :: info
        real(WP)             :: smatval_lotkin
    end function smatval_lotkin

    !------------------------------------------------------------------------------------------------------------------------

    module function dmatval_lotkin &
    (m, n, i, j, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer, intent(in)  :: m
        integer, intent(in)  :: n
        integer, intent(in)  :: i
        integer, intent(in)  :: j
        integer, intent(out) :: info
        real(WP)             :: dmatval_lotkin
    end function dmatval_lotkin

    !------------------------------------------------------------------------------------------------------------------------

    module function smatval_ballistic &
    (m, n, i, j, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer, intent(in)  :: m
        integer, intent(in)  :: n
        integer, intent(in)  :: i
        integer, intent(in)  :: j
        integer, intent(out) :: info
        real(WP)             :: smatval_ballistic
    end function smatval_ballistic

    !------------------------------------------------------------------------------------------------------------------------

    module function dmatval_ballistic &
    (m, n, i, j, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer, intent(in)  :: m
        integer, intent(in)  :: n
        integer, intent(in)  :: i
        integer, intent(in)  :: j
        integer, intent(out) :: info
        real(WP)             :: dmatval_ballistic
    end function dmatval_ballistic

    !------------------------------------------------------------------------------------------------------------------------
end interface

end module maria_access_matrix_mod
