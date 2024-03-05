!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_blas1_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Provides explicit interfaces to Level 1 BLAS subroutines.
!----------------------------------------------------------------------------------------------------------------------------
module maria_blas1_mod
implicit none (type, external)

! Procedures
public :: srotg, &
          drotg, &
          srot,  &
          drot,  &
          sswap, &
          dswap, &
          sscal, &
          dscal, &
          scopy, &
          dcopy, &
          saxpy, &
          daxpy, &
          sdot,  &
          ddot,  &
          snrm2, &
          dnrm2

interface
    !> Generates a plane rotation 
    subroutine srotg &
    (a, b, c, s)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        real(WP), intent(inout) :: a
        real(WP), intent(inout) :: b
        real(WP), intent(out)   :: c
        real(WP), intent(out)   :: s
    end subroutine srotg

    !> Generates a plane rotation 
    subroutine drotg &
    (a, b, c, s)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        real(WP), intent(inout) :: a
        real(WP), intent(inout) :: b
        real(WP), intent(out)   :: c
        real(WP), intent(out)   :: s
    end subroutine drotg

    !> Applies a plane rotation 
    subroutine srot &
    (n, x, incx, y, incy, c, s)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,  intent(in)    :: n
        real(WP), intent(inout) :: x(*)
        integer,  intent(in)    :: incx
        real(WP), intent(inout) :: y(*)
        integer,  intent(in)    :: incy
        real(WP), intent(in)    :: c
        real(WP), intent(in)    :: s
    end subroutine srot

    !> Applies a plane rotation 
    subroutine drot &
    (n, x, incx, y, incy, c, s)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,  intent(in)    :: n
        real(WP), intent(inout) :: x(*)
        integer,  intent(in)    :: incx
        real(WP), intent(inout) :: y(*)
        integer,  intent(in)    :: incy
        real(WP), intent(in)    :: c
        real(WP), intent(in)    :: s
    end subroutine drot

    !> Swaps the entries of two vectors
    subroutine sswap &
    (n, x, incx, y, incy)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,  intent(in)    :: n
        real(WP), intent(inout) :: x(*)
        integer,  intent(in)    :: incx
        real(WP), intent(inout) :: y(*)
        integer,  intent(in)    :: incy
    end subroutine sswap

    !> Swaps the entries of two vectors
    subroutine dswap &
    (n, x, incx, y, incy)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,  intent(in)    :: n
        real(WP), intent(inout) :: x(*)
        integer,  intent(in)    :: incx
        real(WP), intent(inout) :: y(*)
        integer,  intent(in)    :: incy
    end subroutine dswap

    !> Scales a vector by a constant
    subroutine sscal &
    (n, alpha, x, incx)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,  intent(in)    :: n
        real(WP), intent(in)    :: alpha
        real(WP), intent(inout) :: x(*)
        integer,  intent(in)    :: incx
    end subroutine sscal

    !> Scales a vector by a constant
    subroutine dscal &
    (n, alpha, x, incx)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,  intent(in)    :: n
        real(WP), intent(in)    :: alpha
        real(WP), intent(inout) :: x(*)
        integer,  intent(in)    :: incx
    end subroutine dscal

    !> Copies a vector to another vector
    subroutine scopy &
    (n, x, incx, y, incy)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,  intent(in)  :: n
        real(WP), intent(in)  :: x(*)
        integer,  intent(in)  :: incx
        real(WP), intent(out) :: y(*)
        integer,  intent(in)  :: incy
    end subroutine scopy

    !> Copies a vector to another vector
    subroutine dcopy &
    (n, x, incx, y, incy)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,  intent(in)  :: n
        real(WP), intent(in)  :: x(*)
        integer,  intent(in)  :: incx
        real(WP), intent(out) :: y(*)
        integer,  intent(in)  :: incy
    end subroutine dcopy

    !> Adds a constant times a vector to a vector
    subroutine saxpy &
    (n, alpha, x, incx, y, incy)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,  intent(in)    :: n
        real(WP), intent(in)    :: alpha
        real(WP), intent(in)    :: x(*)
        integer,  intent(in)    :: incx
        real(WP), intent(inout) :: y(*)
        integer,  intent(in)    :: incy
    end subroutine saxpy

    !> Adds a constant times a vector to a vector
    subroutine daxpy &
    (n, alpha, x, incx, y, incy)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,  intent(in)    :: n
        real(WP), intent(in)    :: alpha
        real(WP), intent(in)    :: x(*)
        integer,  intent(in)    :: incx
        real(WP), intent(inout) :: y(*)
        integer,  intent(in)    :: incy
    end subroutine daxpy

    !> Computes the dot product of two vectors
    function sdot &
    (n, x, incx, y, incy)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,  intent(in) :: n
        real(WP), intent(in) :: x(*)
        integer,  intent(in) :: incx
        real(WP), intent(in) :: y(*)
        integer,  intent(in) :: incy
        real(WP)             :: sdot
    end function sdot

    !> Computes the dot product of two vectors
    function ddot &
    (n, x, incx, y, incy)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,  intent(in) :: n
        real(WP), intent(in) :: x(*)
        integer,  intent(in) :: incx
        real(WP), intent(in) :: y(*)
        integer,  intent(in) :: incy
        real(WP)             :: ddot
    end function ddot

    !> Computes the Euclidean norm of a vector
    function snrm2 &
    (n, x, incx)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,  intent(in) :: n
        real(WP), intent(in) :: x(*)
        integer,  intent(in) :: incx
        real(WP)             :: snrm2
    end function snrm2

    !> Computes the Euclidean norm of a vector
    function dnrm2 &
    (n, x, incx)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,  intent(in) :: n
        real(WP), intent(in) :: x(*)
        integer,  intent(in) :: incx
        real(WP)             :: dnrm2
    end function dnrm2
end interface
end module maria_blas1_mod
