!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_blas2_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Provides explicit interfaces to Level 2 BLAS subroutines.
!----------------------------------------------------------------------------------------------------------------------------
module maria_blas2_mod
use maria_kinds_mod, only: &
    SP, DP
implicit none (type, external)

! Procedures
public :: sgemv, &
          dgemv, &
          strmv, &
          dtrmv, &
          strsv, &
          dtrsv, &
          sger, &
          dger

interface
    !> Adds a matrix-vector product to a vector 
    subroutine sgemv &
    (trans, m, n, alpha, A, ldA, x, incx, beta, y, incy)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        character(1), intent(in)    :: trans
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        real(WP),     intent(in)    :: alpha
        real(WP),     intent(in)    :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(in)    :: x(*)
        integer,      intent(in)    :: incx
        real(WP),     intent(in)    :: beta
        real(WP),     intent(inout) :: y(*)
        integer,      intent(in)    :: incy
    end subroutine sgemv

    !> Adds a matrix-vector product to a vector 
    subroutine dgemv &
    (trans, m, n, alpha, A, ldA, x, incx, beta, y, incy)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        character(1), intent(in)    :: trans
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        real(WP),     intent(in)    :: alpha
        real(WP),     intent(in)    :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(in)    :: x(*)
        integer,      intent(in)    :: incx
        real(WP),     intent(in)    :: beta
        real(WP),     intent(inout) :: y(*)
        integer,      intent(in)    :: incy
    end subroutine dgemv

    !> Multiplies in-place a vector by a triangular matrix 
    subroutine strmv &
    (uplo, trans, diag, n, A, ldA, x, incx)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        character(1), intent(in)    :: uplo
        character(1), intent(in)    :: trans
        character(1), intent(in)    :: diag
        integer,      intent(in)    :: n
        real(WP),     intent(in)    :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(inout) :: x(*)
        integer,      intent(in)    :: incx
    end subroutine strmv

    !> Multiplies in-place a vector by a triangular matrix 
    subroutine dtrmv &
    (uplo, trans, diag, n, A, ldA, x, incx)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        character(1), intent(in)    :: uplo
        character(1), intent(in)    :: trans
        character(1), intent(in)    :: diag
        integer,      intent(in)    :: n
        real(WP),     intent(in)    :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(inout) :: x(*)
        integer,      intent(in)    :: incx
    end subroutine dtrmv

    !> Solves in-place a linear system with a trinagular matrix
    subroutine strsv &
    (uplo, trans, diag, n, A, ldA, x, incx)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        character(1), intent(in)    :: uplo
        character(1), intent(in)    :: trans
        character(1), intent(in)    :: diag
        integer,      intent(in)    :: n
        real(WP),     intent(in)    :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(inout) :: x(*)
        integer,      intent(in)    :: incx
    end subroutine strsv

    !> Solves in-place a linear system with a trinagular matrix
    subroutine dtrsv &
    (uplo, trans, diag, n, A, ldA, x, incx)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        character(1), intent(in)    :: uplo
        character(1), intent(in)    :: trans
        character(1), intent(in)    :: diag
        integer,      intent(in)    :: n
        real(WP),     intent(in)    :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(inout) :: x(*)
        integer,      intent(in)    :: incx
    end subroutine dtrsv

    !> Adds a rank-1 matrix to a given matrix
    subroutine sger &
    (m, n, alpha, x, incx, y, incy, A, ldA)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        real(WP),     intent(in)    :: alpha
        real(WP),     intent(in)    :: x(*)
        integer,      intent(in)    :: incx
        real(WP),     intent(in)    :: y(*)
        integer,      intent(in)    :: incy
        real(WP),     intent(inout) :: A(*)
        integer,      intent(in)    :: ldA
    end subroutine sger

    !> Adds a rank-1 matrix to a given matrix
    subroutine dger &
    (m, n, alpha, x, incx, y, incy, A, ldA)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        real(WP),     intent(in)    :: alpha
        real(WP),     intent(in)    :: x(*)
        integer,      intent(in)    :: incx
        real(WP),     intent(in)    :: y(*)
        integer,      intent(in)    :: incy
        real(WP),     intent(inout) :: A(*)
        integer,      intent(in)    :: ldA
    end subroutine dger
end interface
end module maria_blas2_mod
