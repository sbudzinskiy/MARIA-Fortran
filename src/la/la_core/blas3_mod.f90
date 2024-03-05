!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_blas3_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Provides explicit interfaces to Level 3 BLAS subroutines.
!----------------------------------------------------------------------------------------------------------------------------
module maria_blas3_mod
use maria_kinds_mod, only: &
    SP, DP
implicit none (type, external)

! Procedures
public :: sgemm, &
          dgemm, &
          strmm, &
          dtrmm, &
          strsm, &
          dtrsm

interface
    !> Adds a matrix-matrix product to a matrix
    subroutine sgemm &
    (transA, transB, m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        character(1), intent(in)    :: transA
        character(1), intent(in)    :: transB
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        integer,      intent(in)    :: k
        real(WP),     intent(in)    :: alpha
        real(WP),     intent(in)    :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(in)    :: B(*)
        integer,      intent(in)    :: ldB
        real(WP),     intent(in)    :: beta
        real(WP),     intent(inout) :: C(*)
        integer,      intent(in)    :: ldC
    end subroutine sgemm

    !> Adds a matrix-matrix product to a matrix
    subroutine dgemm &
    (transA, transB, m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        character(1), intent(in)    :: transA
        character(1), intent(in)    :: transB
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        integer,      intent(in)    :: k
        real(WP),     intent(in)    :: alpha
        real(WP),     intent(in)    :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(in)    :: B(*)
        integer,      intent(in)    :: ldB
        real(WP),     intent(in)    :: beta
        real(WP),     intent(inout) :: C(*)
        integer,      intent(in)    :: ldC
    end subroutine dgemm

    !> Multiplies in-place a matrix by a triangular matrix 
    subroutine strmm &
    (side, uplo, trans, diag, m, n, alpha, A, ldA, B, ldB)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        character(1), intent(in)    :: side
        character(1), intent(in)    :: uplo
        character(1), intent(in)    :: trans
        character(1), intent(in)    :: diag
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        real(WP),     intent(in)    :: alpha
        real(WP),     intent(in)    :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(inout) :: B(*)
        integer,      intent(in)    :: ldB
    end subroutine strmm

    !> Multiplies in-place a matrix by a triangular matrix 
    subroutine dtrmm &
    (side, uplo, trans, diag, m, n, alpha, A, ldA, B, ldB)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        character(1), intent(in)    :: side
        character(1), intent(in)    :: uplo
        character(1), intent(in)    :: trans
        character(1), intent(in)    :: diag
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        real(WP),     intent(in)    :: alpha
        real(WP),     intent(in)    :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(inout) :: B(*)
        integer,      intent(in)    :: ldB
    end subroutine dtrmm

    !> Solves in-place a linear system with a trinagular matrix and multiple rhs
    subroutine strsm &
    (side, uplo, trans, diag, m, n, alpha, A, ldA, B, ldB)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        character(1), intent(in)    :: side
        character(1), intent(in)    :: uplo
        character(1), intent(in)    :: trans
        character(1), intent(in)    :: diag
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        real(WP),     intent(in)    :: alpha
        real(WP),     intent(in)    :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(inout) :: B(*)
        integer,      intent(in)    :: ldB
    end subroutine strsm

    !> Solves in-place a linear system with a trinagular matrix and multiple rhs
    subroutine dtrsm &
    (side, uplo, trans, diag, m, n, alpha, A, ldA, B, ldB)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        character(1), intent(in)    :: side
        character(1), intent(in)    :: uplo
        character(1), intent(in)    :: trans
        character(1), intent(in)    :: diag
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        real(WP),     intent(in)    :: alpha
        real(WP),     intent(in)    :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(inout) :: B(*)
        integer,      intent(in)    :: ldB
    end subroutine dtrsm
end interface
end module maria_blas3_mod
