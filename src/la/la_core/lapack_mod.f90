!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_lapack_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Provides explicit interfaces to certain LAPACK subroutines.
!----------------------------------------------------------------------------------------------------------------------------
module maria_lapack_mod
use maria_kinds_mod, only: &
    SP, DP
implicit none (type, external)

! Procedures
public :: sgetrf, &
          dgetrf, &
          sgeqrf, &
          dgeqrf, &
          sgelqf, &
          dgelqf, &
          sgesdd, &
          dgesdd, &
          strtri, &
          dtrtri, &
          sormqr, &
          dormqr, &
          sorgqr, &
          dorgqr, &
          sormlq, &
          dormlq, &
          sorglq, &
          dorglq, &
          slarf,  &
          dlarf,  &
          slarfg, &
          dlarfg, &
          slacpy, &
          dlacpy, &
          slaset, &
          dlaset, &
          slasrt, &
          dlasrt

interface
    !> Computes an LU decomposition with partial pivoting
    subroutine sgetrf &
    (m, n, A, ldA, ipiv, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,  intent(in)    :: m
        integer,  intent(in)    :: n
        real(WP), intent(inout) :: A(*)
        integer,  intent(in)    :: ldA
        integer,  intent(out)   :: ipiv(*)
        integer,  intent(out)   :: info
    end subroutine sgetrf

    !> Computes an LU decomposition with partial pivoting
    subroutine dgetrf &
    (m, n, A, ldA, ipiv, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,  intent(in)    :: m
        integer,  intent(in)    :: n
        real(WP), intent(inout) :: A(*)
        integer,  intent(in)    :: ldA
        integer,  intent(out)   :: ipiv(*)
        integer,  intent(out)   :: info
    end subroutine dgetrf

    !> Computes a QR decomposition
    subroutine sgeqrf &
    (m, n, A, ldA, tau, work, lwork, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,  intent(in)    :: m
        integer,  intent(in)    :: n
        real(WP), intent(inout) :: A(*)
        integer,  intent(in)    :: ldA
        real(WP), intent(out)   :: tau(*)
        real(WP), intent(out)   :: work(*)
        integer,  intent(in)    :: lwork
        integer,  intent(out)   :: info
    end subroutine sgeqrf

    !> Computes a QR decomposition
    subroutine dgeqrf &
    (m, n, A, ldA, tau, work, lwork, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,  intent(in)    :: m
        integer,  intent(in)    :: n
        real(WP), intent(inout) :: A(*)
        integer,  intent(in)    :: ldA
        real(WP), intent(out)   :: tau(*)
        real(WP), intent(out)   :: work(*)
        integer,  intent(in)    :: lwork
        integer,  intent(out)   :: info
    end subroutine dgeqrf

    !> Computes an LQ decomposition
    subroutine sgelqf &
    (m, n, A, ldA, tau, work, lwork, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,  intent(in)    :: m
        integer,  intent(in)    :: n
        real(WP), intent(inout) :: A(*)
        integer,  intent(in)    :: ldA
        real(WP), intent(out)   :: tau(*)
        real(WP), intent(out)   :: work(*)
        integer,  intent(in)    :: lwork
        integer,  intent(out)   :: info
    end subroutine sgelqf

    !> Computes an LQ decomposition
    subroutine dgelqf &
    (m, n, A, ldA, tau, work, lwork, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,  intent(in)    :: m
        integer,  intent(in)    :: n
        real(WP), intent(inout) :: A(*)
        integer,  intent(in)    :: ldA
        real(WP), intent(out)   :: tau(*)
        real(WP), intent(out)   :: work(*)
        integer,  intent(in)    :: lwork
        integer,  intent(out)   :: info
    end subroutine dgelqf

    !> Computes a singular value decomposition
    subroutine sgesdd &
    (jobZ, m, n, A, ldA, S, U, ldU, VT, ldVT, work, lwork, iwork, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        character(1), intent(in)    :: jobZ
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        real(WP),     intent(inout) :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(out)   :: S(*)
        real(WP),     intent(out)   :: U(*)
        integer,      intent(in)    :: ldU
        real(WP),     intent(out)   :: VT(*)
        integer,      intent(in)    :: ldVT
        real(WP),     intent(out)   :: work(*)
        integer,      intent(in)    :: lwork
        integer,      intent(out)   :: iwork(*)
        integer,      intent(out)   :: info
    end subroutine sgesdd

    !> Computes a singular value decomposition
    subroutine dgesdd &
    (jobZ, m, n, A, ldA, S, U, ldU, VT, ldVT, work, lwork, iwork, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        character(1), intent(in)    :: jobZ
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        real(WP),     intent(inout) :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(out)   :: S(*)
        real(WP),     intent(out)   :: U(*)
        integer,      intent(in)    :: ldU
        real(WP),     intent(out)   :: VT(*)
        integer,      intent(in)    :: ldVT
        real(WP),     intent(out)   :: work(*)
        integer,      intent(in)    :: lwork
        integer,      intent(out)   :: iwork(*)
        integer,      intent(out)   :: info
    end subroutine dgesdd

    !> Computes the inverse of a triangular matrix
    subroutine strtri &
    (uplo, diag, n, A, ldA, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        character(1), intent(in)    :: uplo
        character(1), intent(in)    :: diag
        integer,      intent(in)    :: n
        real(WP),     intent(inout) :: A(*)
        integer,      intent(in)    :: ldA
        integer,      intent(out)   :: info
    end subroutine strtri

    !> Computes the inverse of a triangular matrix
    subroutine dtrtri &
    (uplo, diag, n, A, ldA, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        character(1), intent(in)    :: uplo
        character(1), intent(in)    :: diag
        integer,      intent(in)    :: n
        real(WP),     intent(inout) :: A(*)
        integer,      intent(in)    :: ldA
        integer,      intent(out)   :: info
    end subroutine dtrtri

    !> Multiplies in-place a matrix with an orthogonal matrix 
    !> defined as the product of elementary reflectors
    subroutine sormqr &
    (side, trans, m, n, k, A, ldA, tau, C, ldC, work, lwork, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        character(1), intent(in)    :: side
        character(1), intent(in)    :: trans
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        integer,      intent(in)    :: k
        real(WP),     intent(in)    :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(in)    :: tau(*)
        real(WP),     intent(inout) :: C(*)
        integer,      intent(in)    :: ldC
        real(WP),     intent(out)   :: work(*)
        integer,      intent(in)    :: lwork
        integer,      intent(out)   :: info
    end subroutine sormqr

    !> Multiplies in-place a matrix with an orthogonal matrix 
    !> defined as the product of elementary reflectors
    subroutine dormqr &
    (side, trans, m, n, k, A, ldA, tau, C, ldC, work, lwork, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        character(1), intent(in)    :: side
        character(1), intent(in)    :: trans
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        integer,      intent(in)    :: k
        real(WP),     intent(in)    :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(in)    :: tau(*)
        real(WP),     intent(inout) :: C(*)
        integer,      intent(in)    :: ldC
        real(WP),     intent(out)   :: work(*)
        integer,      intent(in)    :: lwork
        integer,      intent(out)   :: info
    end subroutine dormqr

    !> Generates a matrix with orthonormal columns, which is
    !> defined as the first columns of a product of elementary reflectors
    subroutine sorgqr &
    (m, n, k, A, ldA, tau, work, lwork, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        integer,      intent(in)    :: k
        real(WP),     intent(inout) :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(in)    :: tau(*)
        real(WP),     intent(out)   :: work(*)
        integer,      intent(in)    :: lwork
        integer,      intent(out)   :: info
    end subroutine sorgqr

    !> Generates a matrix with orthonormal columns, which is
    !> defined as the first columns of a product of elementary reflectors
    subroutine dorgqr &
    (m, n, k, A, ldA, tau, work, lwork, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        integer,      intent(in)    :: k
        real(WP),     intent(inout) :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(in)    :: tau(*)
        real(WP),     intent(out)   :: work(*)
        integer,      intent(in)    :: lwork
        integer,      intent(out)   :: info
    end subroutine dorgqr

    !> Multiplies in-place a matrix with an orthogonal matrix 
    !> defined as the product of elementary reflectors
    subroutine sormlq &
    (side, trans, m, n, k, A, ldA, tau, C, ldC, work, lwork, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        character(1), intent(in)    :: side
        character(1), intent(in)    :: trans
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        integer,      intent(in)    :: k
        real(WP),     intent(in)    :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(in)    :: tau(*)
        real(WP),     intent(inout) :: C(*)
        integer,      intent(in)    :: ldC
        real(WP),     intent(out)   :: work(*)
        integer,      intent(in)    :: lwork
        integer,      intent(out)   :: info
    end subroutine sormlq

    !> Multiplies in-place a matrix with an orthogonal matrix 
    !> defined as the product of elementary reflectors
    subroutine dormlq &
    (side, trans, m, n, k, A, ldA, tau, C, ldC, work, lwork, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        character(1), intent(in)    :: side
        character(1), intent(in)    :: trans
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        integer,      intent(in)    :: k
        real(WP),     intent(in)    :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(in)    :: tau(*)
        real(WP),     intent(inout) :: C(*)
        integer,      intent(in)    :: ldC
        real(WP),     intent(out)   :: work(*)
        integer,      intent(in)    :: lwork
        integer,      intent(out)   :: info
    end subroutine dormlq

    !> Generates a matrix with orthonormal rows, which is
    !> defined as the first rows of a product of elementary reflectors
    subroutine sorglq &
    (m, n, k, A, ldA, tau, work, lwork, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        integer,      intent(in)    :: k
        real(WP),     intent(inout) :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(in)    :: tau(*)
        real(WP),     intent(out)   :: work(*)
        integer,      intent(in)    :: lwork
        integer,      intent(out)   :: info
    end subroutine sorglq

    !> Generates a matrix with orthonormal rows, which is
    !> defined as the first rows of a product of elementary reflectors
    subroutine dorglq &
    (m, n, k, A, ldA, tau, work, lwork, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        integer,      intent(in)    :: k
        real(WP),     intent(inout) :: A(*)
        integer,      intent(in)    :: ldA
        real(WP),     intent(in)    :: tau(*)
        real(WP),     intent(out)   :: work(*)
        integer,      intent(in)    :: lwork
        integer,      intent(out)   :: info
    end subroutine dorglq

    !> Applies an elementary reflector to a matrix
    subroutine slarf &
    (side, m, n, v, incv, tau, C, ldC, work)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        character(1), intent(in)    :: side
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        real(WP),     intent(in)    :: v(*)
        integer,      intent(in)    :: incv
        real(WP),     intent(in)    :: tau
        real(WP),     intent(inout) :: C(*)
        integer,      intent(in)    :: ldC
        real(WP),     intent(out)   :: work(*)
    end subroutine slarf

    !> Applies an elementary reflector to a matrix
    subroutine dlarf &
    (side, m, n, v, incv, tau, C, ldC, work)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        character(1), intent(in)    :: side
        integer,      intent(in)    :: m
        integer,      intent(in)    :: n
        real(WP),     intent(in)    :: v(*)
        integer,      intent(in)    :: incv
        real(WP),     intent(in)    :: tau
        real(WP),     intent(inout) :: C(*)
        integer,      intent(in)    :: ldC
        real(WP),     intent(out)   :: work(*)
    end subroutine dlarf

    !> Generates an elementary reflector
    subroutine slarfg &
    (n, alpha, x, incx, tau)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,      intent(in)    :: n
        real(WP),     intent(inout) :: alpha
        real(WP),     intent(inout) :: x(*)
        integer,      intent(in)    :: incx
        real(WP),     intent(out)   :: tau
    end subroutine slarfg

    !> Generates an elementary reflector
    subroutine dlarfg &
    (n, alpha, x, incx, tau)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,      intent(in)    :: n
        real(WP),     intent(inout) :: alpha
        real(WP),     intent(inout) :: x(*)
        integer,      intent(in)    :: incx
        real(WP),     intent(out)   :: tau
    end subroutine dlarfg

    !> Copies a matrix
    subroutine slacpy &
    (uplo, m, n, A, ldA, B, ldB)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        character(1), intent(in)  :: uplo
        integer,      intent(in)  :: m
        integer,      intent(in)  :: n
        real(WP),     intent(in)  :: A(*)
        integer,      intent(in)  :: ldA
        real(WP),     intent(out) :: B(*)
        integer,      intent(in)  :: ldB
    end subroutine slacpy

    !> Copies a matrix
    subroutine dlacpy &
    (uplo, m, n, A, ldA, B, ldB)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        character(1), intent(in)  :: uplo
        integer,      intent(in)  :: m
        integer,      intent(in)  :: n
        real(WP),     intent(in)  :: A(*)
        integer,      intent(in)  :: ldA
        real(WP),     intent(out) :: B(*)
        integer,      intent(in)  :: ldB
    end subroutine dlacpy

    !> Initializes a matrix
    subroutine slaset &
    (uplo, m, n, alpha, beta, A, ldA)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        character(1), intent(in)  :: uplo
        integer,      intent(in)  :: m
        integer,      intent(in)  :: n
        real(WP),     intent(in)  :: alpha
        real(WP),     intent(in)  :: beta
        real(WP),     intent(out) :: A(*)
        integer,      intent(in)  :: ldA
    end subroutine slaset

    !> Initializes a matrix
    subroutine dlaset &
    (uplo, m, n, alpha, beta, A, ldA)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        character(1), intent(in)  :: uplo
        integer,      intent(in)  :: m
        integer,      intent(in)  :: n
        real(WP),     intent(in)  :: alpha
        real(WP),     intent(in)  :: beta
        real(WP),     intent(out) :: A(*)
        integer,      intent(in)  :: ldA
    end subroutine dlaset

    !> Sorts an array
    subroutine slasrt &
    (id, n, D, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        character(1), intent(in)    :: id
        integer,      intent(in)    :: n
        real(WP),     intent(inout) :: D(*)
        integer,      intent(out)   :: info
    end subroutine slasrt

    !> Sorts an array
    subroutine dlasrt &
    (id, n, D, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        character(1), intent(in)    :: id
        integer,      intent(in)    :: n
        real(WP),     intent(inout) :: D(*)
        integer,      intent(out)   :: info
    end subroutine dlasrt

end interface
end module maria_lapack_mod
