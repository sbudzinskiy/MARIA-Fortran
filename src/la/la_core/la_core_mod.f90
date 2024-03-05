!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_la_core_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Provides access to certain BLAS and LAPACK subroutines, their wrappers and extensions.
!----------------------------------------------------------------------------------------------------------------------------
module maria_la_core_mod
use maria_blas1_mod, only: &
    srotg, &
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
use maria_blas2_mod, only: &
    sgemv, &
    dgemv, &
    strmv, &
    dtrmv, &
    strsv, &
    dtrsv, &
    sger, &
    dger
use maria_blas3_mod, only: &
    sgemm, &
    dgemm, &
    strmm, &
    dtrmm, &
    strsm, &
    dtrsm
use maria_lapack_mod, only: &
    sgetrf, &
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
implicit none (type, external)

! Procedures
public :: sgesdd_q,       &
          dgesdd_q,       &
          sroundup_lwork, &
          droundup_lwork, &
          sgescal,        &
          dgescal,        &
          sgedotf,        &
          dgedotf,        &
          sgenrmf,        &
          dgenrmf,        &
          sgenrmc,        &
          dgenrmc,        &
          sdgmm,          &
          ddgmm,          &
          sgepiv,         &
          dgepiv,         &
          sorcangles,     &
          dorcangles

! Interfaces
public :: smatmul, &
          dmatmul

abstract interface
    !------------------------------------------------------------------------------------------------------------------------

    subroutine smatmul &
    (transB, m, n, k, alpha, B, ldB, beta, C, ldC, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        character(1), intent(in)                :: transB
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        integer,      intent(in)                :: k
        real(WP),     intent(in)                :: alpha
        real(WP),     intent(in),    contiguous :: B(:)
        integer,      intent(in)                :: ldB
        real(WP),     intent(in)                :: beta
        real(WP),     intent(inout), contiguous :: C(:)
        integer,      intent(in)                :: ldC
        integer,      intent(out)               :: info
    end subroutine smatmul

    !------------------------------------------------------------------------------------------------------------------------

    subroutine dmatmul &
    (transB, m, n, k, alpha, B, ldB, beta, C, ldC, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        character(1), intent(in)                :: transB
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        integer,      intent(in)                :: k
        real(WP),     intent(in)                :: alpha
        real(WP),     intent(in),    contiguous :: B(:)
        integer,      intent(in)                :: ldB
        real(WP),     intent(in)                :: beta
        real(WP),     intent(inout), contiguous :: C(:)
        integer,      intent(in)                :: ldC
        integer,      intent(out)               :: info
    end subroutine dmatmul

    !------------------------------------------------------------------------------------------------------------------------
end interface

interface
    !------------------------------------------------------------------------------------------------------------------------

    !> Computes a singular value decomposition.
    !>
    !> This is a wrapper over LAPACK's SGESDD, which takes an extra
    !> `liwork` argument and makes a query about the required size of `iwork`.
    module subroutine sgesdd_q &
    (jobZ, m, n, A, ldA, S, U, ldU, VT, ldVT, work, lwork, iwork, liwork, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
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
        !> If `lwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` and `iwork` arrays is calculated and returned
        !> in `work(1)` and `iwork(1)`, respectively.
        integer,      intent(in)                :: lwork
        integer,      intent(out),   contiguous :: iwork(:)
        !> If `liwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` and `iwork` arrays is calculated and returned
        !> in `work(1)` and `iwork(1)`, respectively.
        integer,      intent(in)                :: liwork
        integer,      intent(out)               :: info
    end subroutine sgesdd_q

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes a singular value decomposition.
    !>
    !> This is a wrapper over LAPACK's DGESDD, which takes an extra
    !> `liwork` argument and makes a query about the required size of `iwork`.
    module subroutine dgesdd_q &
    (jobZ, m, n, A, ldA, S, U, ldU, VT, ldVT, work, lwork, iwork, liwork, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
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
        !> If `lwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` and `iwork` arrays is calculated and returned
        !> in `work(1)` and `iwork(1)`, respectively.
        integer,      intent(in)                :: lwork
        integer,      intent(out),   contiguous :: iwork(:)
        !> If `liwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` and `iwork` arrays is calculated and returned
        !> in `work(1)` and `iwork(1)`, respectively.
        integer,      intent(in)                :: liwork
        integer,      intent(out)               :: info
    end subroutine dgesdd_q

    !------------------------------------------------------------------------------------------------------------------------

    !> Makes sure integer `lwork` is rounded up when cast to single precision.
    module function sroundup_lwork &
    (lwork)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer, intent(in) :: lwork
        real(WP)            :: sroundup_lwork
    end function sroundup_lwork

    !------------------------------------------------------------------------------------------------------------------------

    !> Makes sure integer `lwork` is rounded up when cast to double precision.
    module function droundup_lwork &
    (lwork)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer, intent(in) :: lwork
        real(WP)            :: droundup_lwork
    end function droundup_lwork

    !------------------------------------------------------------------------------------------------------------------------

    !> Scales a matrix by a constant:
    !> \[ A \gets \alpha A \]
    module subroutine sgescal &
    (m, n, alpha, A, ldA, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Number of rows in \( A \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,  intent(in)                :: m
        !> Number of columns in \( A \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,  intent(in)                :: n
        !> Scaling constant
        !>
        !> If \( \alpha = 0 \), \( A \) is not referenced.
        real(WP), intent(in)                :: alpha
        !> Matrix of dimensions \( \text{ldA} \times n \)
        real(WP), intent(inout), contiguous :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        integer,  intent(in)                :: ldA
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)               :: info
    end subroutine sgescal

    !------------------------------------------------------------------------------------------------------------------------

    !> Scales a matrix by a constant:
    !> \[ A \gets \alpha A \]
    module subroutine dgescal &
    (m, n, alpha, A, ldA, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Number of rows in \( A \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,  intent(in)                :: m
        !> Number of columns in \( A \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,  intent(in)                :: n
        !> Scaling constant
        !>
        !> If \( \alpha = 0 \), \( A \) is not referenced.
        real(WP), intent(in)                :: alpha
        !> Matrix of dimensions \( \text{ldA} \times n \)
        real(WP), intent(inout), contiguous :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        integer,  intent(in)                :: ldA
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)               :: info
    end subroutine dgescal

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the Frobenius inner product of two matrices:
    !> \[ \langle A, B, \rangle_F = \text{Tr}\big( \text{op}(A)^\top B \big) \]
    !> where 
    !> \[ \text{op}(A) = A \quad \text{or} \quad  \text{op}(A) = A^\top \]
    module function sgedotf &
    (transA, m, n, A, ldA, B, ldB, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Specifies the operation to be performed as follows:
        !> \[ \text{op}(A) = 
        !>    \begin{cases}
        !>      A, & \text{transA} = \text{'N' or 'n'} \\
        !>      A^\top, & \text{transA} = \text{'T' or 't'}
        !>    \end{cases}
        !> \]
        character(1), intent(in)             :: transA
        !> Number of rows in \( B \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,      intent(in)             :: m
        !> Number of columns in \( B \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,      intent(in)             :: n
        !> Matrix of dimensions
        !>
        !> - \( \text{ldA} \times n \) if `transA = "N" or "n"`
        !> - \( \text{ldA} \times m \) if `transA = "T" or "t"`
        real(WP),     intent(in), contiguous :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        !>
        !> - \( \text{ldA} \geq \max(1,m) \) if `transA = "N" or "n"`
        !> - \( \text{ldA} \geq \max(1,n) \) if `transA = "T" or "t"`
        integer,      intent(in)             :: ldA
        !> Matrix of dimensions \( \text{ldB} \times n \)
        real(WP),     intent(in), contiguous :: B(:)
        !> Leading dimension of \( B \)
        !>
        !> Possible values: \( \text{ldB} \geq \max(1, m) \)
        integer,      intent(in)             :: ldB
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,      intent(out)            :: info
        real(WP)                             :: sgedotf
    end function sgedotf

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the Frobenius inner product of two matrices:
    !> \[ \langle A, B, \rangle_F = \text{Tr}\big( \text{op}(A)^\top B \big) \]
    !> where 
    !> \[ \text{op}(A) = A \quad \text{or} \quad  \text{op}(A) = A^\top \]
    module function dgedotf &
    (transA, m, n, A, ldA, B, ldB, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Specifies the operation to be performed as follows:
        !> \[ \text{op}(A) = 
        !>    \begin{cases}
        !>      A, & \text{transA} = \text{'N' or 'n'} \\
        !>      A^\top, & \text{transA} = \text{'T' or 't'}
        !>    \end{cases}
        !> \]
        character(1), intent(in)             :: transA
        !> Number of rows in \( B \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,      intent(in)             :: m
        !> Number of columns in \( B \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,      intent(in)             :: n
        !> Matrix of dimensions
        !>
        !> - \( \text{ldA} \times n \) if `transA = "N" or "n"`
        !> - \( \text{ldA} \times m \) if `transA = "T" or "t"`
        real(WP),     intent(in), contiguous :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        !>
        !> - \( \text{ldA} \geq \max(1,m) \) if `transA = "N" or "n"`
        !> - \( \text{ldA} \geq \max(1,n) \) if `transA = "T" or "t"`
        integer,      intent(in)             :: ldA
        !> Matrix of dimensions \( \text{ldB} \times n \)
        real(WP),     intent(in), contiguous :: B(:)
        !> Leading dimension of \( B \)
        !>
        !> Possible values: \( \text{ldB} \geq \max(1, m) \)
        integer,      intent(in)             :: ldB
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,      intent(out)            :: info
        real(WP)                             :: dgedotf
    end function dgedotf

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the Frobenius norm of a matrix:
    !> \[ \| A \|_F = \text{Tr}\big( A^\top A \big) \]
    module function sgenrmf &
    (m, n, A, ldA, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Number of rows in \( A \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,      intent(in)             :: m
        !> Number of columns in \( A \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,      intent(in)             :: n
        !> Matrix of dimensions
        !>
        !> - \( \text{ldA} \times n \)
        real(WP),     intent(in), contiguous :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        integer,      intent(in)             :: ldA
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,      intent(out)            :: info
        real(WP)                             :: sgenrmf
    end function sgenrmf

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the Frobenius norm of a matrix:
    !> \[ \| A \|_F = \text{Tr}\big( A^\top A \big) \]
    module function dgenrmf &
    (m, n, A, ldA, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Number of rows in \( A \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,      intent(in)             :: m
        !> Number of columns in \( A \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,      intent(in)             :: n
        !> Matrix of dimensions
        !>
        !> - \( \text{ldA} \times n \)
        real(WP),     intent(in), contiguous :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        integer,      intent(in)             :: ldA
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,      intent(out)            :: info
        real(WP)                             :: dgenrmf
    end function dgenrmf

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the Chebyshev (maximum) norm of a matrix and the position of the largest element:
    !> \[ \| A \|_C = \max_{i, j} |A(i,j)| \]
    module function sgenrmc &
    (m, n, A, ldA, info, pos)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Number of rows in \( A \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,      intent(in)             :: m
        !> Number of columns in \( A \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,      intent(in)             :: n
        !> Matrix of dimensions
        !>
        !> - \( \text{ldA} \times n \)
        real(WP),     intent(in), contiguous :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        integer,      intent(in)             :: ldA
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,      intent(out)            :: info
        !> Indices of the largest element
        integer,      intent(out), optional  :: pos(2)
        real(WP)                             :: sgenrmc
    end function sgenrmc

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the Chebyshev (maximum) norm of a matrix and the position of the largest element:
    !> \[ \| A \|_C = \max_{i, j} |A(i,j)| \]
    module function dgenrmc &
    (m, n, A, ldA, info, pos)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Number of rows in \( A \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,      intent(in)             :: m
        !> Number of columns in \( A \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,      intent(in)             :: n
        !> Matrix of dimensions
        !>
        !> - \( \text{ldA} \times n \)
        real(WP),     intent(in), contiguous :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        integer,      intent(in)             :: ldA
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,      intent(out)            :: info
        !> Indices of the largest element
        integer,      intent(out), optional  :: pos(2)
        real(WP)                             :: dgenrmc
    end function dgenrmc

    !------------------------------------------------------------------------------------------------------------------------

    !> Multiplies a matrix by a diagonal matrix on the left or on the right
    !>
    !> \[ A \gets DA \quad \text{or} \quad A \gets AD \]
    module subroutine sdgmm &
    (side, m, n, A, lda, D, incd, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Specifies, whether the diagonal matrix is multiplied on the left or on the right
        !   ='L': computes D*A
        !   ='R': computes A*D
        character(1), intent(in)                :: side
        !> Number of rows in \( A \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,      intent(in)                :: m
        !> Number of columns in \( A \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,      intent(in)                :: n
        !> Matrix of dimensions
        !>
        !> - \( \text{ldA} \times n \)
        real(WP),     intent(inout), contiguous :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        integer,      intent(in)                :: ldA
        !> Array of dimension
        !>  side='L': [1 + (m-1)*abs(incd)]
        !>  side='R': [1 + (n-1)*abs(incd)]
        real(WP),     intent(in),    contiguous :: D(:)
        !> Stride for D
        integer,      intent(in)                :: incD
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,      intent(out)               :: info
    end subroutine sdgmm

    !------------------------------------------------------------------------------------------------------------------------

    !> Multiplies a matrix by a diagonal matrix on the left or on the right
    !>
    !> \[ A \gets DA \quad \text{or} \quad A \gets AD \]
    module subroutine ddgmm &
    (side, m, n, A, lda, D, incd, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Specifies, whether the diagonal matrix is multiplied on the left or on the right
        !   ='L': computes D*A
        !   ='R': computes A*D
        character(1), intent(in)                :: side
        !> Number of rows in \( A \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,      intent(in)                :: m
        !> Number of columns in \( A \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,      intent(in)                :: n
        !> Matrix of dimensions
        !>
        !> - \( \text{ldA} \times n \)
        real(WP),     intent(inout), contiguous :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        integer,      intent(in)                :: ldA
        !> Array of dimension
        !>  side='L': [1 + (m-1)*abs(incd)]
        !>  side='R': [1 + (n-1)*abs(incd)]
        real(WP),     intent(in),    contiguous :: D(:)
        !> Stride for D
        integer,      intent(in)                :: incD
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,      intent(out)               :: info
    end subroutine ddgmm

    !------------------------------------------------------------------------------------------------------------------------

    !> Applies a series of elementary pivot permutations to rows or columns of a matrix
    !> and computes their new ordering.
    !>
    !> The pivots are stored in an array as
    !> \[ \text{ipiv}(1) = j_{\text{from}}, \ldots, \text{ipiv}(\text{to}-\text{from}+1) = j_{\text{to}} \]
    module subroutine sgepiv &
    (what, dir, m, n, A, lda, from, to, ipiv, info, order)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Specifies, whether to permute rows or columns
        !> ='R': rows are permuted
        !> ='C': columns are permuted
        character(1), intent(in)                          :: what
        !> Specifies the order, in which to apply the pivots
        !> ='F': pivots are applied in forward ordering
        !> ='B': pivots are applied in backward ordering
        character(1), intent(in)                          :: dir
        !> Number of rows in \( A \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,      intent(in)                          :: m
        !> Number of columns in \( A \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,      intent(in)                          :: n
        !> Matrix of dimensions
        !>
        !> - \( \text{ldA} \times n \)
        real(WP),     intent(inout), contiguous           :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        integer,      intent(in)                          :: ldA
        !> Index, correspomding to the first pivot in the array
        !>
        !> Possible values: \( \text{from} > 0 \)
        integer,      intent(in)                          :: from
        !> Index, correspomding to the last pivot in the array
        !>
        !> Possible values: \( \text{to} > 0 \)
        !> - \( \text{from} \leq \text{to} \leq m
        !> - \( \text{from} \leq \text{to} \leq n
        integer,      intent(in)                          :: to
        !> Array of pivots of size \( \text{to} - \text{from} + 1 \)
        integer,      intent(in),    contiguous           :: ipiv(:)
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        !> - `info > 0`: ipiv[info] is out of range
        integer,      intent(out)                         :: info
        !> Array of size
        !> - \( m \)
        !> - \( n \)
        !>
        !> On entry, the initial ordering of rows/columns
        !> On exit, the new ordering
        integer,      intent(inout), contiguous, optional :: order(:)
    end subroutine sgepiv

    !------------------------------------------------------------------------------------------------------------------------

    !> Applies a series of elementary pivot permutations to rows or columns of a matrix
    !> and computes their new ordering.
    !>
    !> The pivots are stored in an array as
    !> \[ \text{ipiv}(1) = j_{\text{from}}, \ldots, \text{ipiv}(\text{to}-\text{from}+1) = j_{\text{to}} \]
    module subroutine dgepiv &
    (what, dir, m, n, A, lda, from, to, ipiv, info, order)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Specifies, whether to permute rows or columns
        !> ='R': rows are permuted
        !> ='C': columns are permuted
        character(1), intent(in)                          :: what
        !> Specifies the order, in which to apply the pivots
        !> ='F': pivots are applied in forward ordering
        !> ='B': pivots are applied in backward ordering
        character(1), intent(in)                          :: dir
        !> Number of rows in \( A \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,      intent(in)                          :: m
        !> Number of columns in \( A \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,      intent(in)                          :: n
        !> Matrix of dimensions
        !>
        !> - \( \text{ldA} \times n \)
        real(WP),     intent(inout), contiguous           :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        integer,      intent(in)                          :: ldA
        !> Index, correspomding to the first pivot in the array
        !>
        !> Possible values: \( \text{from} > 0 \)
        integer,      intent(in)                          :: from
        !> Index, correspomding to the last pivot in the array
        !>
        !> Possible values: \( \text{to} > 0 \)
        !> - \( \text{from} \leq \text{to} \leq m
        !> - \( \text{from} \leq \text{to} \leq n
        integer,      intent(in)                          :: to
        !> Array of pivots of size \( \text{to} - \text{from} + 1 \)
        integer,      intent(in),    contiguous           :: ipiv(:)
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        !> - `info > 0`: ipiv[info] is out of range
        integer,      intent(out)                         :: info
        !> Array of size
        !> - \( m \)
        !> - \( n \)
        !>
        !> On entry, the initial ordering of rows/columns
        !> On exit, the new ordering
        integer,      intent(inout), contiguous, optional :: order(:)
    end subroutine dgepiv

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the cosines of the canonical angles between two subspaces
    !> given by their orthonormal bases. 
    module subroutine sorcangles &
    (what, m, n, U, ldU, V, ldV, cosines, work, lwork, iwork, liwork, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Specifies if the subspaces correspond to the column or row spans
        !>
        !> Possible values: 'c', 'C', 'r', 'R'
        character(1), intent(in)              :: what
        !> Number of rows in \( U \) and \( V \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,      intent(in)              :: m
        !> Number of columns in \( U \) and \( V \)
        !>
        !> Possible values: \( n \geq 0 \)
        !> - what='c': \( 0 \leq n \leq m \)
        !> - what='r': \( n \geq m \)
        integer,      intent(in)              :: n
        !> Matrix of dimensions \( m \times n \)
        !> with orthonormal rows or columns
        real(WP),     intent(in),  contiguous :: U(:)
        !> Leading dimension of \( U \)
        !>
        !> Possible values: \( \text{ldU} \geq \max(1, m) \)
        integer,      intent(in)              :: ldU
        !> Matrix of dimensions \( m \times n \)
        !> with orthonormal rows or columns
        real(WP),     intent(in),  contiguous :: V(:)
        !> Leading dimension of \( V \)
        !>
        !> Possible values: \( \text{ldV} \geq \max(1, m) \)
        integer,      intent(in)              :: ldV
        !> The cosines of the canonical angles in decreasing order
        real(WP),     intent(out), contiguous :: cosines(:)
        !> Real workspace of size `lwork`
        real(WP),     intent(out), contiguous :: work(:)
        !> Size of the real workspace.
        !>
        !> If `lwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` and `iwork` arrays is calculated and returned
        !> in `work(1)` and `iwork(1)`, respectively.
        integer,      intent(in)              :: lwork
        !> Integer workspace of size `liwork`
        integer,      intent(out), contiguous :: iwork(:)
        !> Size of the integer workspace
        !>
        !> If `liwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` and `iwork` arrays is calculated and returned
        !> in `work(1)` and `iwork(1)`, respectively.
        integer,      intent(in)              :: liwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        !> - `info > 0`: the computation of SVD did not converge
        integer,      intent(out)             :: info
    end subroutine sorcangles

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the cosines of the canonical angles between two subspaces
    !> given by their orthonormal bases. 
    module subroutine dorcangles &
    (what, m, n, U, ldU, V, ldV, cosines, work, lwork, iwork, liwork, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Specifies if the subspaces correspond to the column or row spans
        !>
        !> Possible values: 'c', 'C', 'r', 'R'
        character(1), intent(in)              :: what
        !> Number of rows in \( U \) and \( V \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,      intent(in)              :: m
        !> Number of columns in \( U \) and \( V \)
        !>
        !> Possible values: \( n \geq 0 \)
        !> - what='c': \( 0 \leq n \leq m \)
        !> - what='r': \( n \geq m \)
        integer,      intent(in)              :: n
        !> Matrix of dimensions \( m \times n \)
        !> with orthonormal rows or columns
        real(WP),     intent(in),  contiguous :: U(:)
        !> Leading dimension of \( U \)
        !>
        !> Possible values: \( \text{ldU} \geq \max(1, m) \)
        integer,      intent(in)              :: ldU
        !> Matrix of dimensions \( m \times n \)
        !> with orthonormal rows or columns
        real(WP),     intent(in),  contiguous :: V(:)
        !> Leading dimension of \( V \)
        !>
        !> Possible values: \( \text{ldV} \geq \max(1, m) \)
        integer,      intent(in)              :: ldV
        !> The cosines of the canonical angles in decreasing order
        real(WP),     intent(out), contiguous :: cosines(:)
        !> Real workspace of size `lwork`
        real(WP),     intent(out), contiguous :: work(:)
        !> Size of the real workspace.
        !>
        !> If `lwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` and `iwork` arrays is calculated and returned
        !> in `work(1)` and `iwork(1)`, respectively.
        integer,      intent(in)              :: lwork
        !> Integer workspace of size `liwork`
        integer,      intent(out), contiguous :: iwork(:)
        !> Size of the integer workspace
        !>
        !> If `liwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` and `iwork` arrays is calculated and returned
        !> in `work(1)` and `iwork(1)`, respectively.
        integer,      intent(in)              :: liwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        !> - `info > 0`: the computation of SVD did not converge
        integer,      intent(out)             :: info
    end subroutine dorcangles

    !------------------------------------------------------------------------------------------------------------------------

end interface

end module maria_la_core_mod
