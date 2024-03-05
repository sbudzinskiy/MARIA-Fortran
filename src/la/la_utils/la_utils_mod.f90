!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_la_utils_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Provides access to utility subroutines for comparing vectors and matrices, generating certain random matrices.
!----------------------------------------------------------------------------------------------------------------------------
module maria_la_utils_mod
implicit none (type, external)

! Interfaces
public :: all_close,   &
          all_const,   &
          geall_close, &
          geall_const
! Procedures
public :: sall_close,   &
          dall_close,   &
          sall_const,   &
          dall_const,   &
          sgeall_close, &
          dgeall_close, &
          sgeall_const, &
          dgeall_const, &
          srandort,     &
          drandort,     &
          srandsvd,     &
          drandsvd

!> Checks if two floating-point vectors are elementwise approximately equal 
!> with given absolute and/or relative tolerance.
interface all_close
    procedure sall_close
    procedure dall_close
end interface all_close

!> Checks if a floating-point vector is elementwise approximately equal to a given constant
!> with given absolute and/or relative tolerance.
interface all_const
    procedure sall_const
    procedure dall_const
end interface all_const

!> Checks if two floating-point matrices are elementwise approximately equal 
!> with given absolute and/or relative tolerance.
interface geall_close
    procedure sgeall_close
    procedure dgeall_close
end interface geall_close

!> Checks if diagonal and off-diagonal elements of a floating-point matrix are approximately equal to two given constants
!> with given absolute and/or relative tolerance.
interface geall_const
    procedure sgeall_const
    procedure dgeall_const
end interface geall_const

interface
    !------------------------------------------------------------------------------------------------------------------------

    !> Checks if two single-precision vectors are elementwise approximately equal
    !> with given absolute and/or relative tolerance:
    !>
    !> \[ |x(i) - y(i)| < \max\big( \text{atol}, \text{rtol} \cdot \max(|x(i)|, |y(i)|)  \big), \quad i = 1, \ldots, n \].
    module function sall_close &
    (n, x, incx, y, incy, info, atol, rtol)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Number of elements of the vectors \( x \) and \( y \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,  intent(in)             :: n
        !> Vector of size
        !> \[ 1 + (n-1)|\text{incx}| \]
        real(WP), intent(in), contiguous :: x(:)
        !> Stride for the elements of \( x \)
        integer,  intent(in)             :: incx
        !> Vector of size
        !> \[ 1 + (n-1)|\text{incy}| \]
        real(WP), intent(in), contiguous :: y(:)
        !> Stride for the elements of \( y \)
        integer,  intent(in)             :: incy
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)            :: info
        !> Absolute tolerance
        !>
        !> Possible values: \( \text{atol} \geq 0 \)
        !>
        !> Default value: [[maria_constants_mod(module):S_MACHTOL(variable)]]
        real(WP), intent(in), optional   :: atol
        !> Relative tolerance
        !>
        !> Possible values: \( \text{rtol} \geq 0 \)
        !>
        !> Default value: \( 0 \)
        real(WP), intent(in), optional   :: rtol
        logical                          :: sall_close
    end function sall_close

    !------------------------------------------------------------------------------------------------------------------------

    !> Checks if two double-precision vectors are elementwise approximately equal
    !> with given absolute and/or relative tolerance:
    !>
    !> \[ |x(i) - y(i)| < \max\big( \text{atol}, \text{rtol} \cdot \max(|x(i)|, |y(i)|)  \big), \quad i = 1, \ldots, n \].
    module function dall_close &
    (n, x, incx, y, incy, info, atol, rtol)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Number of elements of the vectors \( x \) and \( y \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,  intent(in)             :: n
        !> Vector of size
        !> \[ 1 + (n-1)|\text{incx}| \]
        real(WP), intent(in), contiguous :: x(:)
        !> Stride for the elements of \( x \)
        integer,  intent(in)             :: incx
        !> Vector of size
        !> \[ 1 + (n-1)|\text{incy}| \]
        real(WP), intent(in), contiguous :: y(:)
        !> Stride for the elements of \( y \)
        integer,  intent(in)             :: incy
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)            :: info
        !> Absolute tolerance
        !>
        !> Possible values: \( \text{atol} \geq 0 \)
        !>
        !> Default value: [[maria_constants_mod(module):D_MACHTOL(variable)]]
        real(WP), intent(in), optional   :: atol
        !> Relative tolerance
        !>
        !> Possible values: \( \text{rtol} \geq 0 \)
        !>
        !> Default value: \( 0 \)
        real(WP), intent(in), optional   :: rtol
        logical                          :: dall_close
    end function dall_close

    !------------------------------------------------------------------------------------------------------------------------

    !> Checks if a single-precision vector is elementwise approximately equal to a given constant
    !> with given absolute and/or relative tolerance:
    !>
    !> \[ |x(i) - \alpha| < \max\big( \text{atol}, \text{rtol} \cdot \max(|x(i)|, |\alpha|)  \big), \quad i = 1, \ldots, n \].
    module function sall_const &
    (n, x, incx, alpha, info, atol, rtol)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Number of elements of the vector \( x \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,  intent(in)             :: n
        !> Vector of size
        !> \[ 1 + (n-1)|\text{incx}| \]
        real(WP), intent(in), contiguous :: x(:)
        !> Stride for the elements of \( x \)
        integer,  intent(in)             :: incx
        !> Reference value
        real(WP), intent(in)             :: alpha
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)            :: info
        !> Absolute tolerance
        !>
        !> Possible values: \( \text{atol} \geq 0 \)
        !>
        !> Default value: [[maria_constants_mod(module):S_MACHTOL(variable)]]
        real(WP), intent(in), optional   :: atol
        !> Relative tolerance
        !>
        !> Possible values: \( \text{rtol} \geq 0 \)
        !>
        !> Default value: \( 0 \)
        real(WP), intent(in), optional   :: rtol
        logical                          :: sall_const
    end function sall_const
 
    !------------------------------------------------------------------------------------------------------------------------
   
    !> Checks if a double-precision vector is elementwise approximately equal to a given constant
    !> with given absolute tolerance:
    !>
    !> \[ |x(i) - \alpha| < \max\big( \text{atol}, \text{rtol} \cdot \max(|x(i)|, |\alpha|)  \big), \quad i = 1, \ldots, n \].
    module function dall_const &
    (n, x, incx, alpha, info, atol, rtol)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Number of elements of the vector \( x \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,  intent(in)             :: n
        !> Vector of size
        !> \[ 1 + (n-1)|\text{incx}| \]
        real(WP), intent(in), contiguous :: x(:)
        !> Stride for the elements of \( x \)
        integer,  intent(in)             :: incx
        !> Reference value
        real(WP), intent(in)             :: alpha
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)            :: info
        !> Absolute tolerance
        !>
        !> Possible values: \( \text{atol} \geq 0 \)
        !>
        !> Default value: [[maria_constants_mod(module):S_MACHTOL(variable)]]
        real(WP), intent(in), optional   :: atol
        !> Relative tolerance
        !>
        !> Possible values: \( \text{rtol} \geq 0 \)
        !>
        !> Default value: \( 0 \)
        real(WP), intent(in), optional   :: rtol
        logical                          :: dall_const
    end function dall_const

    !------------------------------------------------------------------------------------------------------------------------

    !> Checks if two single-precision matrices are elementwise approximately equal
    !> with given absolute and/or relative tolerance:
    !>
    !> \[ |A(i, j) - B(i, j)| < \max\big( \text{atol}, \text{rtol} \cdot \max(|A(i,j)|, |B(i,j)|)  \big), \quad i = 1, \ldots, m ,
    !   \quad j = 1, \ldots, n\].
    module function sgeall_close &
    (m, n, A, ldA, B, ldB, info, atol, rtol)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Number of rows of the matrices \( A \) and \( B \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,  intent(in)             :: m
        !> Number of columns of the matrices \( A \) and \( B \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,  intent(in)             :: n
        !> Matrix of size
        !> \[ \text{lda} \times n \]
        real(WP), intent(in), contiguous :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        integer,  intent(in)             :: ldA
        !> Matrix of size
        !> \[ \text{ldb} \times n \]
        real(WP), intent(in), contiguous :: B(:)
        !> Leading dimension of \( B \)
        !>
        !> Possible values: \( \text{ldB} \geq \max(1, m) \)
        integer,  intent(in)             :: ldB
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)            :: info
        !> Absolute tolerance
        !>
        !> Possible values: \( \text{atol} \geq 0 \)
        !>
        !> Default value: [[maria_constants_mod(module):S_MACHTOL(variable)]]
        real(WP), intent(in), optional   :: atol
        !> Relative tolerance
        !>
        !> Possible values: \( \text{rtol} \geq 0 \)
        !>
        !> Default value: \( 0 \)
        real(WP), intent(in), optional   :: rtol
        logical                          :: sgeall_close
    end function sgeall_close

    !------------------------------------------------------------------------------------------------------------------------

    !> Checks if two double-precision matrices are elementwise approximately equal
    !> with given absolute and/or relative tolerance:
    !>
    !> \[ |A(i, j) - B(i, j)| < \max\big( \text{atol}, \text{rtol} \cdot \max(|A(i,j)|, |B(i,j)|)  \big), \quad i = 1, \ldots, m ,
    !   \quad j = 1, \ldots, n\].
    module function dgeall_close &
    (m, n, A, ldA, B, ldB, info, atol, rtol)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Number of rows of the matrices \( A \) and \( B \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,  intent(in)             :: m
        !> Number of columns of the matrices \( A \) and \( B \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,  intent(in)             :: n
        !> Matrix of size
        !> \[ \text{lda} \times n \]
        real(WP), intent(in), contiguous :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        integer,  intent(in)             :: ldA
        !> Matrix of size
        !> \[ \text{ldb} \times n \]
        real(WP), intent(in), contiguous :: B(:)
        !> Leading dimension of \( B \)
        !>
        !> Possible values: \( \text{ldB} \geq \max(1, m) \)
        integer,  intent(in)             :: ldB
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)            :: info
        !> Absolute tolerance
        !>
        !> Possible values: \( \text{atol} \geq 0 \)
        !>
        !> Default value: [[maria_constants_mod(module):S_MACHTOL(variable)]]
        real(WP), intent(in), optional   :: atol
        !> Relative tolerance
        !>
        !> Possible values: \( \text{rtol} \geq 0 \)
        !>
        !> Default value: \( 0 \)
        real(WP), intent(in), optional   :: rtol
        logical                          :: dgeall_close
    end function dgeall_close

    !------------------------------------------------------------------------------------------------------------------------

    !> Checks if diagonal and off-diagonal elements of a single-precision matrix
    !> are approximately equal to two given constants with given absolute and/or relative tolerance.
    !>
    !> \[ |A(i, i) - \beta| < \max\big( \text{atol}, \text{rtol} \cdot \max(|A(i,i)|, |\beta|)  \big), \quad i = 1, \ldots, \min(m,n) \].
    !> \[ |A(i, j) - \alpha| < \max\big( \text{atol}, \text{rtol} \cdot \max(|A(i,j)|, |\alpha|)  \big), \quad i = 1, \ldots, m, \quad j = 1, \ldots\ n, \quad j \neq i\].
    module function sgeall_const &
    (m, n, A, ldA, alpha, beta, info, atol, rtol)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Number of rows of the matrices \( A \) and \( B \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,  intent(in)             :: m
        !> Number of columns of the matrices \( A \) and \( B \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,  intent(in)             :: n
        !> Matrix of size
        !> \[ \text{lda} \times n \]
        real(WP), intent(in), contiguous :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        integer,  intent(in)             :: ldA
        !> Reference off-diagonal value
        real(WP), intent(in)             :: alpha
        !> Reference diagonal value
        real(WP), intent(in)             :: beta
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)            :: info
        !> Absolute tolerance
        !>
        !> Possible values: \( \text{atol} \geq 0 \)
        !>
        !> Default value: [[maria_constants_mod(module):S_MACHTOL(variable)]]
        real(WP), intent(in), optional   :: atol
        !> Relative tolerance
        !>
        !> Possible values: \( \text{rtol} \geq 0 \)
        !>
        !> Default value: \( 0 \)
        real(WP), intent(in), optional   :: rtol
        logical                          :: sgeall_const
    end function sgeall_const

    !------------------------------------------------------------------------------------------------------------------------

    !> Checks if diagonal and off-diagonal elements of a double-precision matrix
    !> are approximately equal to two given constants with given absolute and/or relative tolerance.
    !>
    !> \[ |A(i, i) - \beta| < \max\big( \text{atol}, \text{rtol} \cdot \max(|A(i,i)|, |\beta|)  \big), \quad i = 1, \ldots, \min(m,n) \].
    !> \[ |A(i, j) - \alpha| < \max\big( \text{atol}, \text{rtol} \cdot \max(|A(i,j)|, |\alpha|)  \big), \quad i = 1, \ldots, m, \quad j = 1, \ldots\ n, \quad j \neq i\].
    module function dgeall_const &
    (m, n, A, ldA, alpha, beta, info, atol, rtol)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Number of rows of the matrices \( A \) and \( B \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,  intent(in)             :: m
        !> Number of columns of the matrices \( A \) and \( B \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,  intent(in)             :: n
        !> Matrix of size
        !> \[ \text{lda} \times n \]
        real(WP), intent(in), contiguous :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        integer,  intent(in)             :: ldA
        !> Reference off-diagonal value
        real(WP), intent(in)             :: alpha
        !> Reference diagonal value
        real(WP), intent(in)             :: beta
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)            :: info
        !> Absolute tolerance
        !>
        !> Possible values: \( \text{atol} \geq 0 \)
        !>
        !> Default value: [[maria_constants_mod(module):S_MACHTOL(variable)]]
        real(WP), intent(in), optional   :: atol
        !> Relative tolerance
        !>
        !> Possible values: \( \text{rtol} \geq 0 \)
        !>
        !> Default value: \( 0 \)
        real(WP), intent(in), optional   :: rtol
        logical                          :: dgeall_const
    end function dgeall_const

    !------------------------------------------------------------------------------------------------------------------------

    !> Generates a random matrix with orthonormal columns or rows from the uniform distribution
    !> with respect to the Haar measure
    module subroutine srandort &
    (rng, m, n, Q, ldQ, work, lwork, info)
    use maria_prng_mod,  only: &
        prng
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Pseudorandom number generator
        class(prng), intent(in)              :: rng
        !> Number of rows in \( Q \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,     intent(in)              :: m
        !> Number of columns in \( Q \)
        integer,     intent(in)              :: n
        !> Matrix of size \( \text{lda} \times n \)
        real(WP),    intent(out), contiguous :: Q(:)
        !> Leading dimension of \( Q \)
        !>
        !> Possible values: \( \text{ldQ} \geq \max(1, m) \)
        integer,     intent(in)              :: ldQ
        !> Real workspace of size `lwork`
        real(WP),    intent(out), contiguous :: work(:)
        !> Size of the real workspace.
        !>
        !> If `lwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` array is calculated and returned
        !> in `work(1)`.
        integer,     intent(in)              :: lwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,     intent(out)             :: info
    end subroutine srandort

    !------------------------------------------------------------------------------------------------------------------------

    !> Generates a random matrix with orthonormal columns or rows from the uniform distribution
    !> with respect to the Haar measure
    module subroutine drandort &
    (rng, m, n, Q, ldQ, work, lwork, info)
    use maria_prng_mod,  only: &
        prng
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Pseudorandom number generator
        class(prng), intent(in)              :: rng
        !> Number of rows in \( Q \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,     intent(in)              :: m
        !> Number of columns in \( Q \)
        integer,     intent(in)              :: n
        !> Matrix of size \( \text{lda} \times n \)
        real(WP),    intent(out), contiguous :: Q(:)
        !> Leading dimension of \( Q \)
        !>
        !> Possible values: \( \text{ldQ} \geq \max(1, m) \)
        integer,     intent(in)              :: ldQ
        !> Real workspace of size `lwork`
        real(WP),    intent(out), contiguous :: work(:)
        !> Size of the real workspace.
        !>
        !> If `lwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` array is calculated and returned
        !> in `work(1)`.
        integer,     intent(in)              :: lwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,     intent(out)             :: info
    end subroutine drandort

    !------------------------------------------------------------------------------------------------------------------------

    !> Generates a random matrix with prescribed singular values
    module subroutine srandsvd &
    (rng, m, n, A, ldA, S, work, lwork, info)
    use maria_prng_mod,  only: &
        prng
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Pseudorandom number generator
        class(prng), intent(in)              :: rng
        !> Number of rows in \( A \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,     intent(in)              :: m
        !> Number of columns in \( A \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,     intent(in)              :: n
        !> Matrix of size \( \text{lda} \times n \)
        real(WP),    intent(out), contiguous :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        integer,     intent(in)              :: ldA
        !> Array of size \( \min(m,n) \), the singular values
        real(WP),    intent(in),  contiguous :: S(:)
        !> Real workspace of size `lwork`
        real(WP),    intent(out), contiguous :: work(:)
        !> Size of the real workspace.
        !>
        !> If `lwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` array is calculated and returned
        !> in `work(1)`.
        integer,     intent(in)              :: lwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,     intent(out)             :: info
    end subroutine srandsvd

    !------------------------------------------------------------------------------------------------------------------------

    !> Generates a random matrix with prescribed singular values
    module subroutine drandsvd &
    (rng, m, n, A, ldA, S, work, lwork, info)
    use maria_prng_mod,  only: &
        prng
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Pseudorandom number generator
        class(prng), intent(in)              :: rng
        !> Number of rows in \( A \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,     intent(in)              :: m
        !> Number of columns in \( A \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,     intent(in)              :: n
        !> Matrix of size \( \text{lda} \times n \)
        real(WP),    intent(out), contiguous :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        integer,     intent(in)              :: ldA
        !> Array of size \( \min(m,n) \), the singular values
        real(WP),    intent(in),  contiguous :: S(:)
        !> Real workspace of size `lwork`
        real(WP),    intent(out), contiguous :: work(:)
        !> Size of the real workspace.
        !>
        !> If `lwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` array is calculated and returned
        !> in `work(1)`.
        integer,     intent(in)              :: lwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,     intent(out)             :: info
    end subroutine drandsvd

    !------------------------------------------------------------------------------------------------------------------------
end interface

end module maria_la_utils_mod
