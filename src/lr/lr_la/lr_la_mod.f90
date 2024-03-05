!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_lr_la_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Provides access to linear algebra subroutines for factorized matrices.
!----------------------------------------------------------------------------------------------------------------------------
module maria_lr_la_mod
implicit none (type, external)

! Procedures
public :: slrval,       &
          dlrval,       &
          slr2full,     &
          dlr2full,     &
          slrdotf,      &
          dlrdotf,      &
          slrnrmf,      &
          dlrnrmf,      &
          slrnrmf_diff, &
          dlrnrmf_diff, &
          lrort_rank,   &
          slrort,       &
          dlrort,       &
          slrsvd_ort,   &
          dlrsvd_ort

interface
    !------------------------------------------------------------------------------------------------------------------------

    !> Computes an element of a factorized matrix
    module function slrval &
    (m, n, i, j, r, U, ldU, VT, ldVT, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,  intent(in)              :: m
        integer,  intent(in)              :: n
        integer,  intent(in)              :: i
        integer,  intent(in)              :: j
        integer,  intent(in)              :: r
        real(WP), intent(in),  contiguous :: U(:)
        integer,  intent(in)              :: ldU
        real(WP), intent(in),  contiguous :: VT(:)
        integer,  intent(in)              :: ldVT
        integer,  intent(out)             :: info
        real(WP)                          :: slrval
    end function slrval

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes an element of a factorized matrix
    module function dlrval &
    (m, n, i, j, r, U, ldU, VT, ldVT, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,  intent(in)              :: m
        integer,  intent(in)              :: n
        integer,  intent(in)              :: i
        integer,  intent(in)              :: j
        integer,  intent(in)              :: r
        real(WP), intent(in),  contiguous :: U(:)
        integer,  intent(in)              :: ldU
        real(WP), intent(in),  contiguous :: VT(:)
        integer,  intent(in)              :: ldVT
        integer,  intent(out)             :: info
        real(WP)                          :: dlrval
    end function dlrval

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the full representation of a factorized matrix
    module subroutine slr2full &
    (m, n, r, U, ldU, VT, ldVT, A, lda, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,  intent(in)              :: m
        integer,  intent(in)              :: n
        integer,  intent(in)              :: r
        real(WP), intent(in),  contiguous :: U(:)
        integer,  intent(in)              :: ldU
        real(WP), intent(in),  contiguous :: VT(:)
        integer,  intent(in)              :: ldVT
        real(WP), intent(out), contiguous :: A(:)
        integer,  intent(in)              :: ldA
        integer,  intent(out)             :: info
    end subroutine slr2full

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the full representation of a factorized matrix
    module subroutine dlr2full &
    (m, n, r, U, ldU, VT, ldVT, A, lda, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,  intent(in)              :: m
        integer,  intent(in)              :: n
        integer,  intent(in)              :: r
        real(WP), intent(in),  contiguous :: U(:)
        integer,  intent(in)              :: ldU
        real(WP), intent(in),  contiguous :: VT(:)
        integer,  intent(in)              :: ldVT
        real(WP), intent(out), contiguous :: A(:)
        integer,  intent(in)              :: ldA
        integer,  intent(out)             :: info
    end subroutine dlr2full

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the Frobenius inner product of two factorized matrices
    module function slrdotf &
    (m, n, r, U, ldU, VT, ldVT, k, A, ldA, BT, ldBT, work, lwork, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Number of rows in \( U \) and \( A \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,  intent(in)              :: m
        !> Number of columns in \( V^\top \) and \( B^\top \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,  intent(in)              :: n
        !> Factorization rank of \( U V^\top \), i.e. 
        !> number of columns in \( U \) and number of rows in \( V^\top \)
        !>
        !> Possible values: \( r \geq 0 \)
        integer,  intent(in)              :: r
        !> Matrix of size \( \text{ldU} \times r \)
        real(WP), intent(in),  contiguous :: U(:)
        !> Leading dimension of \( U \)
        !>
        !> Possible values: \( \text{ldU} \geq \max(1, m) \)
        integer,  intent(in)              :: ldU
        !> Matrix of size \( \text{ldVT} \times n \)
        real(WP), intent(in),  contiguous :: VT(:)
        !> Leading dimension of \( V^\top \)
        !>
        !> Possible values: \( \text{ldVT} \geq \max(1, r) \)
        integer,  intent(in)              :: ldVT
        !> Factorization rank of \( A B^\top \), i.e. 
        !> number of columns in \( A \) and number of rows in \( B^\top \)
        !>
        !> Possible values: \( k \geq 0 \)
        integer,  intent(in)              :: k
        !> Matrix of size \( \text{ldA} \times k \)
        real(WP), intent(in),  contiguous :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        integer,  intent(in)              :: ldA
        !> Matrix of size \( \text{ldBT} \times n \)
        real(WP), intent(in),  contiguous :: BT(:)
        !> Leading dimension of \( B^\top \)
        !>
        !> Possible values: \( \text{ldBT} \geq \max(1, k) \)
        integer,  intent(in)              :: ldBT
        !> Real workspace of size `lwork`
        real(WP), intent(out), contiguous :: work(:)
        !> Size of the real workspace.
        !>
        !> If `lwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` array is calculated and returned
        !> in `work(1)`.
        integer,  intent(in)              :: lwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)             :: info
        real(WP)                          :: slrdotf
    end function slrdotf

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the Frobenius inner product of two factorized matrices
    module function dlrdotf &
    (m, n, r, U, ldU, VT, ldVT, k, A, ldA, BT, ldBT, work, lwork, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Number of rows in \( U \) and \( A \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,  intent(in)              :: m
        !> Number of columns in \( V^\top \) and \( B^\top \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,  intent(in)              :: n
        !> Factorization rank of \( U V^\top \), i.e. 
        !> number of columns in \( U \) and number of rows in \( V^\top \)
        !>
        !> Possible values: \( r \geq 0 \)
        integer,  intent(in)              :: r
        !> Matrix of size \( \text{ldU} \times r \)
        real(WP), intent(in),  contiguous :: U(:)
        !> Leading dimension of \( U \)
        !>
        !> Possible values: \( \text{ldU} \geq \max(1, m) \)
        integer,  intent(in)              :: ldU
        !> Matrix of size \( \text{ldVT} \times n \)
        real(WP), intent(in),  contiguous :: VT(:)
        !> Leading dimension of \( V^\top \)
        !>
        !> Possible values: \( \text{ldVT} \geq \max(1, r) \)
        integer,  intent(in)              :: ldVT
        !> Factorization rank of \( A B^\top \), i.e. 
        !> number of columns in \( A \) and number of rows in \( B^\top \)
        !>
        !> Possible values: \( k \geq 0 \)
        integer,  intent(in)              :: k
        !> Matrix of size \( \text{ldA} \times k \)
        real(WP), intent(in),  contiguous :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        integer,  intent(in)              :: ldA
        !> Matrix of size \( \text{ldBT} \times n \)
        real(WP), intent(in),  contiguous :: BT(:)
        !> Leading dimension of \( B^\top \)
        !>
        !> Possible values: \( \text{ldBT} \geq \max(1, k) \)
        integer,  intent(in)              :: ldBT
        !> Real workspace of size `lwork`
        real(WP), intent(out), contiguous :: work(:)
        !> Size of the real workspace.
        !>
        !> If `lwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` array is calculated and returned
        !> in `work(1)`.
        integer,  intent(in)              :: lwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)             :: info
        real(WP)                          :: dlrdotf
    end function dlrdotf

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the Frobenius norm of a factorized matrix
    module function slrnrmf &
    (m, n, r, U, ldU, VT, ldVT, work, lwork, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Number of rows in \( U \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,  intent(in)              :: m
        !> Number of columns in \( V^\top \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,  intent(in)              :: n
        !> Factorization rank of \( U V^\top \), i.e. 
        !> number of columns in \( U \) and number of rows in \( V^\top \)
        !>
        !> Possible values: \( r \geq 0 \)
        integer,  intent(in)              :: r
        !> Matrix of size \( \text{ldU} \times r \)
        real(WP), intent(in),  contiguous :: U(:)
        !> Leading dimension of \( U \)
        !>
        !> Possible values: \( \text{ldU} \geq \max(1, m) \)
        integer,  intent(in)              :: ldU
        !> Matrix of size \( \text{ldVT} \times n \)
        real(WP), intent(in),  contiguous :: VT(:)
        !> Leading dimension of \( V^\top \)
        !>
        !> Possible values: \( \text{ldVT} \geq \max(1, r) \)
        integer,  intent(in)              :: ldVT
        !> Real workspace of size `lwork`
        real(WP), intent(out), contiguous :: work(:)
        !> Size of the real workspace.
        !>
        !> If `lwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` array is calculated and returned
        !> in `work(1)`.
        integer,  intent(in)              :: lwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)             :: info
        real(WP)                          :: slrnrmf
    end function slrnrmf

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the Frobenius norm of a factorized matrix
    module function dlrnrmf &
    (m, n, r, U, ldU, VT, ldVT, work, lwork, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Number of rows in \( U \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,  intent(in)              :: m
        !> Number of columns in \( V^\top \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,  intent(in)              :: n
        !> Factorization rank of \( U V^\top \), i.e. 
        !> number of columns in \( U \) and number of rows in \( V^\top \)
        !>
        !> Possible values: \( r \geq 0 \)
        integer,  intent(in)              :: r
        !> Matrix of size \( \text{ldU} \times r \)
        real(WP), intent(in),  contiguous :: U(:)
        !> Leading dimension of \( U \)
        !>
        !> Possible values: \( \text{ldU} \geq \max(1, m) \)
        integer,  intent(in)              :: ldU
        !> Matrix of size \( \text{ldVT} \times n \)
        real(WP), intent(in),  contiguous :: VT(:)
        !> Leading dimension of \( V^\top \)
        !>
        !> Possible values: \( \text{ldVT} \geq \max(1, r) \)
        integer,  intent(in)              :: ldVT
        !> Real workspace of size `lwork`
        real(WP), intent(out), contiguous :: work(:)
        !> Size of the real workspace.
        !>
        !> If `lwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` array is calculated and returned
        !> in `work(1)`.
        integer,  intent(in)              :: lwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)             :: info
        real(WP)                          :: dlrnrmf
    end function dlrnrmf

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the Frobenius norm of a difference of two factorized matrices
    module function slrnrmf_diff &
    (m, n, r, U, ldU, VT, ldVT, k, A, ldA, BT, ldBT, work, lwork, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Number of rows in \( U \) and \( A \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,  intent(in)              :: m
        !> Number of columns in \( V^\top \) and \( B^\top \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,  intent(in)              :: n
        !> Factorization rank of \( U V^\top \), i.e. 
        !> number of columns in \( U \) and number of rows in \( V^\top \)
        !>
        !> Possible values: \( r \geq 0 \)
        integer,  intent(in)              :: r
        !> Matrix of size \( \text{ldU} \times r \)
        real(WP), intent(in),  contiguous :: U(:)
        !> Leading dimension of \( U \)
        !>
        !> Possible values: \( \text{ldU} \geq \max(1, m) \)
        integer,  intent(in)              :: ldU
        !> Matrix of size \( \text{ldVT} \times n \)
        real(WP), intent(in),  contiguous :: VT(:)
        !> Leading dimension of \( V^\top \)
        !>
        !> Possible values: \( \text{ldVT} \geq \max(1, r) \)
        integer,  intent(in)              :: ldVT
        !> Factorization rank of \( A B^\top \), i.e. 
        !> number of columns in \( A \) and number of rows in \( B^\top \)
        !>
        !> Possible values: \( k \geq 0 \)
        integer,  intent(in)              :: k
        !> Matrix of size \( \text{ldA} \times k \)
        real(WP), intent(in),  contiguous :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        integer,  intent(in)              :: ldA
        !> Matrix of size \( \text{ldBT} \times n \)
        real(WP), intent(in),  contiguous :: BT(:)
        !> Leading dimension of \( B^\top \)
        !>
        !> Possible values: \( \text{ldBT} \geq \max(1, k) \)
        integer,  intent(in)              :: ldBT
        !> Real workspace of size `lwork`
        real(WP), intent(out), contiguous :: work(:)
        !> Size of the real workspace.
        !>
        !> If `lwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` array is calculated and returned
        !> in `work(1)`.
        integer,  intent(in)              :: lwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)             :: info
        real(WP)                          :: slrnrmf_diff
    end function slrnrmf_diff

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the Frobenius norm of a difference of two factorized matrices
    module function dlrnrmf_diff &
    (m, n, r, U, ldU, VT, ldVT, k, A, ldA, BT, ldBT, work, lwork, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Number of rows in \( U \) and \( A \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,  intent(in)              :: m
        !> Number of columns in \( V^\top \) and \( B^\top \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,  intent(in)              :: n
        !> Factorization rank of \( U V^\top \), i.e. 
        !> number of columns in \( U \) and number of rows in \( V^\top \)
        !>
        !> Possible values: \( r \geq 0 \)
        integer,  intent(in)              :: r
        !> Matrix of size \( \text{ldU} \times r \)
        real(WP), intent(in),  contiguous :: U(:)
        !> Leading dimension of \( U \)
        !>
        !> Possible values: \( \text{ldU} \geq \max(1, m) \)
        integer,  intent(in)              :: ldU
        !> Matrix of size \( \text{ldVT} \times n \)
        real(WP), intent(in),  contiguous :: VT(:)
        !> Leading dimension of \( V^\top \)
        !>
        !> Possible values: \( \text{ldVT} \geq \max(1, r) \)
        integer,  intent(in)              :: ldVT
        !> Factorization rank of \( A B^\top \), i.e. 
        !> number of columns in \( A \) and number of rows in \( B^\top \)
        !>
        !> Possible values: \( k \geq 0 \)
        integer,  intent(in)              :: k
        !> Matrix of size \( \text{ldA} \times k \)
        real(WP), intent(in),  contiguous :: A(:)
        !> Leading dimension of \( A \)
        !>
        !> Possible values: \( \text{ldA} \geq \max(1, m) \)
        integer,  intent(in)              :: ldA
        !> Matrix of size \( \text{ldBT} \times n \)
        real(WP), intent(in),  contiguous :: BT(:)
        !> Leading dimension of \( B^\top \)
        !>
        !> Possible values: \( \text{ldBT} \geq \max(1, k) \)
        integer,  intent(in)              :: ldBT
        !> Real workspace of size `lwork`
        real(WP), intent(out), contiguous :: work(:)
        !> Size of the real workspace.
        !>
        !> If `lwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` array is calculated and returned
        !> in `work(1)`.
        integer,  intent(in)              :: lwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)             :: info
        real(WP)                          :: dlrnrmf_diff
    end function dlrnrmf_diff

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine lrort_rank &
    (side, m, n, r, newr, info)
    implicit none
        !> Specifies whether to compute the left- or right orthogonal representation
        !>
        !> - side='L': left-orthogonal
        !> - side='R': right-orthogonal
        character(1), intent(in)  :: side
        !> Number of rows in \( U \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,      intent(in)  :: m
        !> Number of columns in \( V^\top \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,      intent(in)  :: n
        !> Factorization rank of \( U V^\top \), i.e. 
        !> number of columns in \( U \) and number of rows in \( V^\top \)
        !>
        !> Possible values: \( r \geq 0 \)
        integer,      intent(in)  :: r
        !> New factorization rank
        integer,      intent(out) :: newr
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,      intent(out) :: info
    end subroutine lrort_rank
       
    !------------------------------------------------------------------------------------------------------------------------

    !> Computes an implicit left-orthogonal or right-orthogonal representation of a factorized matrix.
    !> The resulting factorization is of rank \( \min(m, r) \) or \( \min(n, r) \)
    module subroutine slrort &
    (side, m, n, r, U, ldu, VT, ldvt, tau, work, lwork, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Specifies whether to compute the left- or right orthogonal representation
        !>
        !> - side='L': left-orthogonal
        !> - side='R': right-orthogonal
        character(1), intent(in)                :: side
        !> Number of rows in \( U \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,      intent(in)                :: m
        !> Number of columns in \( V^\top \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,      intent(in)                :: n
        !> Factorization rank of \( U V^\top \), i.e. 
        !> number of columns in \( U \) and number of rows in \( V^\top \)
        !>
        !> Possible values: \( r \geq 0 \)
        integer,      intent(in)                :: r
        !> Matrix of size \( \text{ldU} \times r \)
        !>
        !> If side='L", on exit, contains elementay reflectors
        real(WP),     intent(inout), contiguous :: U(:)
        !> Leading dimension of \( U \)
        !>
        !> Possible values: \( \text{ldU} \geq \max(1, m) \)
        integer,      intent(in)                :: ldU
        !> Matrix of size \( \text{ldVT} \times n \)
        !>
        !> if side='R', on exit, contains elementary reflectors
        real(WP),     intent(inout), contiguous :: VT(:)
        !> Leading dimension of \( V^\top \)
        !>
        !> Possible values: \( \text{ldVT} \geq \max(1, r) \)
        integer,      intent(in)                :: ldVT
        !> Scalar factors of elementary reflectors
        !>
        !> - side='L': \( m, r \)
        !> - side='R': \( n, r \)
        real(WP),     intent(out),   contiguous :: tau(:)
        !> Real workspace of size `lwork`
        real(WP),     intent(out),   contiguous :: work(:)
        !> Size of the real workspace.
        !>
        !> If `lwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` array is calculated and returned
        !> in `work(1)`.
        integer,      intent(in)                :: lwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,      intent(out)               :: info
    end subroutine slrort

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes an implicit left-orthogonal or right-orthogonal representation of a factorized matrix.
    !> The resulting factorization is of rank \( \min(m, r) \) or \( \min(n, r) \)
    module subroutine dlrort &
    (side, m, n, r, U, ldu, VT, ldvt, tau, work, lwork, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Specifies whether to compute the left- or right orthogonal representation
        !>
        !> - side='L': left-orthogonal
        !> - side='R': right-orthogonal
        character(1), intent(in)                :: side
        !> Number of rows in \( U \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,      intent(in)                :: m
        !> Number of columns in \( V^\top \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,      intent(in)                :: n
        !> Factorization rank of \( U V^\top \), i.e. 
        !> number of columns in \( U \) and number of rows in \( V^\top \)
        !>
        !> Possible values: \( r \geq 0 \)
        integer,      intent(in)                :: r
        !> Matrix of size \( \text{ldU} \times r \)
        !>
        !> If side='L", on exit, contains elementay reflectors
        real(WP),     intent(inout), contiguous :: U(:)
        !> Leading dimension of \( U \)
        !>
        !> Possible values: \( \text{ldU} \geq \max(1, m) \)
        integer,      intent(in)                :: ldU
        !> Matrix of size \( \text{ldVT} \times n \)
        !>
        !> if side='R', on exit, contains elementary reflectors
        real(WP),     intent(inout), contiguous :: VT(:)
        !> Leading dimension of \( V^\top \)
        !>
        !> Possible values: \( \text{ldVT} \geq \max(1, r) \)
        integer,      intent(in)                :: ldVT
        !> Scalar factors of elementary reflectors
        !>
        !> - side='L': \( m, r \)
        !> - side='R': \( n, r \)
        real(WP),     intent(out),   contiguous :: tau(:)
        !> Real workspace of size `lwork`
        real(WP),     intent(out),   contiguous :: work(:)
        !> Size of the real workspace.
        !>
        !> If `lwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` array is calculated and returned
        !> in `work(1)`.
        integer,      intent(in)                :: lwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,      intent(out)               :: info
    end subroutine dlrort

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the singular value decomposition of a left- or right-orthogonal factorized matrix
    module subroutine slrsvd_ort &
    (job, side, m, n, r, U, ldu, VT, ldvt, tau, S, Q, ldQ, PT, ldPT, work, lwork, iwork, liwork, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Specifies whether to compute the singular vectors or not:
        !>
        !> - job='S': computes left and right singular vectors 
        !> and stores them in \( Q \) and \( P^\top \)
        !> - job='N': computes only the singular values
        character(1), intent(in)                :: job
        !> Specifies whether the left- or right orthogonal representation is given
        !>
        !> - side='L': left-orthogonal
        !> - side='R': right-orthogonal
        character(1), intent(in)                :: side
        !> Number of rows in \( U \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,      intent(in)                :: m
        !> Number of columns in \( V^\top \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,      intent(in)                :: n
        !> Factorization rank of \( U V^\top \), i.e. 
        !> number of columns in \( U \) and number of rows in \( V^\top \)
        !>
        !> Possible values: \( 0 \leq r \leq \min(m,n) \)
        integer,      intent(in)                :: r
        !> Matrix of size \( \text{ldU} \times r \)
        !>
        !> Contents are destroyed
        real(WP),     intent(inout), contiguous :: U(:)
        !> Leading dimension of \( U \)
        !>
        !> Possible values: \( \text{ldU} \geq \max(1, m) \)
        integer,      intent(in)                :: ldU
        !> Matrix of size \( \text{ldVT} \times n \)
        !>
        !> Contents are destroyed
        real(WP),     intent(inout), contiguous :: VT(:)
        !> Leading dimension of \( V^\top \)
        !>
        !> Possible values: \( \text{ldVT} \geq \max(1, r) \)
        integer,      intent(in)                :: ldVT
        !> Scalar factors of elementary reflectors
        !>
        !> - side='L': \( m, r \)
        !> - side='R': \( n, r \)
        real(WP),     intent(in),    contiguous :: tau(:)
        !> Array of size \( \min\(m, n, r\) \)
        !>
        !> On exit, stores singular values
        real(WP),     intent(out),   contiguous :: S(:)
        !> Matrix of size \( \text{ldU} \times r \)
        !>
        !> If `job="S"`, stores the left singular vectors.
        !> Is not referenced otherwise
        real(WP),     intent(out),   contiguous :: Q(:)
        !> Leading dimension of \( Q \)
        !>
        !> If `job="S"`, \( \text{ldQ} \geq \max(1,m) \)
        !> Else, \( \text{ldQ} \geq 1 \)
        integer,      intent(in)                :: ldQ
        !> Matrix of size \( \text{ldPT} \times n \)
        !>
        !> If `job="S"`, stores the right singular vectors.
        !> Is not referenced otherwise
        real(WP),     intent(out),   contiguous :: PT(:)
        !> Leading dimension of \( P^\top \)
        !>
        !> If `job="S"`, \( \text{ldP} \geq \max(1,r) \)
        !> Else, \( \text{ldPT} \geq 1 \)
        integer,      intent(in)                :: ldPT
        !> Real workspace of size `lwork`
        real(WP),     intent(out),   contiguous :: work(:)
        !> Size of the real workspace.
        !>
        !> If `lwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` and `iwork` arrays is calculated and returned
        !> in `work(1)` and `iwork(1)`, respectively.
        integer,      intent(in)                :: lwork
        !> Integer workspace of size `liwork`
        integer,      intent(out),   contiguous :: iwork(:)
        !> Size of the integer workspace
        !>
        !> If `liwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` and `iwork` arrays is calculated and returned
        !> in `work(1)` and `iwork(1)`, respectively.
        integer,      intent(in)                :: liwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        !> - `info > 0`: SVD didn't converge
        integer,      intent(out)               :: info
    end subroutine slrsvd_ort

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the singular value decomposition of a left- or right-orthogonal factorized matrix
    module subroutine dlrsvd_ort &
    (job, side, m, n, r, U, ldu, VT, ldvt, tau, S, Q, ldQ, PT, ldPT, work, lwork, iwork, liwork, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Specifies whether to compute the singular vectors or not:
        !>
        !> - job='S': computes left and right singular vectors 
        !> and stores them in \( Q \) and \( P^\top \)
        !> - job='N': computes only the singular values
        character(1), intent(in)                :: job
        !> Specifies whether the left- or right orthogonal representation is given
        !>
        !> - side='L': left-orthogonal
        !> - side='R': right-orthogonal
        character(1), intent(in)                :: side
        !> Number of rows in \( U \)
        !>
        !> Possible values: \( m \geq 0 \)
        integer,      intent(in)                :: m
        !> Number of columns in \( V^\top \)
        !>
        !> Possible values: \( n \geq 0 \)
        integer,      intent(in)                :: n
        !> Factorization rank of \( U V^\top \), i.e. 
        !> number of columns in \( U \) and number of rows in \( V^\top \)
        !>
        !> Possible values: \( 0 \leq r \leq \min(m,n) \)
        integer,      intent(in)                :: r
        !> Matrix of size \( \text{ldU} \times r \)
        !>
        !> Contents are destroyed
        real(WP),     intent(inout), contiguous :: U(:)
        !> Leading dimension of \( U \)
        !>
        !> Possible values: \( \text{ldU} \geq \max(1, m) \)
        integer,      intent(in)                :: ldU
        !> Matrix of size \( \text{ldVT} \times n \)
        !>
        !> Contents are destroyed
        real(WP),     intent(inout), contiguous :: VT(:)
        !> Leading dimension of \( V^\top \)
        !>
        !> Possible values: \( \text{ldVT} \geq \max(1, r) \)
        integer,      intent(in)                :: ldVT
        !> Scalar factors of elementary reflectors
        !>
        !> - side='L': \( m, r \)
        !> - side='R': \( n, r \)
        real(WP),     intent(in),    contiguous :: tau(:)
        !> Array of size \( \min\(m, n, r\) \)
        !>
        !> On exit, stores singular values
        real(WP),     intent(out),   contiguous :: S(:)
        !> Matrix of size \( \text{ldU} \times r \)
        !>
        !> If `job="S"`, stores the left singular vectors.
        !> Is not referenced otherwise
        real(WP),     intent(out),   contiguous :: Q(:)
        !> Leading dimension of \( Q \)
        !>
        !> If `job="S"`, \( \text{ldQ} \geq \max(1,m) \)
        !> Else, \( \text{ldQ} \geq 1 \)
        integer,      intent(in)                :: ldQ
        !> Matrix of size \( \text{ldPT} \times n \)
        !>
        !> If `job="S"`, stores the right singular vectors.
        !> Is not referenced otherwise
        real(WP),     intent(out),   contiguous :: PT(:)
        !> Leading dimension of \( P^\top \)
        !>
        !> If `job="S"`, \( \text{ldP} \geq \max(1,r) \)
        !> Else, \( \text{ldPT} \geq 1 \)
        integer,      intent(in)                :: ldPT
        !> Real workspace of size `lwork`
        real(WP),     intent(out),   contiguous :: work(:)
        !> Size of the real workspace.
        !>
        !> If `lwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` and `iwork` arrays is calculated and returned
        !> in `work(1)` and `iwork(1)`, respectively.
        integer,      intent(in)                :: lwork
        !> Integer workspace of size `liwork`
        integer,      intent(out),   contiguous :: iwork(:)
        !> Size of the integer workspace
        !>
        !> If `liwork = -1`, a workspace query is assumed.
        !> The optimal size for the `work` and `iwork` arrays is calculated and returned
        !> in `work(1)` and `iwork(1)`, respectively.
        integer,      intent(in)                :: liwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        !> - `info > 0`: SVD didn't converge
        integer,      intent(out)               :: info
    end subroutine dlrsvd_ort

    !------------------------------------------------------------------------------------------------------------------------
end interface
end module maria_lr_la_mod
