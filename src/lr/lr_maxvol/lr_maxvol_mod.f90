!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_lr_maxvol_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Routines for computing the maximum volume submatrices.
!----------------------------------------------------------------------------------------------------------------------------
module maria_lr_maxvol_mod
implicit none (type, external)

! Procedures
public :: sgevolume,    &
          dgevolume,    &
          sgemaxvol_swap_rows, &
          dgemaxvol_swap_rows, &
          sgemaxvol,    &
          dgemaxvol,    &
          sgemaxvol_rect_add_rows, &
          dgemaxvol_rect_add_rows, &
          sgemaxvol_rect_swap_rows, &
          dgemaxvol_rect_swap_rows, &
          sgemaxvol_rect_swap_cols, &
          dgemaxvol_rect_swap_cols, &
          sgemaxvol_proj,           &
          dgemaxvol_proj

interface
    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the volume of a matrix
    module function sgevolume &
    (m, n, r, A, lda, work, lwork, iwork, liwork, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,  intent(in)                :: m
        integer,  intent(in)                :: n
        integer,  intent(in)                :: r
        real(WP), intent(inout), contiguous :: A(:)
        integer,  intent(in)                :: lda
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: info
        real(WP)                            :: sgevolume
    end function sgevolume

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the volume of a matrix
    module function dgevolume &
    (m, n, r, A, lda, work, lwork, iwork, liwork, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,  intent(in)                :: m
        integer,  intent(in)                :: n
        integer,  intent(in)                :: r
        real(WP), intent(inout), contiguous :: A(:)
        integer,  intent(in)                :: lda
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: info
        real(WP)                            :: dgevolume
    end function dgevolume

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes a dominant square submatrix 
    module subroutine sgemaxvol_swap_rows &
    (ort, m, n, A, lda, irow, thresh, work, lwork, iwork, liwork, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        character(1), intent(in)                :: ort
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        real(WP), intent(inout), contiguous :: A(:)
        integer,  intent(in)                :: lda
        integer,  intent(inout), contiguous :: irow(:)
        real(WP), intent(in)                :: thresh
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: info
    end subroutine sgemaxvol_swap_rows

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes a dominant square submatrix 
    module subroutine dgemaxvol_swap_rows &
    (ort, m, n, A, lda, irow, thresh, work, lwork, iwork, liwork, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        character(1), intent(in)                :: ort
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        real(WP), intent(inout), contiguous :: A(:)
        integer,  intent(in)                :: lda
        integer,  intent(inout), contiguous :: irow(:)
        real(WP), intent(in)                :: thresh
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: info
    end subroutine dgemaxvol_swap_rows

    !------------------------------------------------------------------------------------------------------------------------

    !> Finds the position of matrix of maximum volume
    module subroutine sgemaxvol &
    (m, n, matslc, r, irow, icol, niter, thresh, work, lwork, iwork, liwork, info)
    use maria_kinds_mod, only: &
        WP => SP
    use maria_access_matrix_mod, only: &
        MS => smatslc
    implicit none
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        procedure(MS), intent(in), pointer :: matslc
        integer,  intent(in)                :: r
        integer,  intent(inout), contiguous :: irow(:)
        integer,  intent(inout), contiguous :: icol(:)
        integer,  intent(in)                :: niter
        real(WP), intent(in)                :: thresh
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: info
    end subroutine sgemaxvol

    !------------------------------------------------------------------------------------------------------------------------

    !> Finds the position of matrix of maximum volume
    module subroutine dgemaxvol &
    (m, n, matslc, r, irow, icol, niter, thresh, work, lwork, iwork, liwork, info)
    use maria_kinds_mod, only: &
        WP => DP
    use maria_access_matrix_mod, only: &
        MS => dmatslc
    implicit none
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        procedure(MS), intent(in), pointer :: matslc
        integer,  intent(in)                :: r
        integer,  intent(inout), contiguous :: irow(:)
        integer,  intent(inout), contiguous :: icol(:)
        integer,  intent(in)                :: niter
        real(WP), intent(in)                :: thresh
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: info
    end subroutine dgemaxvol

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes a dominant rectangular submatrix 
    module subroutine sgemaxvol_rect_add_rows &
    (ort, m, n, A, lda, k, irow, work, lwork, iwork, liwork, ierr)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        character(1), intent(in)                :: ort
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        real(WP), intent(inout), contiguous :: A(:)
        integer,  intent(in)                :: lda
        integer,  intent(in)                :: k
        integer,  intent(inout), contiguous :: irow(:)
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: ierr
    end subroutine sgemaxvol_rect_add_rows

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes a dominant rectangular submatrix 
    module subroutine dgemaxvol_rect_add_rows &
    (ort, m, n, A, lda, k, irow, work, lwork, iwork, liwork, ierr)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        character(1), intent(in)                :: ort
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        real(WP), intent(inout), contiguous :: A(:)
        integer,  intent(in)                :: lda
        integer,  intent(in)                :: k
        integer,  intent(inout), contiguous :: irow(:)
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: ierr
    end subroutine dgemaxvol_rect_add_rows

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes a dominant square submatrix 
    module subroutine sgemaxvol_rect_swap_rows &
    (ort, m, n, A, lda, k, irow, thresh, work, lwork, iwork, liwork, ierr)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        character(1), intent(in)                :: ort
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        real(WP), intent(inout), contiguous :: A(:)
        integer,  intent(in)                :: lda
        integer,      intent(in)                :: k
        integer,  intent(inout), contiguous :: irow(:)
        real(WP), intent(in)                :: thresh
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: ierr
    end subroutine sgemaxvol_rect_swap_rows

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes a dominant square submatrix 
    module subroutine dgemaxvol_rect_swap_rows &
    (ort, m, n, A, lda, k, irow, thresh, work, lwork, iwork, liwork, ierr)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        character(1), intent(in)                :: ort
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        real(WP), intent(inout), contiguous :: A(:)
        integer,  intent(in)                :: lda
        integer,      intent(in)                :: k
        integer,  intent(inout), contiguous :: irow(:)
        real(WP), intent(in)                :: thresh
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: ierr
    end subroutine dgemaxvol_rect_swap_rows

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes a dominant square submatrix 
    module subroutine sgemaxvol_rect_swap_cols &
    (m, n, A, lda, k, icol, thresh, work, lwork, iwork, liwork, ierr)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        real(WP), intent(inout), contiguous :: A(:)
        integer,  intent(in)                :: lda
        integer,      intent(in)                :: k
        integer,  intent(inout), contiguous :: icol(:)
        real(WP), intent(in)                :: thresh
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: ierr
    end subroutine sgemaxvol_rect_swap_cols

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes a dominant square submatrix 
    module subroutine dgemaxvol_rect_swap_cols &
    (m, n, A, lda, k, icol, thresh, work, lwork, iwork, liwork, ierr)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        real(WP), intent(inout), contiguous :: A(:)
        integer,  intent(in)                :: lda
        integer,      intent(in)                :: k
        integer,  intent(inout), contiguous :: icol(:)
        real(WP), intent(in)                :: thresh
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: ierr
    end subroutine dgemaxvol_rect_swap_cols

    !------------------------------------------------------------------------------------------------------------------------

    !> Finds the position of matrix of maximum volume
    module subroutine sgemaxvol_proj &
    (m, n, matslc, r, kr, irow, icol_short, kc, icol, irow_short, niter, thresh, work, lwork, iwork, liwork, ierr)
    use maria_kinds_mod, only: &
        WP => SP
    use maria_access_matrix_mod, only: &
        MS => smatslc
    implicit none
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        procedure(MS), intent(in), pointer :: matslc
        integer,  intent(in)                :: r
        integer,  intent(in)        :: kr
        integer,  intent(inout), contiguous :: irow(:)
        integer, intent(inout), contiguous :: icol_short(:)
        integer, intent(in)         :: kc
        integer,  intent(inout), contiguous :: icol(:)
        integer, intent(inout), contiguous :: irow_short(:)
        integer,  intent(in)                :: niter
        real(WP), intent(in)                :: thresh
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: ierr
    end subroutine sgemaxvol_proj

    !------------------------------------------------------------------------------------------------------------------------

    !> Finds the position of matrix of maximum volume
    module subroutine dgemaxvol_proj &
    (m, n, matslc, r, kr, irow, icol_short, kc, icol, irow_short, niter, thresh, work, lwork, iwork, liwork, ierr)
    use maria_kinds_mod, only: &
        WP => DP
    use maria_access_matrix_mod, only: &
        MS => dmatslc
    implicit none
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        procedure(MS), intent(in), pointer :: matslc
        integer,  intent(in)                :: r
        integer,  intent(in)        :: kr
        integer,  intent(inout), contiguous :: irow(:)
        integer, intent(inout), contiguous :: icol_short(:)
        integer, intent(in)         :: kc
        integer,  intent(inout), contiguous :: icol(:)
        integer, intent(inout), contiguous :: irow_short(:)
        integer,  intent(in)                :: niter
        real(WP), intent(in)                :: thresh
        real(WP), intent(out),   contiguous :: work(:)
        integer,  intent(in)                :: lwork
        integer,  intent(out),   contiguous :: iwork(:)
        integer,  intent(in)                :: liwork
        integer,  intent(out)               :: ierr
    end subroutine dgemaxvol_proj

    !------------------------------------------------------------------------------------------------------------------------
end interface
end module maria_lr_maxvol_mod
