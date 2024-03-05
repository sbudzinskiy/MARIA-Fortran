!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_lr_cross_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Cross/pseudoskeleton approximations of matrices.
!----------------------------------------------------------------------------------------------------------------------------
module maria_lr_cross_mod
implicit none (type, external)

! Procedures
public :: smatcross_top, &
          dmatcross_top, &
          smatcross,     &
          dmatcross,     &
          smatcross_aca, &
          dmatcross_aca

interface
    !------------------------------------------------------------------------------------------------------------------------
    
    !> Computes cross approximation from top rows
    module subroutine smatcross_top &
    (m, n, nc, cols, ldc, nr, rows, ldr, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,   only: &
        WP => SP
    implicit none
        integer,       intent(in)                :: m
        integer,       intent(in)                :: n
        integer,       intent(in)                :: nc
        real(WP),      intent(inout), contiguous :: cols(:)
        integer,       intent(in)                :: ldc
        integer,       intent(in)                :: nr
        real(WP),      intent(inout), contiguous :: rows(:)
        integer,       intent(in)                :: ldr
        integer,       intent(in)                :: r
        real(WP),      intent(out),   contiguous :: U(:)
        integer,       intent(in)                :: ldu
        real(WP),      intent(out),   contiguous :: VT(:)
        integer,       intent(in)                :: ldvt
        real(WP),      intent(out),   contiguous :: work(:)
        integer,       intent(in)                :: lwork
        integer,       intent(out),   contiguous :: iwork(:)
        integer,       intent(in)                :: liwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,       intent(out)               :: info
    end subroutine smatcross_top

    !------------------------------------------------------------------------------------------------------------------------
    
    !> Computes cross approximation from top rows
    module subroutine dmatcross_top &
    (m, n, nc, cols, ldc, nr, rows, ldr, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,   only: &
        WP => DP
    implicit none
        integer,       intent(in)                :: m
        integer,       intent(in)                :: n
        integer,       intent(in)                :: nc
        real(WP),      intent(inout), contiguous :: cols(:)
        integer,       intent(in)                :: ldc
        integer,       intent(in)                :: nr
        real(WP),      intent(inout), contiguous :: rows(:)
        integer,       intent(in)                :: ldr
        integer,       intent(in)                :: r
        real(WP),      intent(out),   contiguous :: U(:)
        integer,       intent(in)                :: ldu
        real(WP),      intent(out),   contiguous :: VT(:)
        integer,       intent(in)                :: ldvt
        real(WP),      intent(out),   contiguous :: work(:)
        integer,       intent(in)                :: lwork
        integer,       intent(out),   contiguous :: iwork(:)
        integer,       intent(in)                :: liwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,       intent(out)               :: info
    end subroutine dmatcross_top

    !------------------------------------------------------------------------------------------------------------------------
    
    !> Computes cross approximation from rows and columns
    module subroutine smatcross &
    (m, n, fun, nc, icol, nr, irow, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,   only: &
        WP => SP
    use maria_access_matrix_mod, only: &
        MS => smatslc
    implicit none
        integer,       intent(in)                :: m
        integer,       intent(in)                :: n
        procedure(MS), intent(in),    pointer    :: fun
        integer,       intent(in)                :: nc
        integer,       intent(inout), contiguous :: icol(:)
        integer,       intent(in)                :: nr
        integer,       intent(inout), contiguous :: irow(:)
        integer,       intent(in)                :: r
        real(WP),      intent(out),   contiguous :: U(:)
        integer,       intent(in)                :: ldu
        real(WP),      intent(out),   contiguous :: VT(:)
        integer,       intent(in)                :: ldvt
        real(WP),      intent(out),   contiguous :: work(:)
        integer,       intent(in)                :: lwork
        integer,       intent(out),   contiguous :: iwork(:)
        integer,       intent(in)                :: liwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,       intent(out)               :: info
    end subroutine smatcross

    !------------------------------------------------------------------------------------------------------------------------
    
    !> Computes cross approximation from rows and columns
    module subroutine dmatcross &
    (m, n, fun, nc, icol, nr, irow, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,   only: &
        WP => DP
    use maria_access_matrix_mod, only: &
        MS => dmatslc
    implicit none
        integer,       intent(in)                :: m
        integer,       intent(in)                :: n
        procedure(MS), intent(in),    pointer    :: fun
        integer,       intent(in)                :: nc
        integer,       intent(inout), contiguous :: icol(:)
        integer,       intent(in)                :: nr
        integer,       intent(inout), contiguous :: irow(:)
        integer,       intent(in)                :: r
        real(WP),      intent(out),   contiguous :: U(:)
        integer,       intent(in)                :: ldu
        real(WP),      intent(out),   contiguous :: VT(:)
        integer,       intent(in)                :: ldvt
        real(WP),      intent(out),   contiguous :: work(:)
        integer,       intent(in)                :: lwork
        integer,       intent(out),   contiguous :: iwork(:)
        integer,       intent(in)                :: liwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,       intent(out)               :: info
    end subroutine dmatcross

    !------------------------------------------------------------------------------------------------------------------------
    
    !> Computes adaptive cross approximation with rook pivoting
    module subroutine smatcross_aca &
    (m, n, matval, col_order, maxr, niter_rook, r, U, ldu, VT, ldvt, irow, icol, work, lwork, info, rtol)
    use maria_kinds_mod,   only: &
        WP => SP
    use maria_access_matrix_mod, only: &
        MV => smatval
    implicit none
        integer,       intent(in)                :: m
        integer,       intent(in)                :: n
        procedure(MV), intent(in),    pointer    :: matval
        integer,       intent(in),    contiguous :: col_order(:)
        integer,       intent(in)                :: maxr
        integer,       intent(in)                :: niter_rook
        integer,       intent(out)               :: r
        real(WP),      intent(out),   contiguous :: U(:)
        integer,       intent(in)                :: ldu
        real(WP),      intent(out),   contiguous :: VT(:)
        integer,       intent(in)                :: ldvt
        integer,       intent(out),   contiguous :: irow(:)
        integer,       intent(out),   contiguous :: icol(:)
        real(WP),      intent(out),   contiguous :: work(:)
        integer,       intent(in)                :: lwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,       intent(out)               :: info(2)
        real(WP),      intent(in),   optional    :: rtol
    end subroutine smatcross_aca

    !------------------------------------------------------------------------------------------------------------------------
    
    !> Computes adaptive cross approximation with rook pivoting
    module subroutine dmatcross_aca &
    (m, n, matval, col_order, maxr, niter_rook, r, U, ldu, VT, ldvt, irow, icol, work, lwork, info, rtol)
    use maria_kinds_mod,   only: &
        WP => DP
    use maria_access_matrix_mod, only: &
        MV => dmatval
    implicit none
        integer,       intent(in)                :: m
        integer,       intent(in)                :: n
        procedure(MV), intent(in),    pointer    :: matval
        integer,       intent(in),    contiguous :: col_order(:)
        integer,       intent(in)                :: maxr
        integer,       intent(in)                :: niter_rook
        integer,       intent(out)               :: r
        real(WP),      intent(out),   contiguous :: U(:)
        integer,       intent(in)                :: ldu
        real(WP),      intent(out),   contiguous :: VT(:)
        integer,       intent(in)                :: ldvt
        integer,       intent(out),   contiguous :: irow(:)
        integer,       intent(out),   contiguous :: icol(:)
        real(WP),      intent(out),   contiguous :: work(:)
        integer,       intent(in)                :: lwork
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,       intent(out)               :: info(2)
        real(WP),      intent(in),   optional    :: rtol
    end subroutine dmatcross_aca

    !------------------------------------------------------------------------------------------------------------------------
end interface
end module maria_lr_cross_mod
