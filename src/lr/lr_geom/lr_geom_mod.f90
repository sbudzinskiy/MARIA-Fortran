!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_lr_geom_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Routines for computing Riemannian things
!----------------------------------------------------------------------------------------------------------------------------
module maria_lr_geom_mod
implicit none (type, external)

! Procedures
public :: slrproj_tangent, &
          dlrproj_tangent, &
          slrdotf_tangent, &
          dlrdotf_tangent, &
          slrretr_tangent, &
          dlrretr_tangent

interface
    !------------------------------------------------------------------------------------------------------------------------

    !> Project onto fixed-rank tangent space
    module subroutine slrproj_tangent &
    (m, n, r, U, ldu, VT, ldvt, B2AB, B2BA, pU, ldpu, pVT, ldpvt, UTAV, ldutav, ierr)
    use maria_kinds_mod,   only: &
        WP => SP
    use maria_la_core_mod, only: &
        MM => smatmul
    implicit none
        integer,        intent(in   )               :: m
        integer,        intent(in   )               :: n
        integer,        intent(in   )               :: r
        real(WP),       intent(in   ),  contiguous  :: U(:)
        integer,        intent(in   )               :: ldu
        real(WP),       intent(in   ),  contiguous  :: VT(:)
        integer,        intent(in   )               :: ldvt
        procedure(MM),  intent(in   ),  pointer     :: B2AB
        procedure(MM),  intent(in   ),  pointer     :: B2BA
        real(WP),       intent(  out),  contiguous  :: pU(:)
        integer,        intent(in   )               :: ldpu
        real(WP),       intent(  out),  contiguous  :: pVT(:)
        integer,        intent(in   )               :: ldpvt
        real(WP),       intent(  out),  contiguous  :: UTAV(:)
        integer,        intent(in   )               :: ldutav
        integer,        intent(  out)               :: ierr
    end subroutine slrproj_tangent

    !------------------------------------------------------------------------------------------------------------------------

    !> Project onto fixed-rank tangent space
    module subroutine dlrproj_tangent &
    (m, n, r, U, ldu, VT, ldvt, B2AB, B2BA, pU, ldpu, pVT, ldpvt, UTAV, ldutav, ierr)
    use maria_kinds_mod,   only: &
        WP => DP
    use maria_la_core_mod, only: &
        MM => dmatmul
    implicit none
        integer,        intent(in   )               :: m
        integer,        intent(in   )               :: n
        integer,        intent(in   )               :: r
        real(WP),       intent(in   ),  contiguous  :: U(:)
        integer,        intent(in   )               :: ldu
        real(WP),       intent(in   ),  contiguous  :: VT(:)
        integer,        intent(in   )               :: ldvt
        procedure(MM),  intent(in   ),  pointer     :: B2AB
        procedure(MM),  intent(in   ),  pointer     :: B2BA
        real(WP),       intent(  out),  contiguous  :: pU(:)
        integer,        intent(in   )               :: ldpu
        real(WP),       intent(  out),  contiguous  :: pVT(:)
        integer,        intent(in   )               :: ldpvt
        real(WP),       intent(  out),  contiguous  :: UTAV(:)
        integer,        intent(in   )               :: ldutav
        integer,        intent(  out)               :: ierr
    end subroutine dlrproj_tangent

    !------------------------------------------------------------------------------------------------------------------------

    !> Dot product of tangent vectors
    module function slrdotf_tangent &
    (m, n, r, pU1, ldpu1, pVT1, ldpvt1, C1, ldc1, pU2, ldpu2, pVT2, ldpvt2, C2, ldc2, ierr)
    use maria_kinds_mod,   only: &
        WP => SP
    implicit none
        integer,        intent(in   )               :: m
        integer,        intent(in   )               :: n
        integer,        intent(in   )               :: r
        real(WP),       intent(in   ),  contiguous  :: pU1(:)
        integer,        intent(in   )               :: ldpu1
        real(WP),       intent(in   ),  contiguous  :: pVT1(:)
        integer,        intent(in   )               :: ldpvt1
        real(WP),       intent(in   ),  contiguous  :: C1(:)
        integer,        intent(in   )               :: ldc1
        real(WP),       intent(in   ),  contiguous  :: pU2(:)
        integer,        intent(in   )               :: ldpu2
        real(WP),       intent(in   ),  contiguous  :: pVT2(:)
        integer,        intent(in   )               :: ldpvt2
        real(WP),       intent(in   ),  contiguous  :: C2(:)
        integer,        intent(in   )               :: ldc2
        integer,        intent(  out)               :: ierr
        real(WP)                                    :: slrdotf_tangent
    end function slrdotf_tangent

    !------------------------------------------------------------------------------------------------------------------------

    !> Dot product of tangent vectors
    module function dlrdotf_tangent &
    (m, n, r, pU1, ldpu1, pVT1, ldpvt1, C1, ldc1, pU2, ldpu2, pVT2, ldpvt2, C2, ldc2, ierr)
    use maria_kinds_mod,   only: &
        WP => DP
    implicit none
        integer,        intent(in   )               :: m
        integer,        intent(in   )               :: n
        integer,        intent(in   )               :: r
        real(WP),       intent(in   ),  contiguous  :: pU1(:)
        integer,        intent(in   )               :: ldpu1
        real(WP),       intent(in   ),  contiguous  :: pVT1(:)
        integer,        intent(in   )               :: ldpvt1
        real(WP),       intent(in   ),  contiguous  :: C1(:)
        integer,        intent(in   )               :: ldc1
        real(WP),       intent(in   ),  contiguous  :: pU2(:)
        integer,        intent(in   )               :: ldpu2
        real(WP),       intent(in   ),  contiguous  :: pVT2(:)
        integer,        intent(in   )               :: ldpvt2
        real(WP),       intent(in   ),  contiguous  :: C2(:)
        integer,        intent(in   )               :: ldc2
        integer,        intent(  out)               :: ierr
        real(WP)                                    :: dlrdotf_tangent
    end function dlrdotf_tangent

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes a retraction onto the manfiold of fixed-rank matrices
    module subroutine slrretr_tangent &
    (m, n, r, S, U, ldu, VT, ldvt, alpha, pU, ldpu, pVT, ldpvt, C, ldc, work, lwork, iwork, liwork, ierr)
    use maria_kinds_mod,   only: &
        WP => SP
    implicit none
        integer,        intent(in   )               :: m
        integer,        intent(in   )               :: n
        integer,        intent(in   )               :: r
        real(WP),       intent(inout),  contiguous  :: S(:)
        real(WP),       intent(inout),  contiguous  :: U(:)
        integer,        intent(in   )               :: ldu
        real(WP),       intent(inout),  contiguous  :: VT(:)
        integer,        intent(in   )               :: ldvt
        real(WP),       intent(in   )               :: alpha
        real(WP),       intent(in   ),  contiguous  :: pU(:)
        integer,        intent(in   )               :: ldpu
        real(WP),       intent(in   ),  contiguous  :: pVT(:)
        integer,        intent(in   )               :: ldpvt
        real(WP),       intent(in   ),  contiguous  :: C(:)
        integer,        intent(in   )               :: ldc
        real(WP),       intent(  out),  contiguous  :: work(:)
        integer,        intent(in   )               :: lwork
        integer,        intent(  out),  contiguous  :: iwork(:)
        integer,        intent(in   )               :: liwork
        integer,        intent(  out)               :: ierr
    end subroutine slrretr_tangent

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes a retraction onto the manfiold of fixed-rank matrices
    module subroutine dlrretr_tangent &
    (m, n, r, S, U, ldu, VT, ldvt, alpha, pU, ldpu, pVT, ldpvt, C, ldc, work, lwork, iwork, liwork, ierr)
    use maria_kinds_mod,   only: &
        WP => DP
    implicit none
        integer,        intent(in   )               :: m
        integer,        intent(in   )               :: n
        integer,        intent(in   )               :: r
        real(WP),       intent(inout),  contiguous  :: S(:)
        real(WP),       intent(inout),  contiguous  :: U(:)
        integer,        intent(in   )               :: ldu
        real(WP),       intent(inout),  contiguous  :: VT(:)
        integer,        intent(in   )               :: ldvt
        real(WP),       intent(in   )               :: alpha
        real(WP),       intent(in   ),  contiguous  :: pU(:)
        integer,        intent(in   )               :: ldpu
        real(WP),       intent(in   ),  contiguous  :: pVT(:)
        integer,        intent(in   )               :: ldpvt
        real(WP),       intent(in   ),  contiguous  :: C(:)
        integer,        intent(in   )               :: ldc
        real(WP),       intent(  out),  contiguous  :: work(:)
        integer,        intent(in   )               :: lwork
        integer,        intent(  out),  contiguous  :: iwork(:)
        integer,        intent(in   )               :: liwork
        integer,        intent(  out)               :: ierr
    end subroutine dlrretr_tangent

    !------------------------------------------------------------------------------------------------------------------------
end interface
end module maria_lr_geom_mod
