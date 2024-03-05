!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_tt_tsvd_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> TTSVD
!----------------------------------------------------------------------------------------------------------------------------
module maria_tt_tsvd_mod
implicit none (type, external)

! Procedures
public :: sttsvd, dttsvd

interface
    !------------------------------------------------------------------------------------------------------------------------

    !> Computes TT from full tensor
    module subroutine sttsvd &
    (ort, d, n, A, r, cores, work, lwork, iwork, liwork, ierr, maxr, rtolf, rerrf)
    use maria_kinds_mod,    only:   &
        WP  =>  SP
    use maria_arr_mod,      only:   & 
        AR  =>  sarr
    implicit none
        character(1),   intent(in   )                           :: ort
        integer,        intent(in   )                           :: d
        integer,        intent(in   ),  contiguous              :: n(:)
        real(WP),       intent(inout),  contiguous              :: A(:)
        integer,        intent(  out),  contiguous              :: r(:)
        type(AR),       intent(  out)                           :: cores(:)
        real(WP),       intent(  out),  contiguous              :: work(:)
        integer,        intent(in   )                           :: lwork
        integer,        intent(  out),  contiguous              :: iwork(:)
        integer,        intent(in   )                           :: liwork
        integer,        intent(  out)                           :: ierr
        integer,        intent(in   ),  contiguous, optional    :: maxr(:)
        real(WP),       intent(in   ),              optional    :: rtolf
        real(WP),       intent(  out),              optional    :: rerrf
    end subroutine sttsvd

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes TT from full tensor
    module subroutine dttsvd &
    (ort, d, n, A, r, cores, work, lwork, iwork, liwork, ierr, maxr, rtolf, rerrf)
    use maria_kinds_mod,    only:   &
        WP  =>  DP
    use maria_arr_mod,      only:   & 
        AR  =>  darr
    implicit none
        character(1),   intent(in   )                           :: ort
        integer,        intent(in   )                           :: d
        integer,        intent(in   ),  contiguous              :: n(:)
        real(WP),       intent(inout),  contiguous              :: A(:)
        integer,        intent(  out),  contiguous              :: r(:)
        type(AR),       intent(  out)                           :: cores(:)
        real(WP),       intent(  out),  contiguous              :: work(:)
        integer,        intent(in   )                           :: lwork
        integer,        intent(  out),  contiguous              :: iwork(:)
        integer,        intent(in   )                           :: liwork
        integer,        intent(  out)                           :: ierr
        integer,        intent(in   ),  contiguous, optional    :: maxr(:)
        real(WP),       intent(in   ),              optional    :: rtolf
        real(WP),       intent(  out),              optional    :: rerrf
    end subroutine dttsvd

    !------------------------------------------------------------------------------------------------------------------------
end interface
end module maria_tt_tsvd_mod
