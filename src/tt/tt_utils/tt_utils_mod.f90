!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_tt_utils_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Provides access to TT
!----------------------------------------------------------------------------------------------------------------------------
module maria_tt_utils_mod
implicit none (type, external)

! Procedures
public :: ttrank, stt2full, dtt2full

interface
    !------------------------------------------------------------------------------------------------------------------------

    !> Computes the volume of a matrix
    module subroutine ttrank &
    (d, r, k, rl, rr, ierr, r0, rd)
    implicit none
        integer,    intent(in   )               :: d
        integer,    intent(in   ),  contiguous  :: r(:)
        integer,    intent(in   )               :: k
        integer,    intent(  out)               :: rl
        integer,    intent(  out)               :: rr
        integer,    intent(  out)               :: ierr
        integer,    intent(in   ),  optional    :: r0
        integer,    intent(in   ),  optional    :: rd
    end subroutine ttrank

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes full tensor from TT
    module subroutine stt2full &
    (d, n, r, cores, A, work, lwork, ierr)
    use maria_kinds_mod,    only:   &
        WP  =>  SP
    use maria_arr_mod,      only:   & 
        AR  =>  sarr
    implicit none
        integer,    intent(in   )               :: d
        integer,    intent(in   ),  contiguous  :: n(:)
        integer,    intent(in   ),  contiguous  :: r(:)
        type(AR),   intent(in   )               :: cores(:)
        real(WP),   intent(  out),  contiguous  :: A(:)
        real(WP),   intent(  out),  contiguous  :: work(:)
        integer,    intent(in   )               :: lwork
        integer,    intent(  out)               :: ierr
    end subroutine stt2full

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes full tensor from TT
    module subroutine dtt2full &
    (d, n, r, cores, A, work, lwork, ierr)
    use maria_kinds_mod,    only:   &
        WP  =>  DP
    use maria_arr_mod,      only:   & 
        AR  =>  darr
    implicit none
        integer,    intent(in   )               :: d
        integer,    intent(in   ),  contiguous  :: n(:)
        integer,    intent(in   ),  contiguous  :: r(:)
        type(AR),   intent(in   )               :: cores(:)
        real(WP),   intent(  out),  contiguous  :: A(:)
        real(WP),   intent(  out),  contiguous  :: work(:)
        integer,    intent(in   )               :: lwork
        integer,    intent(  out)               :: ierr
    end subroutine dtt2full

    !------------------------------------------------------------------------------------------------------------------------
end interface
end module maria_tt_utils_mod
