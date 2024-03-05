!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_lr_tsvd_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Routines for computing the truncated SVD.
!----------------------------------------------------------------------------------------------------------------------------
module maria_lr_tsvd_mod
implicit none (type, external)

! Procedures
public :: schop,    &
          dchop,    &
          sgersvd1, &
          dgersvd1, &
          sgersvd2, &
          dgersvd2

interface
    !------------------------------------------------------------------------------------------------------------------------

    !> Given the singular values, computes the truncation rank based on several criteria and the 
    !> corresponding truncation error.
    module function schop &
    (n, S, info, maxr, rtolf, atolf, rtol2, atol2, rerrf, aerrf, rerr2, aerr2)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Number of singular values
        integer,  intent(in)              :: n
        !> The singular values in descending order
        real(WP), intent(in),  contiguous :: S(:)
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)             :: info
        !> Largest possible truncation rank
        integer,  intent(in),  optional   :: maxr
        !> Relative error tolerance in Frobenius norm
        real(WP), intent(in),  optional   :: rtolf
        !> Absolute error tolerance in Frobenius norm
        real(WP), intent(in),  optional   :: atolf
        !> Relative error tolerance in spectral norm
        real(WP), intent(in),  optional   :: rtol2
        !> Absolute error tolerance in spectral norm
        real(WP), intent(in),  optional   :: atol2
        !> Truncation relative error in Frobenius norm
        real(WP), intent(out), optional   :: rerrf
        !> Truncation absolute error in Frobenius norm
        real(WP), intent(out), optional   :: aerrf
        !> Truncation relative error in spectral norm
        real(WP), intent(out), optional   :: rerr2
        !> Truncation absolute error in spectral norm
        real(WP), intent(out), optional   :: aerr2
        integer                           :: schop
    end function schop

    !------------------------------------------------------------------------------------------------------------------------

    !> Given the singular values, computes the truncation rank based on several criteria and the 
    !> corresponding truncation error.
    module function dchop &
    (n, S, info, maxr, rtolf, atolf, rtol2, atol2, rerrf, aerrf, rerr2, aerr2)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Number of singular values
        integer,  intent(in)              :: n
        !> The singular values in descending order
        real(WP), intent(in),  contiguous :: S(:)
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,  intent(out)             :: info
        !> Largest possible truncation rank
        integer,  intent(in),  optional   :: maxr
        !> Relative error tolerance in Frobenius norm
        real(WP), intent(in),  optional   :: rtolf
        !> Absolute error tolerance in Frobenius norm
        real(WP), intent(in),  optional   :: atolf
        !> Relative error tolerance in spectral norm
        real(WP), intent(in),  optional   :: rtol2
        !> Absolute error tolerance in spectral norm
        real(WP), intent(in),  optional   :: atol2
        !> Truncation relative error in Frobenius norm
        real(WP), intent(out), optional   :: rerrf
        !> Truncation absolute error in Frobenius norm
        real(WP), intent(out), optional   :: aerrf
        !> Truncation relative error in spectral norm
        real(WP), intent(out), optional   :: rerr2
        !> Truncation absolute error in spectral norm
        real(WP), intent(out), optional   :: aerr2
        integer                           :: dchop
    end function dchop

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes an approximate singular value decomposition of a matrix using randomized Algorithm 5.1 from 
    !>
    !>   N. Halko, P. G. Martinsson, and J. A. Tropp, “Finding Structure with Randomness: Probabilistic Algorithms for
    !>   Constructing Approximate Matrix Decompositions,” SIAM Review, vol. 53, pp. 217–288, Jan. 2011.
    module subroutine sgersvd1 &
    (m, n, B2AB, B2BA, kc, csk, ldcsk, niter_power, S, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,   only: &
        WP => SP
    use maria_la_core_mod, only: &
        MM => smatmul
    implicit none
        integer,       intent(in)                :: m
        integer,       intent(in)                :: n
        procedure(MM), intent(in),    pointer    :: B2AB
        procedure(MM), intent(in),    pointer    :: B2BA
        integer,       intent(in)                :: kc
        real(WP),      intent(inout), contiguous :: csk(:)
        integer,       intent(in)                :: ldcsk
        integer,       intent(in)                :: niter_power
        real(WP),      intent(out),   contiguous :: S(:)
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
    end subroutine sgersvd1

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes an approximate singular value decomposition of a matrix using randomized Algorithm 5.1 from 
    !>
    !>   N. Halko, P. G. Martinsson, and J. A. Tropp, “Finding Structure with Randomness: Probabilistic Algorithms for
    !>   Constructing Approximate Matrix Decompositions,” SIAM Review, vol. 53, pp. 217–288, Jan. 2011.
    module subroutine dgersvd1 &
    (m, n, B2AB, B2BA, kc, csk, ldcsk, niter_power, S, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,   only: &
        WP => DP
    use maria_la_core_mod, only: &
        MM => dmatmul
    implicit none
        integer,       intent(in)                :: m
        integer,       intent(in)                :: n
        procedure(MM), intent(in),    pointer    :: B2AB
        procedure(MM), intent(in),    pointer    :: B2BA
        integer,       intent(in)                :: kc
        real(WP),      intent(inout), contiguous :: csk(:)
        integer,       intent(in)                :: ldcsk
        integer,       intent(in)                :: niter_power
        real(WP),      intent(out),   contiguous :: S(:)
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
    end subroutine dgersvd1

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes an approximate singular value decomposition of a matrix using randomized Algorithm 4 from 
    !>
    !>  J. A. Tropp, A. Yurtsever, M. Udell, and V. Cevher, “Practical sketching algorithms for low-rank
    !>  matrix approximation,” SIAM Journal on Matrix Analysis and Applications, vol. 38, pp. 1454–1485, Jan. 2017.
    module subroutine sgersvd2 &
    (m, n, kc, csk, ldcsk, test_row, kr, rsk, ldrsk, S, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,   only: &
        WP => SP
    use maria_la_core_mod, only: &
        MM => smatmul
    implicit none
        integer,       intent(in)                :: m
        integer,       intent(in)                :: n
        integer,       intent(in)                :: kc
        real(WP),      intent(inout), contiguous :: csk(:)
        integer,       intent(in)                :: ldcsk
        procedure(MM), intent(in),    pointer    :: test_row
        integer,       intent(in)                :: kr
        real(WP),      intent(inout), contiguous :: rsk(:)
        integer,       intent(in)                :: ldrsk
        real(WP),      intent(out),   contiguous :: S(:)
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
    end subroutine sgersvd2

    !------------------------------------------------------------------------------------------------------------------------

    !> Computes an approximate singular value decomposition of a matrix using randomized Algorithm 4 from 
    !>
    !>  J. A. Tropp, A. Yurtsever, M. Udell, and V. Cevher, “Practical sketching algorithms for low-rank
    !>  matrix approximation,” SIAM Journal on Matrix Analysis and Applications, vol. 38, pp. 1454–1485, Jan. 2017.
    module subroutine dgersvd2 &
    (m, n, kc, csk, ldcsk, test_row, kr, rsk, ldrsk, S, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
    use maria_kinds_mod,   only: &
        WP => DP
    use maria_la_core_mod, only: &
        MM => dmatmul
    implicit none
        integer,       intent(in)                :: m
        integer,       intent(in)                :: n
        integer,       intent(in)                :: kc
        real(WP),      intent(inout), contiguous :: csk(:)
        integer,       intent(in)                :: ldcsk
        procedure(MM), intent(in),    pointer    :: test_row
        integer,       intent(in)                :: kr
        real(WP),      intent(inout), contiguous :: rsk(:)
        integer,       intent(in)                :: ldrsk
        real(WP),      intent(out),   contiguous :: S(:)
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
    end subroutine dgersvd2

    !------------------------------------------------------------------------------------------------------------------------
end interface
end module maria_lr_tsvd_mod
