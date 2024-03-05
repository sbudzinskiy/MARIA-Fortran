!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_prng_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Defines an abstract pseudorandom number generator.
!----------------------------------------------------------------------------------------------------------------------------
module maria_prng_mod
implicit none (type, external)

! Types
public :: prng
private

!> An abstract pseudorandom number generator type.
!>
!> Provides access to subroutines for generating random numbers from multiple distributions:
!>
!> - real-valued normal (Gaussian)
!> - real-valued and integer-valued uniform
type, abstract :: prng
contains
    procedure(init_a),     deferred :: init
    procedure(deinit_a),   deferred :: deinit
    procedure(snormal_a),  deferred :: snormal
    procedure(dnormal_a),  deferred :: dnormal
    procedure(iuniform_a), deferred :: iuniform
    procedure(suniform_a), deferred :: suniform
    procedure(duniform_a), deferred :: duniform
end type prng

abstract interface
    !------------------------------------------------------------------------------------------------------------------------

    !> Initializes a pseudorandom number generator with a given random seed.
    subroutine init_a &
    (this, seed, info)
    import prng
    implicit none
        !> Pseudorandom number generator
        class(prng), intent(in)  :: this        
        !> Random seed
        integer,     intent(in)  :: seed
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info > 0`: a runtime error occured
        integer,     intent(out) :: info
    end subroutine init_a

    !------------------------------------------------------------------------------------------------------------------------

    !> Denitializes a pseudorandom number generator.
    subroutine deinit_a &
    (this, info)
    import prng
    implicit none
        !> Pseudorandom number generator
        class(prng), intent(in)  :: this        
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info > 0`: a runtime error occured
        integer,     intent(out) :: info
    end subroutine deinit_a

    !------------------------------------------------------------------------------------------------------------------------

    !> Fills an array with single-precision random numbers
    !> from the Gaussian distribution \( \mathcal{N}( \text{mean}, \text{std}^2) \)
    !> with mean `mean` and standard deviation `std`.
    subroutine snormal_a &
    (this, n, x, mean, std, info)
    use maria_kinds_mod, only: &
        WP => SP
    import prng
    implicit none
        !> Pseudorandom number generator
        class(prng), intent(in)              :: this        
        !> Size of the array
        !>
        !> Possible values: \( n \geq 0 \)
        integer,     intent(in)              :: n
        !> Array of size \( n \)
        real(WP),    intent(out), contiguous :: x(:)
        !> Mean of the Gaussian distribution
        real(WP),    intent(in)              :: mean
        !> Standard deviation of the Gaussian distibution
        !>
        !> Possible values: \( \text{std} > 0 \)
        real(WP),    intent(in)              :: std
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,     intent(out)             :: info
    end subroutine snormal_a

    !------------------------------------------------------------------------------------------------------------------------

    !> Fills an array with double-precision random numbers
    !> from the Gaussian distribution \( \mathcal{N}( \text{mean}, \text{std}^2) \)
    !> with mean `mean` and standard deviation `std`.
    subroutine dnormal_a &
    (this, n, x, mean, std, info)
    use maria_kinds_mod, only: &
        WP => DP
    import prng
    implicit none
        !> Pseudorandom number generator
        class(prng), intent(in)              :: this        
        !> Size of the array
        !>
        !> Possible values: \( n \geq 0 \)
        integer,     intent(in)              :: n
        !> Array of size \( n \)
        real(WP),    intent(out), contiguous :: x(:)
        !> Mean of the Gaussian distribution
        real(WP),    intent(in)              :: mean
        !> Standard deviation of the Gaussian distibution
        !>
        !> Possible values: \( \text{std} > 0 \)
        real(WP),    intent(in)              :: std
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,     intent(out)             :: info
    end subroutine dnormal_a

    !------------------------------------------------------------------------------------------------------------------------

    !> Fills an array with integer random numbers
    !> from the uniform distribution on \( [a, b] \).
    subroutine iuniform_a &
    (this, n, x, a, b, info)
    import prng
    implicit none
        !> Pseudorandom number generator
        class(prng), intent(in)              :: this       
        !> Size of the array
        !>
        !> Possible values: \( n \geq 0 \)
        integer,     intent(in)              :: n
        !> Array of size \( n \)
        integer,     intent(out), contiguous :: x(:)
        !> Left boundary
        integer,     intent(in)              :: a
        !> Right boundary
        !>
        !> Possible values: \( b > a \)
        integer,     intent(in)              :: b
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,     intent(out)             :: info
    end subroutine iuniform_a

    !------------------------------------------------------------------------------------------------------------------------

    !> Fills an array with single-precision random numbers
    !> from the uniform distribution on \( [a, b] \).
    subroutine suniform_a &
    (this, n, x, a, b, info)
    use maria_kinds_mod, only: &
        WP => SP
    import prng
    implicit none
        !> Pseudorandom number generator
        class(prng), intent(in)              :: this       
        !> Size of the array
        !>
        !> Possible values: \( n \geq 0 \)
        integer,     intent(in)              :: n
        !> Array of size \( n \)
        real(WP),    intent(out), contiguous :: x(:)
        !> Left boundary
        real(WP),    intent(in)              :: a
        !> Right boundary
        !>
        !> Possible values: \( b > a \)
        real(WP),    intent(in)              :: b
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,     intent(out)             :: info
    end subroutine suniform_a

    !------------------------------------------------------------------------------------------------------------------------

    !> Fills an array with double-precision random numbers
    !> from the uniform distribution on \( [a, b] \).
    subroutine duniform_a &
    (this, n, x, a, b, info)
    use maria_kinds_mod, only: &
        WP => DP
    import prng
    implicit none
        !> Pseudorandom number generator
        class(prng), intent(in)              :: this       
        !> Size of the array
        !>
        !> Possible values: \( n \geq 0 \)
        integer,     intent(in)              :: n
        !> Array of size \( n \)
        real(WP),    intent(out), contiguous :: x(:)
        !> Left boundary
        real(WP),    intent(in)              :: a
        !> Right boundary
        !>
        !> Possible values: \( b \geq a \)
        real(WP),    intent(in)              :: b
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,     intent(out)             :: info
    end subroutine duniform_a

    !------------------------------------------------------------------------------------------------------------------------
end interface

end module maria_prng_mod
