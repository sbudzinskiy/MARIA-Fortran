!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_prng_mkl_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Defines a pseudorandom number generator based on Intel MKL VSL
!> that extends the abstract [[maria_prng_mod(module):prng(type)]].
!----------------------------------------------------------------------------------------------------------------------------
module maria_prng_mkl_mod
use maria_prng_mod, only: &
    prng
use mkl_vsl_type,   only: &
    vsl_stream_state
implicit none (type, external)

! Types
public :: prng_mkl
private

!> A pseudorandom number generator type based on Intel VSL MKL.
!>
!> Provides access to subroutines for generating random numbers from multiple distributions:
!> 
!> - real-valued Gaussian
!> - real-valued and integer-valued uniform
type, extends(prng) :: prng_mkl
private
    type(vsl_stream_state) :: stream
contains
    procedure :: init     => init_mkl
    procedure :: deinit   => deinit_mkl
    procedure :: snormal  => snormal_mkl
    procedure :: dnormal  => dnormal_mkl
    procedure :: iuniform => iuniform_mkl
    procedure :: suniform => suniform_mkl
    procedure :: duniform => duniform_mkl
end type prng_mkl

interface
    !------------------------------------------------------------------------------------------------------------------------

    !> Initializes a Mersenne Twister MT19937 pseudorandom number generator with a given random seed.
    module subroutine init_mkl &
    (this, seed, info)
    implicit none
        !> Pseudorandom number generator
        class(prng_mkl), intent(in)  :: this        
        !> Random seed
        integer,         intent(in)  :: seed
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info > 0`: a runtime error occured
        integer,         intent(out) :: info
    end subroutine init_mkl

    !------------------------------------------------------------------------------------------------------------------------

    !> Denitializes a random number generator.
    module subroutine deinit_mkl &
    (this, info)
    implicit none
        !> Pseudorandom number generator
        class(prng_mkl), intent(in)  :: this        
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info > 0`: a runtime error occured
        integer,         intent(out) :: info
    end subroutine deinit_mkl

    !------------------------------------------------------------------------------------------------------------------------

    !> Fills an array with single-precision random numbers
    !> from the Gaussian distribution \( \mathcal{N}( \text{mean}, \text{std}^2) \)
    !> with mean `mean` and standard deviation `std`.
    module subroutine snormal_mkl &
    (this, n, x, mean, std, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Pseudorandom number generator
        class(prng_mkl), intent(in)              :: this        
        !> Size of the array
        !>
        !> Possible values: \( n \geq 0 \)
        integer,         intent(in)              :: n
        !> Array of size \( n \)
        real(WP),        intent(out), contiguous :: x(:)
        !> Mean of the Gaussian distribution
        real(WP),        intent(in)              :: mean
        !> Standard deviation of the Gaussian distibution
        !>
        !> Possible values: \( \text{std} > 0 \)
        real(WP),        intent(in)              :: std
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,         intent(out)             :: info
    end subroutine snormal_mkl

    !------------------------------------------------------------------------------------------------------------------------

    !> Fills an array with double-precision random numbers
    !> from the Gaussian distribution \( \mathcal{N}( \text{mean}, \text{std}^2) \)
    !> with mean `mean` and standard deviation `std`.
    module subroutine dnormal_mkl &
    (this, n, x, mean, std, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Pseudorandom number generator
        class(prng_mkl), intent(in)              :: this        
        !> Size of the array
        !>
        !> Possible values: \( n \geq 0 \)
        integer,         intent(in)              :: n
        !> Array of size \( n \)
        real(WP),        intent(out), contiguous :: x(:)
        !> Mean of the Gaussian distribution
        real(WP),        intent(in)              :: mean
        !> Standard deviation of the Gaussian distibution
        !>
        !> Possible values: \( \text{std} > 0 \)
        real(WP),        intent(in)              :: std
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,         intent(out)             :: info
    end subroutine dnormal_mkl

    !------------------------------------------------------------------------------------------------------------------------

    !> Fills an array with integer random numbers
    !> from the uniform distribution on \( [a, b] \).
    module subroutine iuniform_mkl &
    (this, n, x, a, b, info)
    implicit none
        !> Pseudorandom number generator
        class(prng_mkl), intent(in)              :: this       
        !> Size of the array
        !>
        !> Possible values: \( n \geq 0 \)
        integer,         intent(in)              :: n
        !> Array of size \( n \)
        integer,         intent(out), contiguous :: x(:)
        !> Left boundary
        integer,         intent(in)              :: a
        !> Right boundary
        !>
        !> Possible values: \( b > a \)
        integer,         intent(in)              :: b
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,         intent(out)             :: info
    end subroutine iuniform_mkl

    !------------------------------------------------------------------------------------------------------------------------

    !> Fills an array with single-precision random numbers
    !> from the uniform distribution on \( [a, b] \).
    module subroutine suniform_mkl &
    (this, n, x, a, b, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Pseudorandom number generator
        class(prng_mkl), intent(in)              :: this       
        !> Size of the array
        !>
        !> Possible values: \( n \geq 0 \)
        integer,         intent(in)              :: n
        !> Array of size \( n \)
        real(WP),        intent(out), contiguous :: x(:)
        !> Left boundary
        real(WP),        intent(in)              :: a
        !> Right boundary
        !>
        !> Possible values: \( b > a \)
        real(WP),        intent(in)              :: b
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,         intent(out)             :: info
    end subroutine suniform_mkl

    !------------------------------------------------------------------------------------------------------------------------

    !> Fills an array with double-precision random numbers
    !> from the uniform distribution on \( [a, b] \).
    module subroutine duniform_mkl &
    (this, n, x, a, b, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Pseudorandom number generator
        class(prng_mkl), intent(in)              :: this       
        !> Size of the array
        !>
        !> Possible values: \( n \geq 0 \)
        integer,         intent(in)              :: n
        !> Array of size \( n \)
        real(WP),        intent(out), contiguous :: x(:)
        !> Left boundary
        real(WP),        intent(in)              :: a
        !> Right boundary
        !>
        !> Possible values: \( b > a \)
        real(WP),        intent(in)              :: b
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,         intent(out)             :: info
    end subroutine duniform_mkl

    !------------------------------------------------------------------------------------------------------------------------    
end interface

end module maria_prng_mkl_mod
