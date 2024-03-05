!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_prng_builtin_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Defines a pseudorandom number generator based on Fortran intrinsic `random_number` subroutine
!> that extends the abstract [[maria_prng_mod(module):prng(type)]].
!----------------------------------------------------------------------------------------------------------------------------
module maria_prng_builtin_mod
use maria_prng_mod, only: &
    prng
implicit none (type, external)

! Types
public :: prng_builtin
private

!> A random number generator type based on Fortran intrinsic `random_number` subroutine.
!>
!> Provides access to subroutines for generating random numbers from multiple distributions:
!> 
!> - real-valued Gaussian
!> - real-valued and integer-valued uniform
type, extends(prng) :: prng_builtin
contains
    procedure :: init     => init_builtin
    procedure :: deinit   => deinit_builtin
    procedure :: snormal  => snormal_builtin
    procedure :: dnormal  => dnormal_builtin
    procedure :: iuniform => iuniform_builtin
    procedure :: suniform => suniform_builtin
    procedure :: duniform => duniform_builtin
end type prng_builtin

interface
    !------------------------------------------------------------------------------------------------------------------------    

    !> Initializes a pseudorandom number generator with a given random seed.
    module subroutine init_builtin &
    (this, seed, info)
    implicit none
        !> Pseudorandom number generator
        class(prng_builtin), intent(in)  :: this        
        !> Random seed
        integer,             intent(in)  :: seed
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info > 0`: system-required seed size is larger than assumed
        integer,             intent(out) :: info
    end subroutine init_builtin

    !------------------------------------------------------------------------------------------------------------------------    

    !> Denitializes a random number generator.
    module subroutine deinit_builtin &
    (this, info)
    implicit none
        !> Pseudorandom number generator
        class(prng_builtin), intent(in)  :: this        
        !> Exit code:
        !>
        !> - `info = 0`: successful exit
        integer,             intent(out) :: info
    end subroutine deinit_builtin

    !------------------------------------------------------------------------------------------------------------------------    

    !> Fills an array with single-precision random numbers
    !> from the Gaussian distribution \( \mathcal{N}( \text{mean}, \text{std}^2) \)
    !> with mean `mean` and standard deviation `std`.
    module subroutine snormal_builtin &
    (this, n, x, mean, std, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Pseudorandom number generator
        class(prng_builtin), intent(in)              :: this        
        !> Size of the array
        !>
        !> Possible values: \( n \geq 0 \)
        integer,             intent(in)              :: n
        !> Array of size \( n \)
        real(WP),            intent(out), contiguous :: x(:)
        !> Mean of the Gaussian distribution
        real(WP),            intent(in)              :: mean
        !> Standard deviation of the Gaussian distibution
        !>
        !> Possible values: \( \text{std} > 0 \)
        real(WP),            intent(in)              :: std
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,             intent(out)             :: info
    end subroutine snormal_builtin

    !------------------------------------------------------------------------------------------------------------------------    

    !> Fills an array with double-precision random numbers
    !> from the Gaussian distribution \( \mathcal{N}( \text{mean}, \text{std}^2) \)
    !> with mean `mean` and standard deviation `std`.
    module subroutine dnormal_builtin &
    (this, n, x, mean, std, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Pseudorandom number generator
        class(prng_builtin), intent(in)              :: this        
        !> Size of the array
        !>
        !> Possible values: \( n \geq 0 \)
        integer,             intent(in)              :: n
        !> Array of size \( n \)
        real(WP),            intent(out), contiguous :: x(:)
        !> Mean of the Gaussian distribution
        real(WP),            intent(in)              :: mean
        !> Standard deviation of the Gaussian distibution
        !>
        !> Possible values: \( \text{std} > 0 \)
        real(WP),            intent(in)              :: std
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,             intent(out)             :: info
    end subroutine dnormal_builtin

    !------------------------------------------------------------------------------------------------------------------------    

    !> Fills an array with integer random numbers
    !> from the uniform distribution on \( [a, b] \).
    module subroutine iuniform_builtin &
    (this, n, x, a, b, info)
    implicit none
        !> Pseudorandom number generator
        class(prng_builtin), intent(in)              :: this       
        !> Size of the array
        !>
        !> Possible values: \( n \geq 0 \)
        integer,             intent(in)              :: n
        !> Array of size \( n \)
        integer,             intent(out), contiguous :: x(:)
        !> Left boundary
        integer,             intent(in)              :: a
        !> Right boundary
        !>
        !> Possible values: \( b > a \)
        integer,             intent(in)              :: b
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,             intent(out)             :: info
    end subroutine iuniform_builtin

    !------------------------------------------------------------------------------------------------------------------------    

    !> Fills an array with single-precision random numbers
    !> from the uniform distribution on \( [a, b] \).
    module subroutine suniform_builtin &
    (this, n, x, a, b, info)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        !> Pseudorandom number generator
        class(prng_builtin), intent(in)              :: this       
        !> Size of the array
        !>
        !> Possible values: \( n \geq 0 \)
        integer,             intent(in)              :: n
        !> Array of size \( n \)
        real(WP),            intent(out), contiguous :: x(:)
        !> Left boundary
        real(WP),            intent(in)              :: a
        !> Right boundary
        !>
        !> Possible values: \( b > a \)
        real(WP),            intent(in)              :: b
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,             intent(out)             :: info
    end subroutine suniform_builtin

    !------------------------------------------------------------------------------------------------------------------------    

    !> Fills an array with double-precision random numbers
    !> from the uniform distribution on \( [a, b] \).
    module subroutine duniform_builtin &
    (this, n, x, a, b, info)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        !> Pseudorandom number generator
        class(prng_builtin), intent(in)              :: this       
        !> Size of the array
        !>
        !> Possible values: \( n \geq 0 \)
        integer,             intent(in)              :: n
        !> Array of size \( n \)
        real(WP),            intent(out), contiguous :: x(:)
        !> Left boundary
        real(WP),            intent(in)              :: a
        !> Right boundary
        !>
        !> Possible values: \( b > a \)
        real(WP),            intent(in)              :: b
        !> Exit code:
        !>
        !> - `info = 0`: successful exit,
        !> - `info < 0`: if `info = -i`, the \( i \)-th argument had an illegal value.
        integer,             intent(out)             :: info
    end subroutine duniform_builtin

    !------------------------------------------------------------------------------------------------------------------------    
end interface

end module maria_prng_builtin_mod
