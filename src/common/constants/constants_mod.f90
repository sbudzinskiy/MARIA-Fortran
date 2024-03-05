!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_constants_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Defines various mathematical constants.
!----------------------------------------------------------------------------------------------------------------------------
module maria_constants_mod
use maria_kinds_mod, only: &
    SP,                    &
    DP
implicit none (type, external)

! Variables
public :: S_ONE,     &
          S_ZERO,    &    
          S_HALF,    &
          S_TWO,     &
          S_PI,      &
          S_TWOPI,   &
          S_MACHTOL, &
          D_ONE,     &
          D_ZERO,    &
          D_HALF,    &
          D_TWO,     &
          D_PI,      &
          D_TWOPI,   &
          D_MACHTOL
private      

!> Single-precision 1
real(SP), parameter :: S_ONE  = 1.0_SP
!> Single-precision 0
real(SP), parameter :: S_ZERO = 0.0_SP
!> Single-precision 0.5
real(SP), parameter :: S_HALF = 0.5_SP
!> Single-precision 2
real(SP), parameter :: S_TWO  = 2.0_SP
!> Single-precision \( \pi \)
real(SP), parameter :: S_PI = 4 * atan(S_ONE)
!> Single-precision \( 2\pi \)
real(SP), parameter :: S_TWOPI = S_TWO * S_PI
!> Single-precision unit roundoff
real(SP), parameter :: S_MACHTOL = epsilon(S_ONE)

!> Double-precision 1
real(DP), parameter :: D_ONE  = 1.0_DP
!> Double-precision 0
real(DP), parameter :: D_ZERO = 0.0_DP
!> Double-precision 0.5
real(DP), parameter :: D_HALF = 0.5_DP
!> Double-precision 2
real(DP), parameter :: D_TWO  = 2.0_DP
!> Double-precision \( \pi \)
real(DP), parameter :: D_PI = 4 * atan(D_ONE)
!> Double-precision \( 2\pi \)
real(DP), parameter :: D_TWOPI = D_TWO * D_PI
!> Double-precision unit roundoff
real(DP), parameter :: D_MACHTOL = epsilon(D_ONE)

end module maria_constants_mod
