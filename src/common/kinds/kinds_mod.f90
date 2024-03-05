!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_kinds_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Defines single-precision and double-precision kinds.
!----------------------------------------------------------------------------------------------------------------------------
module maria_kinds_mod
implicit none (type, external)

! Variables
public :: SP, &
          DP
private

!> Single-precision kind
integer, parameter :: SP = selected_real_kind(p=6)
!> Double-precision kind
integer, parameter :: DP = selected_real_kind(p=12)

end module maria_kinds_mod
