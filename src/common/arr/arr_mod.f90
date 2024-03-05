!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_arr_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Defines types for allocatable arrays of floating-point numbers.
!----------------------------------------------------------------------------------------------------------------------------
module maria_arr_mod
use maria_kinds_mod, only: &
    SP, &
    DP
implicit none (type, external)

! Types
public :: sarr, &
          darr, &
          iarr
private

!> A type for allocatable arrays of single-precision numbers.
type sarr
    real(SP), allocatable :: arr(:)
end type sarr

!> A type for allocatable arrays of double-precision numbers.
type darr
    real(DP), allocatable :: arr(:)
end type darr

type iarr
    integer, allocatable :: arr(:)
end type iarr

end module maria_arr_mod
