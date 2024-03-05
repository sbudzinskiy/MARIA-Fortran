!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_assert_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Provides access to an assert subroutine used in unit tests.
!----------------------------------------------------------------------------------------------------------------------------
module maria_assert_mod
implicit none (type, external)

public :: assert_linefile
private

interface
    !------------------------------------------------------------------------------------------------------------------------

    !> Asserts that the input logical argument is `.true.`
    !>
    !> If `flag` is `.true.`, nothing happens.
    !>
    !> If `flag` is `.false.`, an error message is printed with the provided name of the file and line number.
    module subroutine assert_linefile &
    (flag, filename, line_number)
    implicit none
        !> Expression to be asserted
        logical,      intent(in) :: flag
        !> Name of the file
        character(*), intent(in) :: filename
        !> Number of the line in the file
        integer,      intent(in) :: line_number
    end subroutine assert_linefile

    !------------------------------------------------------------------------------------------------------------------------
end interface

end module maria_assert_mod
