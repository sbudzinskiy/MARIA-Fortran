!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the public interface of the [[maria_access_tensor_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Interfaces for accessing matrices via their elements, columns, and rows.
!----------------------------------------------------------------------------------------------------------------------------
module maria_access_tensor_mod
implicit none (type, external)

! Interfaces
public :: stenval, &
          dtenval, &
          stenfib, &
          dtenfib

! Procedures
public ::   mi2i,               &
            i2mi,               &
            stenval2fib,        &
            dtenval2fib,        &
            stenval2full,       &
            dtenval2full,       &
            stenval_hilbert,    &
            dtenval_hilbert,    &
            stenval_sum,        &
            dtenval_sum

abstract interface
    !------------------------------------------------------------------------------------------------------------------------

    function stenval &
    (d, n, mi, ierr)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,    intent(in   )               :: d
        integer,    intent(in   ),  contiguous  :: n(:)
        integer,    intent(in   ),  contiguous  :: mi(:)
        integer,    intent(  out)               :: ierr
        real(WP)                                :: stenval
    end function stenval

    !------------------------------------------------------------------------------------------------------------------------

    function dtenval &
    (d, n, mi, ierr)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,    intent(in   )               :: d
        integer,    intent(in   ),  contiguous  :: n(:)
        integer,    intent(in   ),  contiguous  :: mi(:)
        integer,    intent(  out)               :: ierr
        real(WP)                                :: dtenval
    end function dtenval

    !------------------------------------------------------------------------------------------------------------------------

    subroutine stenfib &
    (d, n, k, li, ri, x, incx, ierr)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,    intent(in   )               :: d
        integer,    intent(in   ),  contiguous  :: n(:)
        integer,    intent(in   )               :: k
        integer,    intent(in   ),  contiguous  :: li(:)
        integer,    intent(in   ),  contiguous  :: ri(:)
        real(WP),   intent(  out),  contiguous  :: x(:)
        integer,    intent(in   )               :: incx
        integer,    intent(  out)               :: ierr
    end subroutine stenfib

    !------------------------------------------------------------------------------------------------------------------------

    subroutine dtenfib &
    (d, n, k, li, ri, x, incx, ierr)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,    intent(in   )               :: d
        integer,    intent(in   ),  contiguous  :: n(:)
        integer,    intent(in   )               :: k
        integer,    intent(in   ),  contiguous  :: li(:)
        integer,    intent(in   ),  contiguous  :: ri(:)
        real(WP),   intent(  out),  contiguous  :: x(:)
        integer,    intent(in   )               :: incx
        integer,    intent(  out)               :: ierr
    end subroutine dtenfib

    !------------------------------------------------------------------------------------------------------------------------
end interface

interface
    !------------------------------------------------------------------------------------------------------------------------

    module subroutine mi2i  &
    (d, n, mi, i, ierr)
    implicit none
        integer,    intent(in   )               :: d
        integer,    intent(in   ),  contiguous  :: n(:)
        integer,    intent(in   ),  contiguous  :: mi(:)
        integer,    intent(  out)               :: i
        integer,    intent(  out)               :: ierr
    end subroutine mi2i

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine i2mi  &
    (d, n, i, mi, ierr)
    implicit none
        integer,    intent(in   )               :: d
        integer,    intent(in   ),  contiguous  :: n(:)
        integer,    intent(in   )               :: i
        integer,    intent(  out),  contiguous  :: mi(:)
        integer,    intent(  out)               :: ierr
    end subroutine i2mi

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine stenval2fib &
    (tenval, d, n, k, li, ri, x, incx, iwork, liwork, ierr)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        procedure(stenval), intent(in   ),  pointer     :: tenval
        integer,            intent(in   )               :: d
        integer,            intent(in   ),  contiguous  :: n(:)
        integer,            intent(in   )               :: k
        integer,            intent(in   ),  contiguous  :: li(:)
        integer,            intent(in   ),  contiguous  :: ri(:)
        real(WP),           intent(  out),  contiguous  :: x(:)
        integer,            intent(in   )               :: incx
        integer,            intent(  out),  contiguous  :: iwork(:)
        integer,            intent(in   )               :: liwork
        integer,            intent(  out)               :: ierr
    end subroutine stenval2fib

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dtenval2fib &
    (tenval, d, n, k, li, ri, x, incx, iwork, liwork, ierr)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        procedure(dtenval), intent(in   ),  pointer     :: tenval
        integer,            intent(in   )               :: d
        integer,            intent(in   ),  contiguous  :: n(:)
        integer,            intent(in   )               :: k
        integer,            intent(in   ),  contiguous  :: li(:)
        integer,            intent(in   ),  contiguous  :: ri(:)
        real(WP),           intent(  out),  contiguous  :: x(:)
        integer,            intent(in   )               :: incx
        integer,            intent(  out),  contiguous  :: iwork(:)
        integer,            intent(in   )               :: liwork
        integer,            intent(  out)               :: ierr
    end subroutine dtenval2fib

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine stenval2full &
    (tenval, d, n, A, iwork, liwork, ierr)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        procedure(stenval), intent(in   ),  pointer     :: tenval
        integer,            intent(in   )               :: d
        integer,            intent(in   ),  contiguous  :: n(:)
        real(WP),           intent(  out),  contiguous  :: A(:)
        integer,            intent(  out),  contiguous  :: iwork(:)
        integer,            intent(in   )               :: liwork
        integer,            intent(  out)               :: ierr
    end subroutine stenval2full

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dtenval2full &
    (tenval, d, n, A, iwork, liwork, ierr)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        procedure(dtenval), intent(in   ),  pointer     :: tenval
        integer,            intent(in   )               :: d
        integer,            intent(in   ),  contiguous  :: n(:)
        real(WP),           intent(  out),  contiguous  :: A(:)
        integer,            intent(  out),  contiguous  :: iwork(:)
        integer,            intent(in   )               :: liwork
        integer,            intent(  out)               :: ierr
    end subroutine dtenval2full

    !------------------------------------------------------------------------------------------------------------------------

    module function stenval_hilbert &
    (d, n, mi, ierr)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,    intent(in   )               :: d
        integer,    intent(in   ),  contiguous  :: n(:)
        integer,    intent(in   ),  contiguous  :: mi(:)
        integer,    intent(  out)               :: ierr
        real(WP)                                :: stenval_hilbert
    end function stenval_hilbert

    !------------------------------------------------------------------------------------------------------------------------

    module function dtenval_hilbert &
    (d, n, mi, ierr)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,    intent(in   )               :: d
        integer,    intent(in   ),  contiguous  :: n(:)
        integer,    intent(in   ),  contiguous  :: mi(:)
        integer,    intent(  out)               :: ierr
        real(WP)                                :: dtenval_hilbert
    end function dtenval_hilbert

    !------------------------------------------------------------------------------------------------------------------------

    module function stenval_sum &
    (d, n, mi, ierr)
    use maria_kinds_mod, only: &
        WP => SP
    implicit none
        integer,    intent(in   )               :: d
        integer,    intent(in   ),  contiguous  :: n(:)
        integer,    intent(in   ),  contiguous  :: mi(:)
        integer,    intent(  out)               :: ierr
        real(WP)                                :: stenval_sum
    end function stenval_sum

    !------------------------------------------------------------------------------------------------------------------------

    module function dtenval_sum &
    (d, n, mi, ierr)
    use maria_kinds_mod, only: &
        WP => DP
    implicit none
        integer,    intent(in   )               :: d
        integer,    intent(in   ),  contiguous  :: n(:)
        integer,    intent(in   ),  contiguous  :: mi(:)
        integer,    intent(  out)               :: ierr
        real(WP)                                :: dtenval_sum
    end function dtenval_sum

    !------------------------------------------------------------------------------------------------------------------------
end interface

end module maria_access_tensor_mod
