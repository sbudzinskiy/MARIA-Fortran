!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the implementation of the [[maria_prng_builtin_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Implements the [[maria_prng_builtin_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
submodule (maria_prng_builtin_mod) maria_prng_builtin_sub
implicit none (type, external)

integer, parameter :: MAX_SEED_SIZE = 16

contains
    !------------------------------------------------------------------------------------------------------------------------

    module subroutine init_builtin &
    (this, seed, info)
    use maria_reports_mod,   only:  &
        MAX_SEED_EXCEEDED_ERR_CODE, &
        report_runtime_err
    !-- Input/output arguments -------------------------------------------------
        class(prng_builtin), intent(in)  :: this        
        integer,             intent(in)  :: seed
        integer,             intent(out) :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "INIT_BUILTIN"
        integer :: seed_size, seed_vec(MAX_SEED_SIZE)

    !-- Executable section -----------------------------------------------------
        info = 0
        call random_seed(size=seed_size)        
        if (seed_size > MAX_SEED_SIZE) then
            info = 1
            call report_runtime_err(SRNAME, MAX_SEED_EXCEEDED_ERR_CODE)
            return
        end if

        seed_vec = seed
        call random_seed(put=seed_vec(1:seed_size))

    !-- Dismisses the unused dummy argument warning ----------------------------
        associate(this => this)
        end associate
    end subroutine init_builtin

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine deinit_builtin &
    (this, info)
    !-- Input/output arguments -------------------------------------------------
        class(prng_builtin), intent(in)  :: this        
        integer,             intent(out) :: info

    !-- Executable section -----------------------------------------------------
        info = 0

    !-- Dismisses the unused dummy argument warning ----------------------------
        associate(this => this)
        end associate
    end subroutine deinit_builtin

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine snormal_builtin &
    (this, n, x, mean, std, info)
    use maria_kinds_mod,     only: &
        WP => SP
    use maria_constants_mod, only: &
        ZERO => S_ZERO,            &
        TWOPI => S_TWOPI
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS,               &
        BAD_IF_SAME
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        class(prng_builtin), intent(in)              :: this        
        integer,             intent(in)              :: n
        real(WP),            intent(out), contiguous :: x(:)
        real(WP),            intent(in)              :: mean
        real(WP),            intent(in)              :: std
        integer,             intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "SNORMAL_BUILTIN"
        integer                 :: mod2, npairs, i
        real(WP)                :: u, v

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, std, ZERO)) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_SAME, std, ZERO)) then
            info = -5
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (n == 0) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        call random_number(x(1:n))

        mod2 = mod(n, 2)
        npairs = (n - mod2) / 2
        do i = 1, npairs
            u = x(2*i - 1)
            v = x(2*i)
            x(2*i - 1) = sqrt(-2 * log(u)) * cos(TWOPI * v)
            x(2*i)     = sqrt(-2 * log(u)) * sin(TWOPI * v)
        end do

        if (mod(n, 2) == 1) then
            call random_number(v)
            x(n) = sqrt(-2 * log(x(n))) * cos(TWOPI * v)
        end if

        x(1:n) = std * x(1:n) + mean

    !-- Dismisses the unused dummy argument warning ----------------------------
        associate(this => this)
        end associate
    end subroutine snormal_builtin

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dnormal_builtin &
    (this, n, x, mean, std, info)
    use maria_kinds_mod,     only: &
        WP => DP
    use maria_constants_mod, only: &
        ZERO => D_ZERO,            &
        TWOPI => D_TWOPI
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS,               &
        BAD_IF_SAME
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        class(prng_builtin), intent(in)              :: this        
        integer,             intent(in)              :: n
        real(WP),            intent(out), contiguous :: x(:)
        real(WP),            intent(in)              :: mean
        real(WP),            intent(in)              :: std
        integer,             intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "DNORMAL_BUILTIN"
        integer                 :: mod2, npairs, i
        real(WP)                :: u, v

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, std, ZERO)) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_SAME, std, ZERO)) then
            info = -5
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (n == 0) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        call random_number(x(1:n))

        mod2 = mod(n, 2)
        npairs = (n - mod2) / 2
        do i = 1, npairs
            u = x(2*i - 1)
            v = x(2*i)
            x(2*i - 1) = sqrt(-2 * log(u)) * cos(TWOPI * v)
            x(2*i)     = sqrt(-2 * log(u)) * sin(TWOPI * v)
        end do

        if (mod(n, 2) == 1) then
            call random_number(v)
            x(n) = sqrt(-2 * log(x(n))) * cos(TWOPI * v)
        end if

        x(1:n) = std * x(1:n) + mean

    !-- Dismisses the unused dummy argument warning ----------------------------
        associate(this => this)
        end associate
    end subroutine dnormal_builtin

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine iuniform_builtin &
    (this, n, x, a, b, info)
    use maria_argcheck_mod, only: &
        arg_is_bad,               &
        BAD_IF_LESS,              &
        BAD_IF_SAME
    use maria_reports_mod,  only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        class(prng_builtin), intent(in)              :: this        
        integer,             intent(in)              :: n
        integer,             intent(out), contiguous :: x(:)
        integer,             intent(in)              :: a
        integer,             intent(in)              :: b
        integer,             intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "IUNIFORM_BUILTIN"
        integer                 :: i
        real                    :: y

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, b, a)) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_SAME, b, a)) then
            info = -5
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (n == 0) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        do i = 1, n
            call random_number(y)
            x(i) = nint(a + (b - a) * y)
        end do

    !-- Dismisses the unused dummy argument warning ----------------------------
        associate(this => this)
        end associate
    end subroutine iuniform_builtin

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine suniform_builtin &
    (this, n, x, a, b, info)
    use maria_kinds_mod,    only: &
        WP => SP
    use maria_argcheck_mod, only: &
        arg_is_bad,               &
        BAD_IF_LESS,              &
        BAD_IF_SAME
    use maria_reports_mod,  only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        class(prng_builtin), intent(in)              :: this        
        integer,             intent(in)              :: n
        real(WP),            intent(out), contiguous :: x(:)
        real(WP),            intent(in)              :: a
        real(WP),            intent(in)              :: b
        integer,             intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "SUNIFORM_BUILTIN"

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, b, a)) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_SAME, b, a)) then
            info = -5
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (n == 0) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        call random_number(x(1:n))
        x(1:n) = a + (b - a) * x(1:n)

    !-- Dismisses the unused dummy argument warning ----------------------------
        associate(this => this)
        end associate
    end subroutine suniform_builtin

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine duniform_builtin &
    (this, n, x, a, b, info)
    use maria_kinds_mod,    only: &
        WP => DP
    use maria_argcheck_mod, only: &
        arg_is_bad,               &
        BAD_IF_LESS,              &
        BAD_IF_SAME
    use maria_reports_mod,  only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        class(prng_builtin), intent(in)              :: this        
        integer,             intent(in)              :: n
        real(WP),            intent(out), contiguous :: x(:)
        real(WP),            intent(in)              :: a
        real(WP),            intent(in)              :: b
        integer,             intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "DUNIFORM_BUILTIN"

    !-- Sanity check -----------------------------------------------------------
    sanity: block
        info = 0
        if (arg_is_bad(BAD_IF_LESS, n, 0)) then
            info = -2
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_LESS, b, a)) then
            info = -5
            exit sanity
        end if
        if (arg_is_bad(BAD_IF_SAME, b, a)) then
            info = -5
            exit sanity
        end if
    end block sanity

    !-- Quick return if possible -----------------------------------------------
    quickr: block
        if (info /= 0) exit quickr

        if (n == 0) return
    end block quickr

    !-- Report bad input -------------------------------------------------------
        if (info /= 0) then
            call report_bad_arg(SRNAME, -info)
            return
        end if

    !-- Executable section -----------------------------------------------------
        call random_number(x(1:n))
        x(1:n) = a + (b - a) * x(1:n)

    !-- Dismisses the unused dummy argument warning ----------------------------
        associate(this => this)
        end associate
    end subroutine duniform_builtin

    !------------------------------------------------------------------------------------------------------------------------
end submodule maria_prng_builtin_sub
