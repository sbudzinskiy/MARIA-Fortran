!----------------------------------------------------------------------------------------------------------------------------
!  MARIA: MAtrix and tensoR Interpolation and Approximation
!----------------------------------------------------------------------------------------------------------------------------
!! Contains the implementation of the [[maria_prng_mkl_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
!> author:  Stanislav Budzinskiy (University of Vienna)
!> version: v0.1
!>
!> Implements the [[maria_prng_mkl_mod(module)]] module.
!----------------------------------------------------------------------------------------------------------------------------
submodule (maria_prng_mkl_mod) maria_prng_mkl_sub
implicit none (type, external)

contains
    !------------------------------------------------------------------------------------------------------------------------

    module subroutine init_mkl &
    (this, seed, info)
    use mkl_vsl,   only: &
        vslNewStream,    &
        VSL_BRNG_MT19937
    !-- Input/output arguments -------------------------------------------------
        class(prng_mkl), intent(in)  :: this        
        integer,         intent(in)  :: seed
        integer,         intent(out) :: info

    !-- Executable section -----------------------------------------------------
        info = vslNewStream(this%stream, VSL_BRNG_MT19937, seed)
    end subroutine init_mkl

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine deinit_mkl &
    (this, info)
    use mkl_vsl,  only: &
        vslDeleteStream
    !-- Input/output arguments -------------------------------------------------
        class(prng_mkl), intent(in)  :: this        
        integer,         intent(out) :: info

    !-- Executable section -----------------------------------------------------
        info = vslDeleteStream(this%stream)
    end subroutine deinit_mkl

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine snormal_mkl &
    (this, n, x, mean, std, info)
    use maria_kinds_mod,       only: &
        WP => SP
    use maria_constants_mod,   only: &
        ZERO => S_ZERO
    use mkl_vsl,               only: &
        vsRngGaussian,               &
        VSL_RNG_METHOD_GAUSSIAN_ICDF
    use maria_argcheck_mod,    only: &
        arg_is_bad,                  &
        BAD_IF_LESS,                 &
        BAD_IF_SAME
    use maria_reports_mod,     only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        class(prng_mkl), intent(in)              :: this        
        integer,         intent(in)              :: n
        real(WP),        intent(out), contiguous :: x(:)
        real(WP),        intent(in)              :: mean
        real(WP),        intent(in)              :: std
        integer,         intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "SNORMAL_MKL"

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
        info = vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, this%stream, n, x, mean, std)
    end subroutine snormal_mkl

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine dnormal_mkl &
    (this, n, x, mean, std, info)
    use maria_kinds_mod,       only: &
        WP => DP
    use maria_constants_mod,   only: &
        ZERO => D_ZERO
    use mkl_vsl,               only: &
        vdRngGaussian,               &
        VSL_RNG_METHOD_GAUSSIAN_ICDF
    use maria_argcheck_mod,    only: &
        arg_is_bad,                  &
        BAD_IF_LESS,                 &
        BAD_IF_SAME
    use maria_reports_mod,     only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        class(prng_mkl), intent(in)              :: this        
        integer,         intent(in)              :: n
        real(WP),        intent(out), contiguous :: x(:)
        real(WP),        intent(in)              :: mean
        real(WP),        intent(in)              :: std
        integer,         intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "DNORMAL_MKL"

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
        info = vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, this%stream, n, x, mean, std)
    end subroutine dnormal_mkl

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine iuniform_mkl &
    (this, n, x, a, b, info)
    use mkl_vsl,             only: &
        viRngUniform,              &
        VSL_RNG_METHOD_UNIFORM_STD
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS,               &
        BAD_IF_SAME
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        class(prng_mkl), intent(in)              :: this        
        integer,         intent(in)              :: n
        integer,         intent(out), contiguous :: x(:)
        integer,         intent(in)              :: a
        integer,         intent(in)              :: b
        integer,         intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "IUNIFORM_MKL"

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
        info = viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, this%stream, n, x, a, b + 1)
    end subroutine iuniform_mkl

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine suniform_mkl &
    (this, n, x, a, b, info)
    use maria_kinds_mod,     only: &
        WP => SP
    use mkl_vsl,             only: &
        vsRngUniform,              &
        VSL_RNG_METHOD_UNIFORM_STD
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS,               &
        BAD_IF_SAME
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        class(prng_mkl), intent(in)              :: this        
        integer,         intent(in)              :: n
        real(WP),        intent(out), contiguous :: x(:)
        real(WP),        intent(in)              :: a
        real(WP),        intent(in)              :: b
        integer,         intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "SUNIFORM_MKL"

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
        info = vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, this%stream, n, x, a, b)
    end subroutine suniform_mkl

    !------------------------------------------------------------------------------------------------------------------------

    module subroutine duniform_mkl &
    (this, n, x, a, b, info)
    use maria_kinds_mod,     only: &
        WP => DP
    use mkl_vsl,             only: &
        vdRngUniform,              &
        VSL_RNG_METHOD_UNIFORM_STD
    use maria_argcheck_mod,  only: &
        arg_is_bad,                &
        BAD_IF_LESS,               &
        BAD_IF_SAME
    use maria_reports_mod,   only: &
        report_bad_arg
    !-- Input/output arguments -------------------------------------------------
        class(prng_mkl), intent(in)              :: this        
        integer,         intent(in)              :: n
        real(WP),        intent(out), contiguous :: x(:)
        real(WP),        intent(in)              :: a
        real(WP),        intent(in)              :: b
        integer,         intent(out)             :: info

    !-- Inner variables --------------------------------------------------------
        character(*), parameter :: SRNAME = "DUNIFORM_MKL"

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
        info = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, this%stream, n, x, a, b)
    end subroutine duniform_mkl

    !------------------------------------------------------------------------------------------------------------------------
end submodule maria_prng_mkl_sub
