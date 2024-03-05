program slrval_test
#include "maria_assert_mod.h"
#ifdef MARIA_MKL
use maria_prng_mkl_mod,     only: &
    prng => prng_mkl
#else
use maria_prng_builtin_mod, only: &
    prng => prng_builtin
#endif
use maria_kinds_mod,        only: &
    WP => SP
use maria_constants_mod,    only: &
    ONE => S_ONE,                 &
    ZERO => S_ZERO
use maria_comparison_mod,   only: &
    safe_eq
use maria_utils_mod,        only: &
    are_close
use maria_la_core_mod,      only: &
    sgemm,                        &
    sgenrmc
use maria_lr_la_mod,        only: &
    slrval
implicit none (type, external)

    character(100) :: seed_str
    integer        :: seed

    call get_command_argument(1, seed_str)
    read(seed_str, *) seed

    call test_sanity()
    call test_quick()
    call test(seed)

contains
    subroutine test_sanity()
        integer    :: info
        real(WP)   :: val, foo(1)

        val = slrval(0, 0, 0, 0, -1, foo, 0, foo, 0, info)
        ASSERT(info == -1)

        val = slrval(1, 0, 0, 0, -1, foo, 0, foo, 0, info)
        ASSERT(info == -2)

        val = slrval(1, 2, 0, 0, -1, foo, 0, foo, 0, info)
        ASSERT(info == -3)

        val = slrval(1, 2, 2, 0, -1, foo, 0, foo, 0, info)
        ASSERT(info == -3)

        val = slrval(2, 1, 1, 0, -1, foo, 0, foo, 0, info)
        ASSERT(info == -4)

        val = slrval(2, 1, 1, 2, -1, foo, 0, foo, 0, info)
        ASSERT(info == -4)

        val = slrval(1, 1, 1, 1, -1, foo, 0, foo, 0, info)
        ASSERT(info == -5)

        val = slrval(1, 1, 1, 1, 1, foo, 0, foo, 0, info)
        ASSERT(info == -7)

        val = slrval(1, 1, 1, 1, 1, foo, 1, foo, 0, info)
        ASSERT(info == -9)
    end subroutine test_sanity

    subroutine test_quick()
        integer    :: info
        real(WP)   :: val, foo(1)

        val = slrval(1, 1, 1, 1, 0, foo, 1, foo, 1, info)
        ASSERT(info == 0)
        ASSERT(safe_eq(val, ZERO))
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer               :: info, m, n, i, j, r, ldU, ldVT
        real(WP)              :: val, nrmc
        type(prng)            :: rng       
        real(WP), allocatable :: U(:), VT(:), A(:)
 
        call rng%init(seed, info)
        ASSERT(info == 0)

        m = 20
        n = 10
        r = 15
        ldU = m + 1
        ldVT = r + 1

        allocate(U(ldU*r), VT(ldVT*n), A(m*n))
        call rng%snormal(ldU*r, U, ZERO, ONE, info)
        call rng%snormal(ldVT*n, VT, ZERO, ONE, info)
        call sgemm('n', 'n', m, n, r, ONE, U, ldu, VT, ldVT, zero, A, m)
        nrmc = sgenrmc(m, n, A, m, info)
        do j = 1, n
            do i = 1, m
                val = slrval(m, n, i, j, r, U, ldu, VT, ldvt, info)
                ASSERT(info == 0)
                ASSERT(are_close(val, A(i + (j-1)*m), info, atol=1e-5_WP * nrmc))
            end do
        end do
    end subroutine test
end program slrval_test
