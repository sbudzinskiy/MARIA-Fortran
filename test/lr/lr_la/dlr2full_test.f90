program dlr2full_test
#include "maria_assert_mod.h"
#ifdef MARIA_MKL
use maria_prng_mkl_mod,     only: &
    prng => prng_mkl
#else
use maria_prng_builtin_mod, only: &
    prng => prng_builtin
#endif
use maria_kinds_mod,        only: &
    WP => DP
use maria_constants_mod,    only: &
    ONE => D_ONE,                 &
    ZERO => D_ZERO
use maria_utils_mod,        only: &
    are_close
use maria_la_utils_mod,     only: &
    dgeall_const
use maria_la_core_mod,      only: &
    dgenrmc
use maria_lr_la_mod,        only: &
    dlrval,                       &
    dlr2full
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
        real(WP)   :: foo(1)

        call dlr2full(-1, -1, -1, foo, 0, foo, 0, foo, 0, info)
        ASSERT(info == -1)

        call dlr2full(1, -1, -1, foo, 0, foo, 0, foo, 0, info)
        ASSERT(info == -2)

        call dlr2full(1, 1, -1, foo, 0, foo, 0, foo, 0, info)
        ASSERT(info == -3)

        call dlr2full(1, 1, 1, foo, 0, foo, 0, foo, 0, info)
        ASSERT(info == -5)

        call dlr2full(1, 1, 1, foo, 1, foo, 0, foo, 0, info)
        ASSERT(info == -7)

        call dlr2full(1, 1, 1, foo, 1, foo, 1, foo, 0, info)
        ASSERT(info == -9)
    end subroutine test_sanity

    subroutine test_quick()
        integer    :: info
        real(WP)   :: foo(1)

        call dlr2full(0, 1, 1, foo, 1, foo, 1, foo, 1, info)
        ASSERT(info == 0)

        call dlr2full(1, 0, 1, foo, 1, foo, 1, foo, 1, info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer               :: info, m, n, i, j, r, ldU, ldVT, lda
        real(WP)              :: val, nrmc, foo(1)
        type(prng)            :: rng       
        real(WP), allocatable :: U(:), VT(:), A(:)
 
        call rng%init(seed, info)
        ASSERT(info == 0)

    !-- r > 0
        m = 20
        n = 10
        r = 15
        ldU = m + 1
        ldVT = r + 1
        lda = m + 1

        allocate(U(ldU*r), VT(ldVT*n), A(lda*n))
        call rng%dnormal(ldU*r, U, ZERO, ONE, info)
        call rng%dnormal(ldVT*n, VT, ZERO, ONE, info)
        call dlr2full(m, n, r, U, ldu, VT, ldvt, A, lda, info)
        ASSERT(info == 0)
        nrmc = dgenrmc(m, n, A, lda, info)
        do j = 1, n
            do i = 1, m
                val = dlrval(m, n, i, j, r, U, ldu, VT, ldvt, info)
                ASSERT(are_close(val, A(i + (j-1)*lda), info, atol=1e-12_WP * nrmc))
            end do
        end do

        deallocate(U, VT, A)

    !-- r = 0
        m = 20
        n = 10
        r = 0
        ldU = m + 1
        ldVT = r + 1
        lda = m + 1

        allocate(A(lda*n))
        call dlr2full(m, n, r, foo, ldu, foo, ldvt, A, lda, info)
        ASSERT(info == 0)
        ASSERT(dgeall_const(m, n, A, lda, ZERO, ZERO, info))
        deallocate(A)
    end subroutine test
end program dlr2full_test
