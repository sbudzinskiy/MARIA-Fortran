program dgevolume_test
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
use maria_comparison_mod,   only: &
    safe_eq
use maria_utils_mod,        only: &
    are_close
use maria_la_utils_mod,     only: &
    drandsvd
use maria_la_core_mod,      only: &
    dlasrt
use maria_lr_maxvol_mod,    only: &
    dgevolume
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
        integer                :: info, ifoo(1)
        real(WP)               :: foo(1), vol

        vol = dgevolume(-1, -1, -1, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -1)

        vol = dgevolume(1, -1, -1, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -2)

        vol = dgevolume(1, 1, -1, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -3)

        vol = dgevolume(1, 1, 2, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -3)

        vol = dgevolume(1, 1, 1, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -5)
    end subroutine test_sanity

    subroutine test_quick()
        integer                :: info, ifoo(1)
        real(WP)               :: foo(1), vol

        vol = dgevolume(0, 1, 0, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)
        ASSERT(safe_eq(vol, ZERO))

        vol = dgevolume(1, 0, 0, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)
        ASSERT(safe_eq(vol, ZERO))

        vol = dgevolume(1, 1, 0, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)
        ASSERT(safe_eq(vol, ZERO))
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer                :: m, n, r, lda, lwork, liwork, info, ifoo(1), i
        real(WP)               :: vol, volprod, foo(1)
        type(prng)             :: rng
        integer,  allocatable  :: iwork(:)
        real(WP), allocatable  :: A(:), S(:), work(:)
        
        call rng%init(seed, info)
        ASSERT(info == 0)

    !-- m = n = r
        m = 30
        n = m
        ldA = m + 1
        r = m

        call drandsvd(rng, m, n, foo, lda, foo, foo, -1, info)
        lwork = int(foo(1))
        vol = dgevolume(m, n, r, foo, lda, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = ifoo(1)

        allocate(A(lda*n), S(min(m,n)), work(lwork), iwork(liwork))
        
        call rng%duniform(min(m,n), S, 0.5_WP, 1.5_WP, info)
        call dlasrt('d', min(m,n), S, info)
        call drandsvd(rng, m, n, A, lda, S, work, lwork, info)

        volprod = ONE
        do i = 1, r
            volprod = volprod * S(i)
        end do
        vol = dgevolume(m, n, r, A, lda, work, lwork, iwork, liwork, info)
        ASSERT(are_close(vol, volprod, info, rtol=1e-12_WP))
        deallocate(A, S, work, iwork)

    !-- r < min(m,n)
        m = 30
        n = 12
        ldA = m + 1
        r = 5

        call drandsvd(rng, m, n, foo, lda, foo, foo, -1, info)
        lwork = int(foo(1))
        vol = dgevolume(m, n, r, foo, lda, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = ifoo(1)

        allocate(A(lda*n), S(min(m,n)), work(lwork), iwork(liwork))
        
        call rng%duniform(min(m,n), S, 0.5_WP, 1.5_WP, info)
        call dlasrt('d', min(m,n), S, info)
        call drandsvd(rng, m, n, A, lda, S, work, lwork, info)

        volprod = ONE
        do i = 1, r
            volprod = volprod * S(i)
        end do
        vol = dgevolume(m, n, r, A, lda, work, lwork, iwork, liwork, info)
        ASSERT(are_close(vol, volprod, info, rtol=1e-12_WP))
        deallocate(A, S, work, iwork)
    end subroutine test
end program dgevolume_test
