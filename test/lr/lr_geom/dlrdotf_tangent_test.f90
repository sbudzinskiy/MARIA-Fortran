program dlrdotf_tangent_test
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
    safe_leq,                     &
    safe_less
use maria_utils_mod,        only: &
    are_close
use maria_la_utils_mod,     only: &
    drandort
use maria_la_core_mod,      only: &
    dgemm,                        &
    dgenrmf,                      &
    dgedotf,                      &
    MM => dmatmul, dlaset
use maria_lr_geom_mod,  only:   &
    dlrdotf_tangent, dlrproj_tangent
implicit none (type, external)

    character(100) :: seed_str
    integer        :: seed, lda, ldb

    real(WP), allocatable :: A(:), B(:)

    call get_command_argument(1, seed_str)
    read(seed_str, *) seed

    call test_sanity()
    call test_quick()
    call test(seed)

contains
    subroutine test_sanity()
        integer                 :: ierr
        real(WP)                :: foo(1), dot
        
        dot = dlrdotf_tangent(-1, -1, -1, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -1 )

        dot = dlrdotf_tangent(2, -1, -1, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -2 )

        dot = dlrdotf_tangent(2, 2, -1, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -3 )

        dot = dlrdotf_tangent(2, 2, 3, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -3 )

        dot = dlrdotf_tangent(2, 2, 1, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -5 )

        dot = dlrdotf_tangent(2, 2, 1, foo, 2, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -7 )

        dot = dlrdotf_tangent(2, 2, 1, foo, 2, foo, 1, foo, 0, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -9 )

        dot = dlrdotf_tangent(2, 2, 1, foo, 2, foo, 1, foo, 1, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -11 )

        dot = dlrdotf_tangent(2, 2, 1, foo, 2, foo, 1, foo, 1, foo, 2, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -13 )

        dot = dlrdotf_tangent(2, 2, 1, foo, 2, foo, 1, foo, 1, foo, 2, foo, 1, foo, 0, ierr)
        ASSERT( ierr == -15 )
    end subroutine test_sanity

    subroutine test_quick()
        integer                :: ierr
        real(WP)               :: foo(1), dot

        dot = dlrdotf_tangent(0, 1, 0, foo, 1, foo, 1, foo, 1, foo, 1, foo, 1, foo, 1, ierr)
        ASSERT( ierr == 0 )

        dot = dlrdotf_tangent(1, 0, 0, foo, 1, foo, 1, foo, 1, foo, 1, foo, 1, foo, 1, ierr)
        ASSERT( ierr == 0 )

        dot = dlrdotf_tangent(1, 1, 0, foo, 1, foo, 1, foo, 1, foo, 1, foo, 1, foo, 1, ierr)
        ASSERT( ierr == 0 )
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer                :: m, n, r, ldpu1, ldpvt1, ldc1, ldpu2, ldpvt2, ldc2, lwork, ierr
        real(WP)               :: dot1, dot2, foo(1)
        type(prng)             :: rng
        procedure(MM), pointer   :: X2XY, X2YX
        real(WP), allocatable  :: U(:), VT(:), pU1(:), pVT1(:), C1(:), pU2(:), pVT2(:), C2(:), work(:), tmp(:)
        
        call rng%init(seed, ierr)
        ASSERT( ierr == 0 )

        m = 20
        n = 25
        r = 15
        lda = m + 1
        ldb = m + 1
        ldpu1 = m + 1
        ldpvt1 = r + 1
        ldc1 = r + 1
        ldpu2 = m + 1
        ldpvt2 = r + 1
        ldc2 = r + 1

        call drandort(rng, m, r, foo, m, foo, -1, ierr)
        lwork = int(foo(1))
        call drandort(rng, r, n, foo, r, foo, -1, ierr)
        lwork = max(lwork, int(foo(1)))

        allocate( U(m*r), VT(r*n), pU1(ldpu1*r), pVT1(ldpvt1*n), C1(ldc1*r), &
            pU2(ldpu2*r), pVT2(ldpvt2*n), C2(ldc2*r), work(lwork), A(m*n), B(m*n), tmp(m*r) )

        call drandort(rng, m, r, U, m, work, lwork, ierr)
        call drandort(rng, r, n, VT, r, work, lwork, ierr)

        call rng%dnormal(m*n, A, ZERO, ONE, ierr)
        X2XY => mm_A_right
        X2YX => mm_A_left
        call dlrproj_tangent(m, n, r, U, m, VT, r, X2YX, X2XY, pU1, ldpu1, pVT1, ldpvt1, C1, ldc1, ierr)

        call rng%dnormal(m*n, B, ZERO, ONE, ierr)
        X2XY => mm_B_right
        X2YX => mm_B_left
        call dlrproj_tangent(m, n, r, U, m, VT, r, X2YX, X2XY, pU2, ldpu2, pVT2, ldpvt2, C2, ldc2, ierr)

        call dgemm('n', 'n', m, n, r, ONE, U, m, pVT1, ldpvt1, ZERO, A, m)
        call dgemm('n', 'n', m, n, r, ONE, pU1, ldpu1, VT, r, ONE, A, m)
        call dgemm('n', 'n', m, r, r, ONE, U, m, C1, ldc1, ZERO, tmp, m)
        call dgemm('n', 'n', m, n, r, ONE, tmp, m, VT, r, ONE, A, m)

        call dgemm('n', 'n', m, n, r, ONE, U, m, pVT2, ldpvt2, ZERO, B, m)
        call dgemm('n', 'n', m, n, r, ONE, pU2, ldpu2, VT, r, ONE, B, m)
        call dgemm('n', 'n', m, r, r, ONE, U, m, C2, ldc2, ZERO, tmp, m)
        call dgemm('n', 'n', m, n, r, ONE, tmp, m, VT, r, ONE, B, m)

        dot1 = dlrdotf_tangent(m, n, r, pU1, ldpu1, pVT1, ldpvt1, C1, ldc1, pU2, ldpu2, pVT2, ldpvt2, C2, ldc2, ierr)
        dot2 = dgedotf('n', m, n, A, m, B, m, ierr)

        ASSERT( are_close(dot1, dot2, ierr, rtol=1e-12_WP) )
    end subroutine test

    subroutine mm_A_left &
    (transX, m, n, k, alpha, X, ldx, beta, Y, ldy, info)
    use maria_kinds_mod, only: &
        WP => DP
    use maria_la_core_mod,   only: &
        dgemm
    implicit none
        character(1), intent(in)                :: transX
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        integer,      intent(in)                :: k
        real(WP),     intent(in)                :: alpha
        real(WP),     intent(in),    contiguous :: X(:)
        integer,      intent(in)                :: ldX
        real(WP),     intent(in)                :: beta
        real(WP),     intent(inout), contiguous :: Y(:)
        integer,      intent(in)                :: ldY
        integer,      intent(out)               :: info

        call dgemm('n', transX, m, n, k, alpha, A, lda, X, ldx, beta, Y, ldy)
        info = 0
    end subroutine mm_A_left

    subroutine mm_A_right &
    (transX, m, n, k, alpha, X, ldX, beta, Y, ldy, info)
    use maria_kinds_mod, only: &
        WP => DP
    use maria_la_core_mod,   only: &
        dgemm
    implicit none
        character(1), intent(in)                :: transX
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        integer,      intent(in)                :: k
        real(WP),     intent(in)                :: alpha
        real(WP),     intent(in),    contiguous :: X(:)
        integer,      intent(in)                :: ldX
        real(WP),     intent(in)                :: beta
        real(WP),     intent(inout), contiguous :: Y(:)
        integer,      intent(in)                :: ldy
        integer,      intent(out)               :: info

        call dgemm(transX, 'n', m, n, k, alpha, X, ldx, A, ldA, beta, Y, ldy)
        info = 0
    end subroutine mm_A_right

    subroutine mm_B_left &
    (transX, m, n, k, alpha, X, ldx, beta, Y, ldy, info)
    use maria_kinds_mod, only: &
        WP => DP
    use maria_la_core_mod,   only: &
        dgemm
    implicit none
        character(1), intent(in)                :: transX
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        integer,      intent(in)                :: k
        real(WP),     intent(in)                :: alpha
        real(WP),     intent(in),    contiguous :: X(:)
        integer,      intent(in)                :: ldX
        real(WP),     intent(in)                :: beta
        real(WP),     intent(inout), contiguous :: Y(:)
        integer,      intent(in)                :: ldY
        integer,      intent(out)               :: info

        call dgemm('n', transX, m, n, k, alpha, B, ldb, X, ldx, beta, Y, ldy)
        info = 0
    end subroutine mm_B_left

    subroutine mm_B_right &
    (transX, m, n, k, alpha, X, ldX, beta, Y, ldy, info)
    use maria_kinds_mod, only: &
        WP => DP
    use maria_la_core_mod,   only: &
        dgemm
    implicit none
        character(1), intent(in)                :: transX
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        integer,      intent(in)                :: k
        real(WP),     intent(in)                :: alpha
        real(WP),     intent(in),    contiguous :: X(:)
        integer,      intent(in)                :: ldX
        real(WP),     intent(in)                :: beta
        real(WP),     intent(inout), contiguous :: Y(:)
        integer,      intent(in)                :: ldy
        integer,      intent(out)               :: info

        call dgemm(transX, 'n', m, n, k, alpha, X, ldx, B, ldB, beta, Y, ldy)
        info = 0
    end subroutine mm_B_right

end program dlrdotf_tangent_test
