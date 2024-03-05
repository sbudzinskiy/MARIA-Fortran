program slrdotf_tangent_test
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
    safe_leq,                     &
    safe_less
use maria_utils_mod,        only: &
    are_close
use maria_la_utils_mod,     only: &
    srandort
use maria_la_core_mod,      only: &
    sgemm,                        &
    sgenrmf,                      &
    sgedotf,                      &
    MM => smatmul, slaset
use maria_lr_geom_mod,  only:   &
    slrdotf_tangent, slrproj_tangent
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
        
        dot = slrdotf_tangent(-1, -1, -1, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -1 )

        dot = slrdotf_tangent(2, -1, -1, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -2 )

        dot = slrdotf_tangent(2, 2, -1, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -3 )

        dot = slrdotf_tangent(2, 2, 3, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -3 )

        dot = slrdotf_tangent(2, 2, 1, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -5 )

        dot = slrdotf_tangent(2, 2, 1, foo, 2, foo, 0, foo, 0, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -7 )

        dot = slrdotf_tangent(2, 2, 1, foo, 2, foo, 1, foo, 0, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -9 )

        dot = slrdotf_tangent(2, 2, 1, foo, 2, foo, 1, foo, 1, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -11 )

        dot = slrdotf_tangent(2, 2, 1, foo, 2, foo, 1, foo, 1, foo, 2, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -13 )

        dot = slrdotf_tangent(2, 2, 1, foo, 2, foo, 1, foo, 1, foo, 2, foo, 1, foo, 0, ierr)
        ASSERT( ierr == -15 )
    end subroutine test_sanity

    subroutine test_quick()
        integer                :: ierr
        real(WP)               :: foo(1), dot

        dot = slrdotf_tangent(0, 1, 0, foo, 1, foo, 1, foo, 1, foo, 1, foo, 1, foo, 1, ierr)
        ASSERT( ierr == 0 )

        dot = slrdotf_tangent(1, 0, 0, foo, 1, foo, 1, foo, 1, foo, 1, foo, 1, foo, 1, ierr)
        ASSERT( ierr == 0 )

        dot = slrdotf_tangent(1, 1, 0, foo, 1, foo, 1, foo, 1, foo, 1, foo, 1, foo, 1, ierr)
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

        call srandort(rng, m, r, foo, m, foo, -1, ierr)
        lwork = int(foo(1))
        call srandort(rng, r, n, foo, r, foo, -1, ierr)
        lwork = max(lwork, int(foo(1)))

        allocate( U(m*r), VT(r*n), pU1(ldpu1*r), pVT1(ldpvt1*n), C1(ldc1*r), &
            pU2(ldpu2*r), pVT2(ldpvt2*n), C2(ldc2*r), work(lwork), A(m*n), B(m*n), tmp(m*r) )

        call srandort(rng, m, r, U, m, work, lwork, ierr)
        call srandort(rng, r, n, VT, r, work, lwork, ierr)

        call rng%snormal(m*n, A, ZERO, ONE, ierr)
        X2XY => mm_A_right
        X2YX => mm_A_left
        call slrproj_tangent(m, n, r, U, m, VT, r, X2YX, X2XY, pU1, ldpu1, pVT1, ldpvt1, C1, ldc1, ierr)

        call rng%snormal(m*n, B, ZERO, ONE, ierr)
        X2XY => mm_B_right
        X2YX => mm_B_left
        call slrproj_tangent(m, n, r, U, m, VT, r, X2YX, X2XY, pU2, ldpu2, pVT2, ldpvt2, C2, ldc2, ierr)

        call sgemm('n', 'n', m, n, r, ONE, U, m, pVT1, ldpvt1, ZERO, A, m)
        call sgemm('n', 'n', m, n, r, ONE, pU1, ldpu1, VT, r, ONE, A, m)
        call sgemm('n', 'n', m, r, r, ONE, U, m, C1, ldc1, ZERO, tmp, m)
        call sgemm('n', 'n', m, n, r, ONE, tmp, m, VT, r, ONE, A, m)

        call sgemm('n', 'n', m, n, r, ONE, U, m, pVT2, ldpvt2, ZERO, B, m)
        call sgemm('n', 'n', m, n, r, ONE, pU2, ldpu2, VT, r, ONE, B, m)
        call sgemm('n', 'n', m, r, r, ONE, U, m, C2, ldc2, ZERO, tmp, m)
        call sgemm('n', 'n', m, n, r, ONE, tmp, m, VT, r, ONE, B, m)

        dot1 = slrdotf_tangent(m, n, r, pU1, ldpu1, pVT1, ldpvt1, C1, ldc1, pU2, ldpu2, pVT2, ldpvt2, C2, ldc2, ierr)
        dot2 = sgedotf('n', m, n, A, m, B, m, ierr)

        ASSERT( are_close(dot1, dot2, ierr, rtol=1e-5_WP) )
    end subroutine test

    subroutine mm_A_left &
    (transX, m, n, k, alpha, X, ldx, beta, Y, ldy, info)
    use maria_kinds_mod, only: &
        WP => SP
    use maria_la_core_mod,   only: &
        sgemm
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

        call sgemm('n', transX, m, n, k, alpha, A, lda, X, ldx, beta, Y, ldy)
        info = 0
    end subroutine mm_A_left

    subroutine mm_A_right &
    (transX, m, n, k, alpha, X, ldX, beta, Y, ldy, info)
    use maria_kinds_mod, only: &
        WP => SP
    use maria_la_core_mod,   only: &
        sgemm
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

        call sgemm(transX, 'n', m, n, k, alpha, X, ldx, A, ldA, beta, Y, ldy)
        info = 0
    end subroutine mm_A_right

    subroutine mm_B_left &
    (transX, m, n, k, alpha, X, ldx, beta, Y, ldy, info)
    use maria_kinds_mod, only: &
        WP => SP
    use maria_la_core_mod,   only: &
        sgemm
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

        call sgemm('n', transX, m, n, k, alpha, B, ldb, X, ldx, beta, Y, ldy)
        info = 0
    end subroutine mm_B_left

    subroutine mm_B_right &
    (transX, m, n, k, alpha, X, ldX, beta, Y, ldy, info)
    use maria_kinds_mod, only: &
        WP => SP
    use maria_la_core_mod,   only: &
        sgemm
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

        call sgemm(transX, 'n', m, n, k, alpha, X, ldx, B, ldB, beta, Y, ldy)
        info = 0
    end subroutine mm_B_right

end program slrdotf_tangent_test
