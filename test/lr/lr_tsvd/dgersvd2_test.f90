program dgersvd2_test
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
use maria_access_matrix_mod, only: &
    dmatval_ballistic
use maria_la_core_mod,      only: &
    dgemm,                        &
    dgenrmf,                      &
    ddgmm,                        &
    dgesdd_q,                     &
    MM => dmatmul
use maria_lr_tsvd_mod,      only: &
    dchop,                        &
    dgersvd2
implicit none (type, external)

    character(100) :: seed_str
    integer        :: seed
    real(WP), allocatable :: gauss_left(:)

    call get_command_argument(1, seed_str)
    read(seed_str, *) seed

    call test_sanity()
    call test_quick()
    call test(seed)

contains
    subroutine test_sanity()
        integer                :: info, ifoo(1)
        real(WP)               :: foo(1)
        procedure(MM), pointer :: fun

        call dgersvd2(-1, -1, -1, foo, 0, fun, -1, foo, 0, foo, foo, 0, foo, 0, foo, -1, ifoo, -1, info)
        ASSERT(info == -1)

        call dgersvd2(0, -1, -1, foo, 0, fun, -1, foo, 0, foo, foo, 0, foo, 0, foo, -1, ifoo, -1, info)
        ASSERT(info == -2)

        call dgersvd2(1, 1, -1, foo, 0, fun, -1, foo, 0, foo, foo, 0, foo, 0, foo, -1, ifoo, -1, info)
        ASSERT(info == -3)

        call dgersvd2(1, 1, 2, foo, 0, fun, -1, foo, 0, foo, foo, 0, foo, 0, foo, -1, ifoo, -1, info)
        ASSERT(info == -3)

        call dgersvd2(1, 1, 0, foo, 0, fun, -1, foo, 0, foo, foo, 0, foo, 0, foo, -1, ifoo, -1, info)
        ASSERT(info == -5)

        call dgersvd2(1, 1, 1, foo, 1, fun, 0, foo, 0, foo, foo, 0, foo, 0, foo, -1, ifoo, -1, info)
        ASSERT(info == -7)

        call dgersvd2(1, 1, 1, foo, 1, fun, 2, foo, 0, foo, foo, 0, foo, 0, foo, -1, ifoo, -1, info)
        ASSERT(info == -7)

        call dgersvd2(1, 1, 1, foo, 1, fun, 1, foo, 0, foo, foo, 0, foo, 0, foo, -1, ifoo, -1, info)
        ASSERT(info == -9)

        call dgersvd2(1, 1, 1, foo, 1, fun, 1, foo, 1, foo, foo, 0, foo, 0, foo, -1, ifoo, -1, info)
        ASSERT(info == -12)

        call dgersvd2(1, 1, 1, foo, 1, fun, 1, foo, 1, foo, foo, 1, foo, 0, foo, -1, ifoo, -1, info)
        ASSERT(info == -14)
    end subroutine test_sanity

    subroutine test_quick()
        integer                :: info, ifoo(1)
        real(WP)               :: foo(1)
        procedure(MM), pointer :: fun

        call dgersvd2(0, 1, 0, foo, 1, fun, 0, foo, 1, foo, foo, 1, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)

        call dgersvd2(1, 0, 0, foo, 1, fun, 1, foo, 1, foo, foo, 1, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer                :: i, j, m, n, lda, r, kc, ldcsk, kr, ldrsk, ldu, ldvt, lwork, liwork, info, ifoo(1)
        real(WP)               :: err_rsvd0, err_tsvd, foo(1)
        type(prng)             :: rng
        procedure(MM), pointer :: fun
        integer,  allocatable  :: iwork(:)
        real(WP), allocatable  :: A(:), gauss_right(:), csk(:), rsk(:), S(:), U(:), VT(:), work(:)
        
        call rng%init(seed, info)
        ASSERT(info == 0)

        fun => sketch_rowspace

        m = 30
        n = 40
        ldA = m + 1
        r = 5
        kc = r + 5
        ldcsk = m + 1
        kr = kc + 2
        ldrsk = kr + 1
        ldu = m + 1
        ldvt = kc + 1

        call dgersvd2(m, n, kc, foo, ldcsk, fun, kr, foo, ldrsk, foo, foo, ldu, foo, ldvt, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)
        call dgesdd_q('n', m, n, foo, ldA, foo, foo, ldu, foo, ldvt, foo, -1, ifoo, -1, info)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        allocate(A(lda*n), gauss_right(n*kc), csk(ldcsk*kc), gauss_left(kr*m), rsk(ldrsk*n), &
            S(min(m,n)), U(ldu*kc), VT(ldvt*n), work(lwork), iwork(liwork))
        
        do i = 1, m
            do j = 1, n
                A(i + (j-1)*lda) = dmatval_ballistic(m, n, i, j, info)
            end do
        end do
        call rng%dnormal(n*kc, gauss_right, ZERO, ONE, info)
        call rng%dnormal(kr*m, gauss_left, ZERO, ONE, info)
        call dgemm('n', 'n', m, kc, n, ONE, A, lda, gauss_right, n, ZERO, csk, ldcsk)
        call dgemm('n', 'n', kr, n, m, ONE, gauss_left, kr, A, lda, ZERO, rsk, ldrsk)

        call dgersvd2(m, n, kc, csk, ldcsk, fun, kr, rsk, ldrsk, S, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
        call ddgmm('r', m, kc, U, ldu, S, 1, info)
        call dgemm('n', 'n', m, n, r, -ONE, U, ldu, VT, ldvt, ONE, A, lda)
        err_rsvd0 = dgenrmf(m, n, A, lda, info)

        do i = 1, m
            do j = 1, n
                A(i + (j-1)*lda) = dmatval_ballistic(m, n, i, j, info)
            end do
        end do
        call dgesdd_q('n', m, n, A, lda, S, foo, ldu, foo, ldvt, work, lwork, iwork, liwork, info)
        ifoo(1) = dchop(min(m,n), S, info, maxr=r, aerrf=err_tsvd)
        ASSERT(safe_leq(err_tsvd, err_rsvd0))
        ASSERT(safe_leq(err_rsvd0, 2.0_WP * err_tsvd))
        deallocate(A, gauss_right, csk, gauss_left, rsk, S, U, VT, work, iwork)
    end subroutine test

    subroutine sketch_rowspace &
    (transB, m, n, k, alpha, B, ldb, beta, C, ldc, info)
    use maria_kinds_mod, only: &
        WP => DP
    use maria_la_core_mod,   only: &
        dgemm
    implicit none
        character(1), intent(in)                :: transB
        integer,      intent(in)                :: m
        integer,      intent(in)                :: n
        integer,      intent(in)                :: k
        real(WP),     intent(in)                :: alpha
        real(WP),     intent(in),    contiguous :: B(:)
        integer,      intent(in)                :: ldB
        real(WP),     intent(in)                :: beta
        real(WP),     intent(inout), contiguous :: C(:)
        integer,      intent(in)                :: ldC
        integer,      intent(out)               :: info

        call dgemm('n', transB, m, n, k, alpha, gauss_left, m, B, ldB, beta, C, ldc)
        info = 0
    end subroutine sketch_rowspace

end program dgersvd2_test
