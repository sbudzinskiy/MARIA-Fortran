program slrretr_tangent_test
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
use maria_la_utils_mod,     only: &
    srandort,                       &
    all_close,                      &
    all_const
use maria_la_core_mod,      only: &
    sgemm,                        &
    slacpy,                     &
    sdgmm,                      &
    saxpy,                      &
    sgesdd_q,                   &
    sorcangles,                 &
    MM => smatmul
use maria_lr_geom_mod,  only:   &
    slrproj_tangent,            &
    slrretr_tangent
implicit none (type, external)

    character(100) :: seed_str
    integer        :: seed, lda
    real(WP), allocatable :: A(:)

    call get_command_argument(1, seed_str)
    read(seed_str, *) seed

    call test_sanity()
    call test_quick()
    call test(seed)

contains
    subroutine test_sanity()
        integer                 :: ierr, ifoo(1)
        real(WP)                :: foo(1)

        call slrretr_tangent(-1, -1, -1, foo, foo, 0, foo, 0, ONE, foo, 0, foo, 0, foo, 0, foo, 0, ifoo, 0, ierr)
        ASSERT( ierr == -1 )

        call slrretr_tangent(2, -1, -1, foo, foo, 0, foo, 0, ONE, foo, 0, foo, 0, foo, 0, foo, 0, ifoo, 0, ierr)
        ASSERT( ierr == -2 )

        call slrretr_tangent(2, 2, -1, foo, foo, 0, foo, 0, ONE, foo, 0, foo, 0, foo, 0, foo, 0, ifoo, 0, ierr)
        ASSERT( ierr == -3 )

        call slrretr_tangent(2, 2, 3, foo, foo, 0, foo, 0, ONE, foo, 0, foo, 0, foo, 0, foo, 0, ifoo, 0, ierr)
        ASSERT( ierr == -3 )

        call slrretr_tangent(2, 2, 1, foo, foo, 0, foo, 0, ONE, foo, 0, foo, 0, foo, 0, foo, 0, ifoo, 0, ierr)
        ASSERT( ierr == -6 )

        call slrretr_tangent(2, 2, 1, foo, foo, 2, foo, 0, ONE, foo, 0, foo, 0, foo, 0, foo, 0, ifoo, 0, ierr)
        ASSERT( ierr == -8 )

        call slrretr_tangent(2, 2, 1, foo, foo, 2, foo, 1, ONE, foo, 0, foo, 0, foo, 0, foo, 0, ifoo, 0, ierr)
        ASSERT( ierr == -11 )

        call slrretr_tangent(2, 2, 1, foo, foo, 2, foo, 1, ONE, foo, 2, foo, 0, foo, 0, foo, 0, ifoo, 0, ierr)
        ASSERT( ierr == -13 )

        call slrretr_tangent(2, 2, 1, foo, foo, 2, foo, 1, ONE, foo, 2, foo, 1, foo, 0, foo, 0, ifoo, 0, ierr)
        ASSERT( ierr == -15 )
    end subroutine test_sanity

    subroutine test_quick()
        integer                :: ierr, ifoo(1)
        real(WP)               :: foo(1)

        call slrretr_tangent(0, 1, 0, foo, foo, 1, foo, 1, ONE, foo, 1, foo, 1, foo, 1, foo, 0, ifoo, 0, ierr)
        ASSERT( ierr == 0 )

        call slrretr_tangent(1, 0, 0, foo, foo, 1, foo, 1, ONE, foo, 1, foo, 1, foo, 1, foo, 0, ifoo, 0, ierr)
        ASSERT( ierr == 0 )

        call slrretr_tangent(1, 1, 0, foo, foo, 1, foo, 1, ONE, foo, 1, foo, 1, foo, 1, foo, 0, ifoo, 0, ierr)
        ASSERT( ierr == 0 )
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer                :: m, n, r, ldu, ldvt, ldpu, ldpvt, ldc, lwork, liwork, ierr, ifoo(1)
        real(WP)               :: foo(1), alpha
        type(prng)             :: rng
        procedure(MM), pointer :: fun, B2AB, B2BA
        integer,    allocatable :: iwork(:)
        real(WP), allocatable  :: U(:), VT(:), S(:), pU(:), pVT(:), C(:), work(:), tmp(:), UU(:), VVT(:), SS(:), X(:), cosines(:)
        
        call rng%init(seed, ierr)
        ASSERT( ierr == 0 )

        B2AB => mm_A_left
        B2BA => mm_A_right

    !-- m > 2r, n > 2r
        m = 20
        n = 25
        r = 8
        lda = m
        ldu = m + 1
        ldvt = r + 1
        ldpu = m + 1
        ldpvt = r + 1
        ldc = r + 1
        alpha = 2 * ONE

        call srandort(rng, m, r, foo, ldu, foo, -1, ierr)
        lwork = int(foo(1))
        call srandort(rng, r, n, foo, ldvt, foo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        call slrretr_tangent(m, n, r, foo, foo, ldu, foo, ldvt, alpha, foo, ldpu, foo, ldpvt, foo, ldc, foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = ifoo(1)
        call sgesdd_q('s', m, n, foo, m, foo, foo, m, foo, min(m,n), foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
        call sorcangles('c', m, r, foo, ldu, foo, m, foo, foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
        call sorcangles('r', r, n, foo, ldvt, foo, r, foo, foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        allocate( A(m*n), X(m*n), U(ldu*r), VT(ldvt*n), S(r), pU(ldpu*r), pVT(ldpvt*n), C(ldc*r), work(lwork), tmp(n*r), & 
            UU(m*min(m,n)), VVT(min(m,n)*n), SS(min(m,n)), cosines(r), iwork(liwork) )

        ! Generate X as SVD, form it explicitly
        call srandort(rng, m, r, U, ldu, work, lwork, ierr)
        call srandort(rng, r, n, VT, ldvt, work, lwork, ierr)
        call rng%suniform(r, S, ONE, 2*ONE, ierr)
        call slacpy('a', r, n, VT, ldvt, tmp, r)
        call sdgmm('l', r, n, tmp, r, S, 1, ierr)
        call sgemm('n', 'n', m, n, r, ONE, U, ldu, tmp, r, ZERO, X, m)

        ! Generate random A, project it onto tangent space at X, form the projected A explicitly
        call rng%snormal(m*n, A, ZERO, ONE, ierr)
        call slrproj_tangent(m, n, r, U, ldu, VT, ldvt, B2AB, B2BA, pU, ldpu, pVT, ldpvt, C, ldc, ierr)
        call sgemm('n', 'n', m, n, r, ONE, U, ldu, pVT, ldpVT, ZERO, A, m)
        call sgemm('n', 'n', m, n, r, ONE, pU, ldpu, VT, ldvt, ONE, A, m)
        call sgemm('n', 'n', r, n, r, ONE, C, ldc, VT, ldvt, ZERO, tmp, r)
        call sgemm('n', 'n', m, n, r, ONE, U, ldu, tmp, r, ONE, A, m)

        ! Do retraction of X + alpha*A
        call slrretr_tangent(m, n, r, S, U, ldu, VT, ldvt, &
            alpha, pU, ldpu, pVT, ldpvt, C, ldc, work, lwork, iwork, liwork, ierr)
        ASSERT( ierr == 0 )

        ! Add X + alpha*A explicitly, truncated SVD
        call saxpy(m*n, alpha, A, 1, X, 1)
        call sgesdd_q('s', m, n, X, m, SS, UU, m, VVT, min(m,n), work, lwork, iwork, liwork, ierr)
        
        ! Compare
        ASSERT( all_close(r, S, 1, SS, 1, ierr, atol=1e-5_WP) )
        call sorcangles('c', m, r, U, ldu, UU, m, cosines, work, lwork, iwork, liwork, ierr)
        ASSERT( all_const(r, cosines, 1, ONE, ierr, atol=1e-5_WP) )
        call sorcangles('r', r, n, VT, ldvt, VVT, min(m,n), cosines, work, lwork, iwork, liwork, ierr)
        ASSERT( all_const(r, cosines, 1, ONE, ierr, atol=1e-5_WP) )

        deallocate( A, X, U, VT, S, pU, pVT, C, work, tmp, UU, VVT, SS, cosines, iwork )

    !-- m < 2r, n > 2r
        m = 12
        n = 25
        r = 8
        lda = m
        ldu = m + 1
        ldvt = r + 1
        ldpu = m + 1
        ldpvt = r + 1
        ldc = r + 1
        alpha = 2 * ONE

        call srandort(rng, m, r, foo, ldu, foo, -1, ierr)
        lwork = int(foo(1))
        call srandort(rng, r, n, foo, ldvt, foo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        call slrretr_tangent(m, n, r, foo, foo, ldu, foo, ldvt, alpha, foo, ldpu, foo, ldpvt, foo, ldc, foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = ifoo(1)
        call sgesdd_q('s', m, n, foo, m, foo, foo, m, foo, min(m,n), foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
        call sorcangles('c', m, r, foo, ldu, foo, m, foo, foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
        call sorcangles('r', r, n, foo, ldvt, foo, r, foo, foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        allocate( A(m*n), X(m*n), U(ldu*r), VT(ldvt*n), S(r), pU(ldpu*r), pVT(ldpvt*n), C(ldc*r), work(lwork), tmp(n*r), & 
            UU(m*min(m,n)), VVT(min(m,n)*n), SS(min(m,n)), cosines(r), iwork(liwork) )

        ! Generate X as SVD, form it explicitly
        call srandort(rng, m, r, U, ldu, work, lwork, ierr)
        call srandort(rng, r, n, VT, ldvt, work, lwork, ierr)
        call rng%suniform(r, S, ONE, 2*ONE, ierr)
        call slacpy('a', r, n, VT, ldvt, tmp, r)
        call sdgmm('l', r, n, tmp, r, S, 1, ierr)
        call sgemm('n', 'n', m, n, r, ONE, U, ldu, tmp, r, ZERO, X, m)

        ! Generate random A, project it onto tangent space at X, form the projected A explicitly
        call rng%snormal(m*n, A, ZERO, ONE, ierr)
        call slrproj_tangent(m, n, r, U, ldu, VT, ldvt, B2AB, B2BA, pU, ldpu, pVT, ldpvt, C, ldc, ierr)
        call sgemm('n', 'n', m, n, r, ONE, U, ldu, pVT, ldpVT, ZERO, A, m)
        call sgemm('n', 'n', m, n, r, ONE, pU, ldpu, VT, ldvt, ONE, A, m)
        call sgemm('n', 'n', r, n, r, ONE, C, ldc, VT, ldvt, ZERO, tmp, r)
        call sgemm('n', 'n', m, n, r, ONE, U, ldu, tmp, r, ONE, A, m)

        ! Do retraction of X + alpha*A
        call slrretr_tangent(m, n, r, S, U, ldu, VT, ldvt, &
            alpha, pU, ldpu, pVT, ldpvt, C, ldc, work, lwork, iwork, liwork, ierr)
        ASSERT( ierr == 0 )

        ! Add X + alpha*A explicitly, truncated SVD
        call saxpy(m*n, alpha, A, 1, X, 1)
        call sgesdd_q('s', m, n, X, m, SS, UU, m, VVT, min(m,n), work, lwork, iwork, liwork, ierr)
        
        ! Compare
        ASSERT( all_close(r, S, 1, SS, 1, ierr, atol=1e-5_WP) )
        call sorcangles('c', m, r, U, ldu, UU, m, cosines, work, lwork, iwork, liwork, ierr)
        ASSERT( all_const(r, cosines, 1, ONE, ierr, atol=1e-5_WP) )
        call sorcangles('r', r, n, VT, ldvt, VVT, min(m,n), cosines, work, lwork, iwork, liwork, ierr)
        ASSERT( all_const(r, cosines, 1, ONE, ierr, atol=1e-5_WP) )

        deallocate( A, X, U, VT, S, pU, pVT, C, work, tmp, UU, VVT, SS, cosines, iwork )

    !-- m > 2r, n < 2r
        m = 20
        n = 13
        r = 8
        lda = m
        ldu = m + 1
        ldvt = r + 1
        ldpu = m + 1
        ldpvt = r + 1
        ldc = r + 1
        alpha = 2 * ONE

        call srandort(rng, m, r, foo, ldu, foo, -1, ierr)
        lwork = int(foo(1))
        call srandort(rng, r, n, foo, ldvt, foo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        call slrretr_tangent(m, n, r, foo, foo, ldu, foo, ldvt, alpha, foo, ldpu, foo, ldpvt, foo, ldc, foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = ifoo(1)
        call sgesdd_q('s', m, n, foo, m, foo, foo, m, foo, min(m,n), foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
        call sorcangles('c', m, r, foo, ldu, foo, m, foo, foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
        call sorcangles('r', r, n, foo, ldvt, foo, r, foo, foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        allocate( A(m*n), X(m*n), U(ldu*r), VT(ldvt*n), S(r), pU(ldpu*r), pVT(ldpvt*n), C(ldc*r), work(lwork), tmp(n*r), & 
            UU(m*min(m,n)), VVT(min(m,n)*n), SS(min(m,n)), cosines(r), iwork(liwork) )

        ! Generate X as SVD, form it explicitly
        call srandort(rng, m, r, U, ldu, work, lwork, ierr)
        call srandort(rng, r, n, VT, ldvt, work, lwork, ierr)
        call rng%suniform(r, S, ONE, 2*ONE, ierr)
        call slacpy('a', r, n, VT, ldvt, tmp, r)
        call sdgmm('l', r, n, tmp, r, S, 1, ierr)
        call sgemm('n', 'n', m, n, r, ONE, U, ldu, tmp, r, ZERO, X, m)

        ! Generate random A, project it onto tangent space at X, form the projected A explicitly
        call rng%snormal(m*n, A, ZERO, ONE, ierr)
        call slrproj_tangent(m, n, r, U, ldu, VT, ldvt, B2AB, B2BA, pU, ldpu, pVT, ldpvt, C, ldc, ierr)
        call sgemm('n', 'n', m, n, r, ONE, U, ldu, pVT, ldpVT, ZERO, A, m)
        call sgemm('n', 'n', m, n, r, ONE, pU, ldpu, VT, ldvt, ONE, A, m)
        call sgemm('n', 'n', r, n, r, ONE, C, ldc, VT, ldvt, ZERO, tmp, r)
        call sgemm('n', 'n', m, n, r, ONE, U, ldu, tmp, r, ONE, A, m)

        ! Do retraction of X + alpha*A
        call slrretr_tangent(m, n, r, S, U, ldu, VT, ldvt, &
            alpha, pU, ldpu, pVT, ldpvt, C, ldc, work, lwork, iwork, liwork, ierr)
        ASSERT( ierr == 0 )

        ! Add X + alpha*A explicitly, truncated SVD
        call saxpy(m*n, alpha, A, 1, X, 1)
        call sgesdd_q('s', m, n, X, m, SS, UU, m, VVT, min(m,n), work, lwork, iwork, liwork, ierr)
        
        ! Compare
        ASSERT( all_close(r, S, 1, SS, 1, ierr, atol=1e-5_WP) )
        call sorcangles('c', m, r, U, ldu, UU, m, cosines, work, lwork, iwork, liwork, ierr)
        ASSERT( all_const(r, cosines, 1, ONE, ierr, atol=1e-5_WP) )
        call sorcangles('r', r, n, VT, ldvt, VVT, min(m,n), cosines, work, lwork, iwork, liwork, ierr)
        ASSERT( all_const(r, cosines, 1, ONE, ierr, atol=1e-5_WP) )

        deallocate( A, X, U, VT, S, pU, pVT, C, work, tmp, UU, VVT, SS, cosines, iwork )

    !-- m < 2r, n < 2r
        m = 11
        n = 13
        r = 8
        lda = m
        ldu = m + 1
        ldvt = r + 1
        ldpu = m + 1
        ldpvt = r + 1
        ldc = r + 1
        alpha = 2 * ONE

        call srandort(rng, m, r, foo, ldu, foo, -1, ierr)
        lwork = int(foo(1))
        call srandort(rng, r, n, foo, ldvt, foo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        call slrretr_tangent(m, n, r, foo, foo, ldu, foo, ldvt, alpha, foo, ldpu, foo, ldpvt, foo, ldc, foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = ifoo(1)
        call sgesdd_q('s', m, n, foo, m, foo, foo, m, foo, min(m,n), foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
        call sorcangles('c', m, r, foo, ldu, foo, m, foo, foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
        call sorcangles('r', r, n, foo, ldvt, foo, r, foo, foo, -1, ifoo, -1, ierr)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        allocate( A(m*n), X(m*n), U(ldu*r), VT(ldvt*n), S(r), pU(ldpu*r), pVT(ldpvt*n), C(ldc*r), work(lwork), tmp(n*r), & 
            UU(m*min(m,n)), VVT(min(m,n)*n), SS(min(m,n)), cosines(r), iwork(liwork) )

        ! Generate X as SVD, form it explicitly
        call srandort(rng, m, r, U, ldu, work, lwork, ierr)
        call srandort(rng, r, n, VT, ldvt, work, lwork, ierr)
        call rng%suniform(r, S, ONE, 2*ONE, ierr)
        call slacpy('a', r, n, VT, ldvt, tmp, r)
        call sdgmm('l', r, n, tmp, r, S, 1, ierr)
        call sgemm('n', 'n', m, n, r, ONE, U, ldu, tmp, r, ZERO, X, m)

        ! Generate random A, project it onto tangent space at X, form the projected A explicitly
        call rng%snormal(m*n, A, ZERO, ONE, ierr)
        call slrproj_tangent(m, n, r, U, ldu, VT, ldvt, B2AB, B2BA, pU, ldpu, pVT, ldpvt, C, ldc, ierr)
        call sgemm('n', 'n', m, n, r, ONE, U, ldu, pVT, ldpVT, ZERO, A, m)
        call sgemm('n', 'n', m, n, r, ONE, pU, ldpu, VT, ldvt, ONE, A, m)
        call sgemm('n', 'n', r, n, r, ONE, C, ldc, VT, ldvt, ZERO, tmp, r)
        call sgemm('n', 'n', m, n, r, ONE, U, ldu, tmp, r, ONE, A, m)

        ! Do retraction of X + alpha*A
        call slrretr_tangent(m, n, r, S, U, ldu, VT, ldvt, &
            alpha, pU, ldpu, pVT, ldpvt, C, ldc, work, lwork, iwork, liwork, ierr)
        ASSERT( ierr == 0 )

        ! Add X + alpha*A explicitly, truncated SVD
        call saxpy(m*n, alpha, A, 1, X, 1)
        call sgesdd_q('s', m, n, X, m, SS, UU, m, VVT, min(m,n), work, lwork, iwork, liwork, ierr)
        
        ! Compare
        ASSERT( all_close(r, S, 1, SS, 1, ierr, atol=1e-5_WP) )
        call sorcangles('c', m, r, U, ldu, UU, m, cosines, work, lwork, iwork, liwork, ierr)
        ASSERT( all_const(r, cosines, 1, ONE, ierr, atol=1e-5_WP) )
        call sorcangles('r', r, n, VT, ldvt, VVT, min(m,n), cosines, work, lwork, iwork, liwork, ierr)
        ASSERT( all_const(r, cosines, 1, ONE, ierr, atol=1e-5_WP) )
    end subroutine test

    subroutine mm_A_left &
    (transB, m, n, k, alpha, B, ldb, beta, C, ldc, info)
    use maria_kinds_mod, only: &
        WP => SP
    use maria_la_core_mod,   only: &
        sgemm
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

        call sgemm('n', transB, m, n, k, alpha, A, lda, B, ldB, beta, C, ldc)
        info = 0
    end subroutine mm_A_left

    subroutine mm_A_right &
    (transB, m, n, k, alpha, B, ldb, beta, C, ldc, info)
    use maria_kinds_mod, only: &
        WP => SP
    use maria_la_core_mod,   only: &
        sgemm
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

        call sgemm(transB, 'n', m, n, k, alpha, B, ldB, A, lda, beta, C, ldc)
        info = 0
    end subroutine mm_A_right

end program slrretr_tangent_test
