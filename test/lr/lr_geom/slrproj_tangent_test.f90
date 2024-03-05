program slrproj_tangent_test
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
    srandort
use maria_la_core_mod,      only: &
    sgemm,                        &
    sgenrmf,                      &
    slaset,                       &
    MM => smatmul
use maria_lr_geom_mod,  only:   &
    slrproj_tangent
implicit none (type, external)

    character(100) :: seed_str
    integer        :: seed, ldA
    real(WP), allocatable :: A(:)

    call get_command_argument(1, seed_str)
    read(seed_str, *) seed

    call test_sanity()
    call test_quick()
    call test(seed)

contains
    subroutine test_sanity()
        integer                 :: ierr
        real(WP)                :: foo(1)
        procedure(MM),  pointer :: fun

        call slrproj_tangent(-1, -1, -1, foo, 0, foo, 0, fun, fun, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -1 )

        call slrproj_tangent(3, -1, -1, foo, 0, foo, 0, fun, fun, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -2 )

        call slrproj_tangent(3, 3, -1, foo, 0, foo, 0, fun, fun, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -3 )

        call slrproj_tangent(3, 3, 4, foo, 0, foo, 0, fun, fun, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -3 )

        call slrproj_tangent(3, 3, 1, foo, 0, foo, 0, fun, fun, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -5 )

        call slrproj_tangent(3, 3, 1, foo, 3, foo, 0, fun, fun, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -7 )

        call slrproj_tangent(3, 3, 1, foo, 3, foo, 3, fun, fun, foo, 0, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -11 )

        call slrproj_tangent(3, 3, 1, foo, 3, foo, 3, fun, fun, foo, 3, foo, 0, foo, 0, ierr)
        ASSERT( ierr == -13 )

        call slrproj_tangent(3, 3, 1, foo, 3, foo, 3, fun, fun, foo, 3, foo, 3, foo, 0, ierr)
        ASSERT( ierr == -15 )
    end subroutine test_sanity

    subroutine test_quick()
        integer                :: ierr
        real(WP)               :: foo(1)
        procedure(MM), pointer :: fun

        call slrproj_tangent(0, 1, 0, foo, 1, foo, 1, fun, fun, foo, 1, foo, 1, foo, 1, ierr)
        ASSERT( ierr == 0 )

        call slrproj_tangent(1, 0, 0, foo, 1, foo, 1, fun, fun, foo, 1, foo, 1, foo, 1, ierr)
        ASSERT( ierr == 0 )

        call slrproj_tangent(1, 1, 0, foo, 1, foo, 1, fun, fun, foo, 1, foo, 1, foo, 1, ierr)
        ASSERT( ierr == 0 )
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer                :: m, n, r, ldu, ldvt, ldpu, ldpvt, ldutav, lwork, ierr
        real(WP)               :: err, nrm, foo(1)
        type(prng)             :: rng
        procedure(MM), pointer :: fun, B2AB, B2BA
        real(WP), allocatable  :: U(:), VT(:), pU(:), pVT(:), UTAV(:), work(:), tmp(:), pUV(:)
        
        call rng%init(seed, ierr)
        ASSERT( ierr == 0 )

        B2AB => mm_A_left
        B2BA => mm_A_right

    !-- Same column-spaces
        m = 20
        n = 25
        r = 15
        lda = m + 1
        ldu = m + 1
        ldvt = r + 1
        ldpu = m + 1
        ldpvt = r + 1
        ldutav = r + 1

        call srandort(rng, m, r, foo, ldu, foo, -1, ierr)
        lwork = int(foo(1))
        call srandort(rng, r, n, foo, ldvt, foo, -1, ierr)
        lwork = max(lwork, int(foo(1)))

        allocate( A(lda*n), U(ldu*r), VT(ldvt*n), pU(ldpu*r), pVT(ldpvt*n), UTAV(ldutav*r), work(lwork), tmp(n*r) )

        call srandort(rng, m, r, U, ldu, work, lwork, ierr)
        call srandort(rng, r, n, VT, ldvt, work, lwork, ierr)

        call rng%snormal(r*n, tmp, ZERO, ONE, ierr)
        call sgemm('n', 'n', m, n, r, ONE, U, ldu, tmp, r, ZERO, A, lda)
        nrm = sgenrmf(m, n, A, lda, ierr)

        call slrproj_tangent(m, n, r, U, ldu, VT, ldvt, B2AB, B2BA, pU, ldpu, pVT, ldpvt, UTAV, ldutav, ierr)
        call sgemm('n', 'n', m, n, r, -ONE, U, ldu, pVT, ldpvt, ONE, A, lda)
        call sgemm('n', 'n', m, n, r, -ONE, pU, ldpu, VT, ldvt, ONE, A, lda)
        call sgemm('n', 'n', r, n, r, ONE, UTAV, ldutav, VT, ldvt, ZERO, tmp, r)
        call sgemm('n', 'n', m, n, r, -ONE, U, ldu, tmp, r, ONE, A, lda)
        err = sgenrmf(m, n, A, lda, ierr)
        ASSERT( safe_leq(err/nrm, 1e-5_WP) )

        deallocate(A, U, VT, pU, pVT, UTAV, work, tmp)

    !-- Same row-spaces
        m = 25
        n = 20
        r = 15
        lda = m + 1
        ldu = m + 1
        ldvt = r + 1
        ldpu = m + 1
        ldpvt = r + 1
        ldutav = r + 1

        call srandort(rng, m, r, foo, ldu, foo, -1, ierr)
        lwork = int(foo(1))
        call srandort(rng, r, n, foo, ldvt, foo, -1, ierr)
        lwork = max(lwork, int(foo(1)))

        allocate( A(lda*n), U(ldu*r), VT(ldvt*n), pU(ldpu*r), pVT(ldpvt*n), UTAV(ldutav*r), work(lwork), tmp(m*r) )

        call srandort(rng, m, r, U, ldu, work, lwork, ierr)
        call srandort(rng, r, n, VT, ldvt, work, lwork, ierr)

        call rng%snormal(m*r, tmp, ZERO, ONE, ierr)
        call sgemm('n', 'n', m, n, r, ONE, tmp, m, VT, ldvt, ZERO, A, lda)
        nrm = sgenrmf(m, n, A, lda, ierr)

        call slrproj_tangent(m, n, r, U, ldu, VT, ldvt, B2AB, B2BA, pU, ldpu, pVT, ldpvt, UTAV, ldutav, ierr)
        call sgemm('n', 'n', m, n, r, -ONE, U, ldu, pVT, ldpvt, ONE, A, lda)
        call sgemm('n', 'n', m, n, r, -ONE, pU, ldpu, VT, ldvt, ONE, A, lda)
        call sgemm('n', 'n', m, r, r, ONE, U, ldu, UTAV, ldutav, ZERO, tmp, m)
        call sgemm('n', 'n', m, n, r, -ONE, tmp, m, VT, ldvt, ONE, A, lda)
        err = sgenrmf(m, n, A, lda, ierr)
        ASSERT( safe_leq(err/nrm, 1e-5_WP) )

        deallocate(A, U, VT, pU, pVT, UTAV, work, tmp)

    !-- Orthogonal column- and row-spaces
        m = 20
        n = 25
        r = 15
        lda = m + 1
        ldu = m + 1
        ldvt = r + 1
        ldpu = m + 1
        ldpvt = r + 1
        ldutav = r + 1

        call srandort(rng, m, r, foo, ldu, foo, -1, ierr)
        lwork = int(foo(1))
        call srandort(rng, r, n, foo, ldvt, foo, -1, ierr)
        lwork = max(lwork, int(foo(1)))

        allocate( A(lda*n), U(ldu*r), VT(ldvt*n), pU(ldpu*r), pVT(ldpvt*n), UTAV(ldutav*r), work(lwork), tmp(m*n), pUV(max(m,n)**2) )

        call srandort(rng, m, r, U, ldu, work, lwork, ierr)
        call srandort(rng, r, n, VT, ldvt, work, lwork, ierr)

        call rng%snormal(lda*n, A, ZERO, ONE, ierr)

        call slaset('a', m, m, ZERO, ONE, pUV, m)
        call sgemm('n', 't', m, m, r, -ONE, U, ldu, U, ldu, ONE, pUV, m)
        call sgemm('n', 'n', m, n, m, ONE, pUV, m, A, lda, ZERO, tmp, m)
        call slaset('a', n, n, ZERO, ONE, pUV, n)
        call sgemm('t', 'n', n, n, r, -ONE, VT, ldvt, VT, ldvt, ONE, pUV, n)
        call sgemm('n', 'n', m, n, n, ONE, tmp, m, pUV, n, ZERO, A, lda)
        nrm = sgenrmf(m, n, A, lda, ierr)

        call slrproj_tangent(m, n, r, U, ldu, VT, ldvt, B2AB, B2BA, pU, ldpu, pVT, ldpvt, UTAV, ldutav, ierr)
        call sgemm('n', 'n', m, n, r, ONE, U, ldu, pVT, ldpvt, ZERO, A, lda)
        call sgemm('n', 'n', m, n, r, ONE, pU, ldpu, VT, ldvt, ONE, A, lda)
        call sgemm('n', 'n', m, r, r, ONE, U, ldu, UTAV, ldutav, ZERO, tmp, m)
        call sgemm('n', 'n', m, n, r, ONE, tmp, m, VT, ldvt, ONE, A, lda)
        err = sgenrmf(m, n, A, lda, ierr)
        ASSERT( safe_leq(err/nrm, 1e-5_WP) )
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

        call sgemm('n', transB, m, n, k, alpha, A, ldA, B, ldB, beta, C, ldc)
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

        call sgemm(transB, 'n', m, n, k, alpha, B, ldB, A, ldA, beta, C, ldc)
        info = 0
    end subroutine mm_A_right

end program slrproj_tangent_test
