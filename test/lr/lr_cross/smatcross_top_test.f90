program smatcross_top_test
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
    sgeall_const
use maria_la_core_mod,      only: &
    slacpy,                       &
    slaset,                       &
    sgemm,                        &
    sgenrmc
use maria_lr_cross_mod,     only: &
    smatcross_top
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
        real(WP)               :: foo(1)

        call smatcross_top(-1, -1, -1, foo, 0, -1, foo, 0, -1, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -1)

        call smatcross_top(1, -1, -1, foo, 0, -1, foo, 0, -1, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -2)

        call smatcross_top(1, 1, -1, foo, 0, -1, foo, 0, -1, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -3)

        call smatcross_top(1, 1, 2, foo, 0, -1, foo, 0, -1, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -3)

        call smatcross_top(1, 1, 1, foo, 0, -1, foo, 0, -1, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -5)

        call smatcross_top(1, 1, 1, foo, 1, -1, foo, 0, -1, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -6)

        call smatcross_top(1, 1, 1, foo, 1, 2, foo, 0, -1, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -6)

        call smatcross_top(1, 1, 1, foo, 1, 1, foo, 0, -1, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -8)

        call smatcross_top(1, 1, 1, foo, 1, 1, foo, 1, -1, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -9)

        call smatcross_top(1, 1, 1, foo, 1, 1, foo, 1, 2, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -9)

        call smatcross_top(1, 1, 1, foo, 1, 1, foo, 1, 1, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -11)

        call smatcross_top(1, 1, 1, foo, 1, 1, foo, 1, 1, foo, 1, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -13)
    end subroutine test_sanity

    subroutine test_quick()
        integer                :: info, ifoo(1)
        real(WP)               :: foo(1)

        call smatcross_top(0, 1, 1, foo, 1, 0, foo, 1, 0, foo, 1, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)

        call smatcross_top(1, 0, 0, foo, 1, 1, foo, 1, 0, foo, 1, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)

        call smatcross_top(1, 1, 0, foo, 1, 1, foo, 1, 0, foo, 1, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)

        call smatcross_top(1, 1, 1, foo, 1, 0, foo, 1, 0, foo, 1, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)

        call smatcross_top(1, 1, 1, foo, 1, 1, foo, 1, 0, foo, 1, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer                :: m, n, lda, nc, ldc, nr, ldr, r, ldu, ldvt, lwork, liwork, info, ifoo(1)
        real(WP)               :: nrmc, foo(1)
        type(prng)             :: rng
        integer,  allocatable  :: iwork(:)
        real(WP), allocatable  :: A(:), cols(:), rows(:), U(:), VT(:), work(:)
        
        call rng%init(seed, info)
        ASSERT(info == 0)

    !-- nc = nr = r
        m = 15
        n = 20
        ldA = m + 1
        nc = 4
        ldc = m + 1
        nr = nc
        ldr = nr + 1
        r = nc
        ldu = m + 1
        ldvt = r + 1

        call smatcross_top(m, n, nc, foo, ldc, nr, foo, ldr, r, foo, ldu, foo, ldvt, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)

        allocate(A(lda*n), cols(ldc*nc), rows(ldr*n), U(ldu*r), VT(ldvt*n), work(lwork), iwork(liwork))
        
        call rng%snormal(ldu*r, U, ZERO, ONE, info)
        call rng%snormal(ldvt*n, VT, ZERO, ONE, info)

        call sgemm('n', 'n', m, n, r, ONE, U, ldu, VT, ldvt, ZERO, A, lda)
        nrmc = sgenrmc(m, n, A, lda, info)
        call slacpy('a', m, nc, A, lda, cols, ldc)
        call slacpy('a', nr, n, A, lda, rows, ldr)
        call smatcross_top(m, n, nc, cols, ldc, nr, rows, ldr, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
        ASSERT(info == 0)
        call sgemm('n', 'n', m, n, r, -ONE, U, ldu, VT, ldvt, ONE, A, lda)
        ASSERT(sgeall_const(m, n, A, lda, ZERO, ZERO, info, atol=1e-5_WP * nrmc))

        deallocate(A, cols, rows,U, VT, work, iwork)

    !-- r < nc < nr
        m = 15
        n = 20
        ldA = m + 1
        nc = 4
        ldc = m + 1
        nr = nc + 1
        ldr = nr + 1
        r = nc - 1
        ldu = m + 1
        ldvt = r + 1

        call smatcross_top(m, n, nc, foo, ldc, nr, foo, ldr, r, foo, ldu, foo, ldvt, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)

        allocate(A(lda*n), cols(ldc*nc), rows(ldr*n), U(ldu*r), VT(ldvt*n), work(lwork), iwork(liwork))
        
        call rng%snormal(ldu*r, U, ZERO, ONE, info)
        call rng%snormal(ldvt*n, VT, ZERO, ONE, info)

        call sgemm('n', 'n', m, n, r, ONE, U, ldu, VT, ldvt, ZERO, A, lda)
        nrmc = sgenrmc(m, n, A, lda, info)
        call slacpy('a', m, nc, A, lda, cols, ldc)
        call slacpy('a', nr, n, A, lda, rows, ldr)
        call smatcross_top(m, n, nc, cols, ldc, nr, rows, ldr, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
        ASSERT(info == 0)
        call sgemm('n', 'n', m, n, r, -ONE, U, ldu, VT, ldvt, ONE, A, lda)
        ASSERT(sgeall_const(m, n, A, lda, ZERO, ZERO, info, atol=1e-5_WP * nrmc))

        deallocate(A, cols, rows,U, VT, work, iwork)

    !-- nc = nr = r, runtimer error
        m = 15
        n = 20
        ldA = m + 1
        nc = 4
        ldc = m + 1
        nr = nc
        ldr = nr + 1
        r = nc
        ldu = m + 1
        ldvt = r + 1

        call smatcross_top(m, n, nc, foo, ldc, nr, foo, ldr, r, foo, ldu, foo, ldvt, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)

        allocate(A(lda*n), cols(ldc*nc), rows(ldr*n), U(ldu*r), VT(ldvt*n), work(lwork), iwork(liwork))
        
        call rng%snormal(ldc*nc, cols, ZERO, ONE, info)
        call rng%snormal(ldr*n, rows, ZERO, ONE, info)
        call slaset('a', nr, nc, ZERO, ZERO, cols, ldc)

        call smatcross_top(m, n, nc, cols, ldc, nr, rows, ldr, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
        ASSERT(info > 0)

        deallocate(A, cols, rows,U, VT, work, iwork)

    !-- r < nc < nr, runtimer error
        m = 15
        n = 20
        ldA = m + 1
        nc = 4
        ldc = m + 1
        nr = nc + 1
        ldr = nr + 1
        r = nc - 1
        ldu = m + 1
        ldvt = r + 1

        call smatcross_top(m, n, nc, foo, ldc, nr, foo, ldr, r, foo, ldu, foo, ldvt, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        liwork = ifoo(1)

        allocate(A(lda*n), cols(ldc*nc), rows(ldr*n), U(ldu*r), VT(ldvt*n), work(lwork), iwork(liwork))
        
        call rng%snormal(ldc*nc, cols, ZERO, ONE, info)
        call rng%snormal(ldr*n, rows, ZERO, ONE, info)
        call slaset('a', nr, nc, ZERO, ZERO, cols, ldc)

        call smatcross_top(m, n, nc, cols, ldc, nr, rows, ldr, r, U, ldu, VT, ldvt, work, lwork, iwork, liwork, info)
        ASSERT(info > 0)

        deallocate(A, cols, rows,U, VT, work, iwork)
    end subroutine test

end program smatcross_top_test
