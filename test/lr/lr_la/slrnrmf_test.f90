program slrnrmf_test
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
use maria_utils_mod,        only: &
    are_close
use maria_la_core_mod,      only: &
    sgemm,                        &
    sgenrmf
use maria_lr_la_mod,        only: &
    slrnrmf
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
        real(WP)   :: nrm, foo(1)

        nrm = slrnrmf(-1, -1, -1, foo, 0, foo, 0, foo, 0, info)
        ASSERT(info == -1)

        nrm = slrnrmf(1, -1, -1, foo, 0, foo, 0, foo, 0, info)
        ASSERT(info == -2)

        nrm = slrnrmf(1, 1, -1, foo, 0, foo, 0, foo, 0, info)
        ASSERT(info == -3)

        nrm = slrnrmf(1, 1, 1, foo, 0, foo, 0, foo, 0, info)
        ASSERT(info == -5)

        nrm = slrnrmf(1, 1, 1, foo, 1, foo, 0, foo, 0, info)
        ASSERT(info == -7)

        nrm = slrnrmf(1, 1, 1, foo, 1, foo, 1, foo, 0, info)
        ASSERT(info == -9)
    end subroutine test_sanity

    subroutine test_quick()
        integer    :: info
        real(WP)   :: nrm, foo(1)

        nrm = slrnrmf(0, 1, 1, foo, 1, foo, 1, foo, 0, info)
        ASSERT(info == 0)

        nrm = slrnrmf(1, 0, 1, foo, 1, foo, 1, foo, 0, info)
        ASSERT(info == 0)

        nrm = slrnrmf(1, 1, 0, foo, 1, foo, 1, foo, 0, info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        integer               :: info, m, n, r, ldU, ldVT, lwork
        real(WP)              :: foo(1), nrm1, nrm2
        type(prng)            :: rng       
        real(WP), allocatable :: U(:), VT(:), work(:), A(:)
 
        call rng%init(seed, info)
        ASSERT(info == 0)

        m = 20
        n = 10
        r = 15
        ldU = m + 1
        ldVT = r + 1
        nrm1 = slrnrmf(m, n, r, foo, ldU, foo, ldVT, foo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))
        allocate(U(ldU*r), VT(ldVT*n), A(m*n), work(lwork))
        call rng%snormal(ldU*r, U, ZERO, ONE, info)
        call rng%snormal(ldVT*n, VT, ZERO, ONE, info)
        nrm1 = slrnrmf(m, n, r, U, ldU, VT, ldVT, work, lwork, info)
        ASSERT(info == 0)
        call sgemm('n', 'n', m, n, r, ONE, U, ldU, VT, ldVT, ZERO, A, m)
        nrm2 = sgenrmf(m, n, A, m, info)
        ASSERT(are_close(nrm1, nrm2, info, rtol=1e-5_WP))

        call rng%deinit(info)
        ASSERT(info == 0)
    end subroutine test
end program slrnrmf_test
