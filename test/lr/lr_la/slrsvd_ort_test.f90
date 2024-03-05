program slrsvd_ort_test
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
use maria_la_utils_mod,     only: &
    all_const,                    &
    all_close
use maria_la_core_mod,      only: &
    sgemm,                        &
    sgesdd_q,                     &
    sorcangles
use maria_lr_la_mod,        only: &
    slrort,                       &
    slrsvd_ort,                   &
    lrort_rank
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
        integer    :: ifoo(1), info
        real(WP)   :: foo(1)

        call slrsvd_ort('x', 'x', -1, -1, -1, foo, 0, foo, 0, foo, foo, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -1)

        call slrsvd_ort('s', 'x', -1, -1, -1, foo, 0, foo, 0, foo, foo, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -2)

        call slrsvd_ort('s', 'l', -1, -1, -1, foo, 0, foo, 0, foo, foo, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -3)

        call slrsvd_ort('s', 'l', 0, -1, -1, foo, 0, foo, 0, foo, foo, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -4)

        call slrsvd_ort('s', 'l', 0, 0, -1, foo, 0, foo, 0, foo, foo, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -5)

        call slrsvd_ort('s', 'l', 0, 1, 1, foo, 0, foo, 0, foo, foo, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -5)

        call slrsvd_ort('s', 'r', 1, 0, 1, foo, 0, foo, 0, foo, foo, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -5)

        call slrsvd_ort('s', 'l', 0, 0, 0, foo, 0, foo, 0, foo, foo, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -7)

        call slrsvd_ort('s', 'l', 0, 0, 0, foo, 1, foo, 0, foo, foo, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -9)

        call slrsvd_ort('s', 'l', 0, 0, 0, foo, 1, foo, 1, foo, foo, foo, 0, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -13)

        call slrsvd_ort('s', 'l', 2, 0, 0, foo, 2, foo, 1, foo, foo, foo, 1, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -13)

        call slrsvd_ort('s', 'l', 0, 0, 0, foo, 1, foo, 1, foo, foo, foo, 1, foo, 0, foo, 0, ifoo, 0, info)
        ASSERT(info == -15)

        call slrsvd_ort('s', 'l', 5, 3, 2, foo, 5, foo, 2, foo, foo, foo, 5, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == -15)

        call slrsvd_ort('s', 'l', 5, 2, 3, foo, 5, foo, 3, foo, foo, foo, 5, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == -15)

        call slrsvd_ort('s', 'r', 2, 4, 3, foo, 2, foo, 3, foo, foo, foo, 3, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == -15)

        call slrsvd_ort('s', 'r', 3, 4, 2, foo, 3, foo, 2, foo, foo, foo, 3, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == -15)
    end subroutine test_sanity

    subroutine test_quick()
        integer    :: ifoo(1), info
        real(WP)   :: foo(1)

        call slrsvd_ort('s', 'l', 0, 1, 0, foo, 1, foo, 1, foo, foo, foo, 1, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)

        call slrsvd_ort('s', 'l', 1, 0, 0, foo, 1, foo, 1, foo, foo, foo, 1, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)

        call slrsvd_ort('s', 'l', 1, 1, 0, foo, 1, foo, 1, foo, foo, foo, 1, foo, 1, foo, 0, ifoo, 0, info)
        ASSERT(info == 0)
    end subroutine test_quick

    subroutine test(seed)
        integer, intent(in) :: seed

        character(1)          :: job, side
        integer               :: info, m, n, r, ldu, ldvt, ldu1, ldv1t, ldu2, ldv2t, lwork, liwork, ifoo(1)
        real(WP)              :: foo(1)
        type(prng)            :: rng     
        integer,  allocatable :: iwork(:)  
        real(WP), allocatable :: U(:), VT(:), tau(:), A(:), U1(:), S1(:), V1T(:), &
                                 U2(:), S2(:), V2T(:), work(:)
 
        call rng%init(seed, info)
        ASSERT(info == 0)

    !-- job = 's', side = 'l', r < n ------------------------------------------------------
        job = 's'
        side = 'l'
        m = 20
        n = 10
        r = 8
        ldU = m + 1
        ldVT = r + 1
        ldU1 = m + 1
        ldV1T = min(r,n) + 1
        ldU2 = m + 1
        ldV2T = min(m,n) + 1

        call slrort(side, m, n, r, foo, ldU, foo, ldVT, foo, foo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))

        call slrsvd_ort(job, side, m, n, r, foo, ldu, foo, ldvt, foo, foo, foo, ldu1, foo, ldv1t, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = ifoo(1)

        call sgesdd_q(job, m, n, foo, m, foo, foo, ldu2, foo, ldv2t, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
 
        call sorcangles('c', m, min(r,n), foo, ldu1, foo, ldu2, foo, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        call sorcangles('r', min(r,n), n, foo, ldv1t, foo, ldv2t, foo, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
       
        allocate(U(ldU*r), VT(ldVT*n), tau(r), A(m*n), U1(ldu1*min(r,n)), S1(min(r,n)), V1T(ldV1T*n), &
            U2(ldu2*min(m,n)), S2(min(m,n)), V2T(ldv2t*n), work(lwork), iwork(liwork))
        call rng%snormal(ldU*r, U, ZERO, ONE, info)
        call rng%snormal(ldVT*n, VT, ZERO, ONE, info)
        call sgemm('n', 'n', m, n, r, ONE, U, ldU, VT, ldVT, ZERO, A, m)
        call slrort(side, m, n, r, U, ldU, VT, ldVT, tau, work, lwork, info)
        call slrsvd_ort(job, side, m, n, r, U, ldu, VT, ldvt, tau, &
            S1, U1, ldu1, V1T, ldv1t, work, lwork, iwork, liwork, info)
        call sgesdd_q(job, m, n, A, m, S2, U2, ldu2, V2T, ldv2t, work, lwork, iwork, liwork, info)
       
        ASSERT(all_close(min(r,n), S1, 1, S2, 1, info, rtol=1e-5_WP))
        ASSERT(all_const(min(m,n)-min(r,n), S2(min(r,n)+1:), 1, ZERO, info, atol=1e-5_WP))

        call sorcangles('c', m, min(r,n), U1, ldu1, U2, ldu2, S1, work, lwork, iwork, liwork, info)
        ASSERT(all_const(min(r,n), S1, 1, ONE, info, atol=1e-5_WP))
        call sorcangles('r', min(r,n), n, V1T, ldv1t, V2T, ldv2t, S1, work, lwork, iwork, liwork, info)
        ASSERT(all_const(min(r,n), S1, 1, ONE, info, atol=1e-5_WP))

        deallocate(U, VT, tau, A, U1, S1, V1T, U2, S2, V2T, work, iwork)

    !-- job = 's', side = 'l', r > n ------------------------------------------------------
        job = 's'
        side = 'l'
        m = 20
        n = 10
        r = 12
        ldU = m + 1
        ldVT = r + 1
        ldU1 = m + 1
        ldV1T = min(r,n) + 1
        ldU2 = m + 1
        ldV2T = min(m,n) + 1

        call slrort(side, m, n, r, foo, ldU, foo, ldVT, foo, foo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))

        call slrsvd_ort(job, side, m, n, r, foo, ldu, foo, ldvt, foo, foo, foo, ldu1, foo, ldv1t, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = ifoo(1)

        call sgesdd_q(job, m, n, foo, m, foo, foo, ldu2, foo, ldv2t, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        call sorcangles('c', m, min(r,n), foo, ldu1, foo, ldu2, foo, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        call sorcangles('r', min(r,n), n, foo, ldv1t, foo, ldv2t, foo, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        allocate(U(ldU*r), VT(ldVT*n), tau(r), A(m*n), U1(ldu1*min(r,n)), S1(min(r,n)), V1T(ldV1T*n), &
            U2(ldu2*min(m,n)), S2(min(m,n)), V2T(ldv2t*n), work(lwork), iwork(liwork))
        call rng%snormal(ldU*r, U, ZERO, ONE, info)
        call rng%snormal(ldVT*n, VT, ZERO, ONE, info)
        call sgemm('n', 'n', m, n, r, ONE, U, ldU, VT, ldVT, ZERO, A, m)
        call slrort(side, m, n, r, U, ldU, VT, ldVT, tau, work, lwork, info)
        call slrsvd_ort(job, side, m, n, r, U, ldu, VT, ldvt, tau, &
            S1, U1, ldu1, V1T, ldv1t, work, lwork, iwork, liwork, info)
        call sgesdd_q(job, m, n, A, m, S2, U2, ldu2, V2T, ldv2t, work, lwork, iwork, liwork, info)

        ASSERT(all_close(min(r,n), S1, 1, S2, 1, info, rtol=1e-5_WP))
        ASSERT(all_const(min(m,n)-min(r,n), S2(min(r,n)+1:), 1, ZERO, info, atol=1e-5_WP))

        call sorcangles('c', m, min(r,n), U1, ldu1, U2, ldu2, S1, work, lwork, iwork, liwork, info)
        ASSERT(all_const(min(r,n), S1, 1, ONE, info, atol=1e-5_WP))
        call sorcangles('r', min(r,n), n, V1T, ldv1t, V2T, ldv2t, S1, work, lwork, iwork, liwork, info)
        ASSERT(all_const(min(r,n), S1, 1, ONE, info, atol=1e-5_WP))

        deallocate(U, VT, tau, A, U1, S1, V1T, U2, S2, V2T, work, iwork)

    !-- job = 's', side = 'r', r < m ------------------------------------------------------
        job = 's'
        side = 'r'
        m = 3
        n = 4
        r = 2
        ldU = m + 1
        ldVT = r + 1
        ldU1 = m + 1
        ldV1T = min(m,r) + 1
        ldU2 = m + 1
        ldV2T = min(m,n) + 1

        call slrort(side, m, n, r, foo, ldU, foo, ldVT, foo, foo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))

        call slrsvd_ort(job, side, m, n, r, foo, ldu, foo, ldvt, foo, foo, foo, ldu1, foo, ldv1t, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = ifoo(1)

        call sgesdd_q(job, m, n, foo, m, foo, foo, ldu2, foo, ldv2t, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
 
        call sorcangles('c', m, min(m,r), foo, ldu1, foo, ldu2, foo, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        call sorcangles('r', min(m,r), n, foo, ldv1t, foo, ldv2t, foo, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
       
        allocate(U(ldU*r), VT(ldVT*n), tau(r), A(m*n), U1(ldu1*min(m,r)), S1(min(m,r)), V1T(ldV1T*n), &
            U2(ldu2*min(m,n)), S2(min(m,n)), V2T(ldv2t*n), work(lwork), iwork(liwork))
        call rng%snormal(ldU*r, U, ZERO, ONE, info)
        call rng%snormal(ldVT*n, VT, ZERO, ONE, info)
        call sgemm('n', 'n', m, n, r, ONE, U, ldU, VT, ldVT, ZERO, A, m)
        call slrort(side, m, n, r, U, ldU, VT, ldVT, tau, work, lwork, info)
        call slrsvd_ort(job, side, m, n, r, U, ldu, VT, ldvt, tau, &
            S1, U1, ldu1, V1T, ldv1t, work, lwork, iwork, liwork, info)
        call sgesdd_q(job, m, n, A, m, S2, U2, ldu2, V2T, ldv2t, work, lwork, iwork, liwork, info)
       
        ASSERT(all_close(min(m,r), S1, 1, S2, 1, info, rtol=1e-5_WP))
        ASSERT(all_const(min(m,n)-min(m,r), S2(min(m,r)+1:), 1, ZERO, info, atol=1e-5_WP))

        call sorcangles('c', m, min(m,r), U1, ldu1, U2, ldu2, S1, work, lwork, iwork, liwork, info)
        ASSERT(all_const(min(m,r), S1, 1, ONE, info, atol=1e-5_WP))

        print '( 2(1pe12.5) )', tau
        print '( 4(1pe12.5) )', VT(1), VT(1 + ldvt), VT(1 + 2*ldvt), VT(1 + 3*ldvt)
        print '( 4(1pe12.5) )', VT(2), VT(2 + ldvt), VT(2 + 2*ldvt), VT(2 + 3*ldvt)
        print '( 4(1pe12.5) )', V1T(1), V1T(1 + ldv1t), V1T(1 + 2*ldv1t), V1T(1 + 3*ldv1t)
        print '( 4(1pe12.5) )', V1T(2), V1T(2 + ldv1t), V1T(2 + 2*ldv1t), V1T(2 + 3*ldv1t)

        call sorcangles('r', min(m,r), n, V1T, ldv1t, V2T, ldv2t, S1, work, lwork, iwork, liwork, info)
        ASSERT(all_const(min(m,r), S1, 1, ONE, info, atol=1e-5_WP))

        deallocate(U, VT, tau, A, U1, S1, V1T, U2, S2, V2T, work, iwork)

    !-- job = 's', side = 'r', r > m ------------------------------------------------------
        job = 's'
        side = 'r'
        m = 10
        n = 20
        r = 15
        ldU = m + 1
        ldVT = r + 1
        ldU1 = m + 1
        ldV1T = min(m,r) + 1
        ldU2 = m + 1
        ldV2T = min(m,n) + 1

        call slrort(side, m, n, r, foo, ldU, foo, ldVT, foo, foo, -1, info)
        ASSERT(info == 0)
        lwork = int(foo(1))

        call slrsvd_ort(job, side, m, n, r, foo, ldu, foo, ldvt, foo, foo, foo, ldu1, foo, ldv1t, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = ifoo(1)

        call sgesdd_q(job, m, n, foo, m, foo, foo, ldu2, foo, ldv2t, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
 
        call sorcangles('c', m, min(m,r), foo, ldu1, foo, ldu2, foo, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))

        call sorcangles('r', min(m,r), n, foo, ldv1t, foo, ldv2t, foo, foo, -1, ifoo, -1, info)
        ASSERT(info == 0)
        lwork = max(lwork, int(foo(1)))
        liwork = max(liwork, ifoo(1))
       
        allocate(U(ldU*r), VT(ldVT*n), tau(r), A(m*n), U1(ldu1*min(m,r)), S1(min(m,r)), V1T(ldV1T*n), &
            U2(ldu2*min(m,n)), S2(min(m,n)), V2T(ldv2t*n), work(lwork), iwork(liwork))
        call rng%snormal(ldU*r, U, ZERO, ONE, info)
        call rng%snormal(ldVT*n, VT, ZERO, ONE, info)
        call sgemm('n', 'n', m, n, r, ONE, U, ldU, VT, ldVT, ZERO, A, m)
        call slrort(side, m, n, r, U, ldU, VT, ldVT, tau, work, lwork, info)
        call slrsvd_ort(job, side, m, n, r, U, ldu, VT, ldvt, tau, &
            S1, U1, ldu1, V1T, ldv1t, work, lwork, iwork, liwork, info)
        call sgesdd_q(job, m, n, A, m, S2, U2, ldu2, V2T, ldv2t, work, lwork, iwork, liwork, info)
       
        ASSERT(all_close(min(m,r), S1, 1, S2, 1, info, rtol=1e-5_WP))
        ASSERT(all_const(min(m,n)-min(m,r), S2(min(m,r)+1:), 1, ZERO, info, atol=1e-5_WP))

        call sorcangles('c', m, min(m,r), U1, ldu1, U2, ldu2, S1, work, lwork, iwork, liwork, info)
        ASSERT(all_const(min(m,r), S1, 1, ONE, info, atol=1e-5_WP))
        call sorcangles('r', min(m,r), n, V1T, ldv1t, V2T, ldv2t, S1, work, lwork, iwork, liwork, info)
        ASSERT(all_const(min(m,r), S1, 1, ONE, info, atol=1e-5_WP))

        deallocate(U, VT, tau, A, U1, S1, V1T, U2, S2, V2T, work, iwork)
    end subroutine test
end program slrsvd_ort_test
