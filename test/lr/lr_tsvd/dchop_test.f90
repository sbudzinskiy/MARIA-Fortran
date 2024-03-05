program dchop_test
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
    safe_eq,                      &
    safe_leq,                     &
    safe_less
use maria_utils_mod,        only: &
    are_close
use maria_la_core_mod,      only: &
    sgemm,                        &
    sgenrmf
use maria_lr_tsvd_mod,      only: &
    dchop
implicit none (type, external)

    character(100) :: seed_str
    integer        :: seed

    call get_command_argument(1, seed_str)
    read(seed_str, *) seed

    call test_sanity()
    call test_quick()
    call test()

contains
    subroutine test_sanity()
        integer    :: info, rk
        real(WP)   :: foo(1)

        foo(1) = ZERO

        rk = dchop(-1, foo, info)
        ASSERT(info == -1)

        rk = dchop(1, foo, info, maxr=-1)
        ASSERT(info == -4)

        rk = dchop(1, foo, info, rtolf=-ONE)
        ASSERT(info == -5)

        rk = dchop(1, foo, info, atolf=-ONE)
        ASSERT(info == -6)

        rk = dchop(1, foo, info, rtol2=-ONE)
        ASSERT(info == -7)

        rk = dchop(1, foo, info, atol2=-ONE)
        ASSERT(info == -8)
    end subroutine test_sanity

    subroutine test_quick()
        integer    :: info, rk
        real(WP)   :: foo(1), rerrf, aerrf, rerr2, aerr2

        foo(1) = ZERO

        rk = dchop(0, foo, info, rerrf=rerrf, aerrf=aerrf, rerr2=rerr2, aerr2=aerr2)
        ASSERT(info == 0)
        ASSERT(rk == 0)
        ASSERT(safe_eq(rerrf, ZERO))
        ASSERT(safe_eq(aerrf, ZERO))
        ASSERT(safe_eq(rerr2, ZERO))
        ASSERT(safe_eq(aerr2, ZERO))

        rk = dchop(1, foo, info, rerrf=rerrf, aerrf=aerrf, rerr2=rerr2, aerr2=aerr2)
        ASSERT(info == 0)
        ASSERT(rk == 0)
        ASSERT(safe_eq(rerrf, ZERO))
        ASSERT(safe_eq(aerrf, ZERO))
        ASSERT(safe_eq(rerr2, ZERO))
        ASSERT(safe_eq(aerr2, ZERO))
    end subroutine test_quick

    subroutine test()
        integer               :: rk, n, info, i, maxr
        real(WP)              :: rtolf, atolf, rtol2, atol2, rerrf, aerrf, rerr2, aerr2, summ
        real(WP), allocatable :: S(:)
 
        n = 4
        allocate(S(n))
        do i = 1, n
            S(i) = 10**(n - i)
        end do
        
    !-- maxr -------------------------------------------------------------------
        maxr = 2
        rk = dchop(n, S, info, maxr=maxr, rerrf=rerrf, aerrf=aerrf, rerr2=rerr2, aerr2=aerr2)
        ASSERT(info == 0)
        ASSERT(rk == maxr)
        ASSERT(safe_eq(aerr2, S(rk+1)))
        ASSERT(safe_eq(rerr2, S(rk+1)/S(1)))
        summ = sum(S(rk+1:n)**2)
        ASSERT(safe_eq(aerrf, sqrt(summ)))
        ASSERT(safe_eq(rerrf, sqrt(summ/sum(S**2))))

    !-- atol2 ------------------------------------------------------------------
        atol2 = 2 * ONE
        rk = dchop(n, S, info, atol2=atol2, aerr2=aerr2)
        ASSERT(info == 0)
        ASSERT(rk == 3)
        ASSERT(safe_leq(aerr2, atol2))
        ASSERT(safe_less(atol2, S(rk)))

    !-- rtol2 ------------------------------------------------------------------
        rtol2 = 0.1_WP
        rk = dchop(n, S, info, rtol2=rtol2, rerr2=rerr2)
        ASSERT(info == 0)
        ASSERT(rk == 1)
        ASSERT(safe_leq(rerr2, rtol2))
        ASSERT(safe_less(rtol2, S(rk)/S(1)))

    !-- atolf ------------------------------------------------------------------
        atolf = 11 * ONE
        rk = dchop(n, S, info, atolf=atolf, aerrf=aerrf)
        ASSERT(info == 0)
        ASSERT(rk == 2)
        ASSERT(safe_leq(aerrf, atolf))
        ASSERT(safe_less(atolf, sqrt(sum(S(rk:n)**2))))

    !-- rtolf ------------------------------------------------------------------
        rtolf = 1e-4_WP
        rk = dchop(n, S, info, rtolf=rtolf, rerrf=rerrf)
        ASSERT(info == 0)
        ASSERT(rk == 4)
        ASSERT(safe_leq(rerrf, rtolf))
        ASSERT(safe_less(rtolf, sqrt(sum(S(rk:n)**2 / sum(S(1:n)**2)))))

    !-- atolf, atol2 -----------------------------------------------------------
        rk = dchop(n, S, info, atolf=atolf, atol2=atol2, aerrf=aerrf, aerr2=aerr2)
        ASSERT(info == 0)
        ASSERT(rk == 3)
        ASSERT(safe_leq(aerrf, atolf))
        ASSERT(safe_leq(aerr2, atol2))

    !-- maxr, rtol2 ------------------------------------------------------------
        rk = dchop(n, S, info, maxr=maxr, rtol2=rtol2, rerr2=rerr2)
        ASSERT(info == 0)
        ASSERT(rk == 1)
        ASSERT(safe_leq(rerr2, rtol2))

    end subroutine test
end program dchop_test
