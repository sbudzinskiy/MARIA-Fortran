set(MARIA_CTEST_LABEL_UTILS maria_utils)

add_new_test(loptional_val
    LINKS maria_assert maria_utils
    LABELS ${MARIA_CTEST_LABEL_UTILS} ${MARIA_CTEST_LABEL_L})

add_new_test(ioptional_val
    LINKS maria_assert maria_utils
    LABELS ${MARIA_CTEST_LABEL_UTILS} ${MARIA_CTEST_LABEL_I})

add_new_test(soptional_val
    LINKS maria_assert maria_kinds maria_comparison maria_utils
    LABELS ${MARIA_CTEST_LABEL_UTILS} ${MARIA_CTEST_LABEL_S})

add_new_test(doptional_val
    LINKS maria_assert maria_kinds maria_comparison maria_utils
    LABELS ${MARIA_CTEST_LABEL_UTILS} ${MARIA_CTEST_LABEL_D})

add_new_test(sare_close
    LINKS maria_assert maria_kinds maria_constants maria_utils
    LABELS ${MARIA_CTEST_LABEL_UTILS} ${MARIA_CTEST_LABEL_S})

add_new_test(dare_close
    LINKS maria_assert maria_kinds maria_constants maria_utils
    LABELS ${MARIA_CTEST_LABEL_UTILS} ${MARIA_CTEST_LABEL_D})

add_new_test(iswap_pair
    LINKS maria_assert maria_utils
    LABELS ${MARIA_CTEST_LABEL_UTILS} ${MARIA_CTEST_LABEL_I})

add_new_test(sswap_pair
    LINKS maria_assert maria_kinds maria_comparison maria_utils
    LABELS ${MARIA_CTEST_LABEL_UTILS} ${MARIA_CTEST_LABEL_S})

add_new_test(dswap_pair
    LINKS maria_assert maria_kinds maria_comparison maria_utils
    LABELS ${MARIA_CTEST_LABEL_UTILS} ${MARIA_CTEST_LABEL_D})

add_new_test(slinspace
    LINKS maria_assert maria_kinds maria_constants maria_utils
    LABELS ${MARIA_CTEST_LABEL_UTILS} ${MARIA_CTEST_LABEL_S})

add_new_test(dlinspace
    LINKS maria_assert maria_kinds maria_constants maria_utils
    LABELS ${MARIA_CTEST_LABEL_UTILS} ${MARIA_CTEST_LABEL_D})

add_new_test(iarange
    LINKS maria_assert maria_utils
    LABELS ${MARIA_CTEST_LABEL_UTILS} ${MARIA_CTEST_LABEL_I})

add_new_test(sarange
    LINKS maria_assert maria_kinds maria_constants maria_utils
    LABELS ${MARIA_CTEST_LABEL_UTILS} ${MARIA_CTEST_LABEL_S})

add_new_test(darange
    LINKS maria_assert maria_kinds maria_constants maria_utils
    LABELS ${MARIA_CTEST_LABEL_UTILS} ${MARIA_CTEST_LABEL_D})

add_new_test(imedian_of_3
    LINKS maria_assert maria_utils
    LABELS ${MARIA_CTEST_LABEL_UTILS} ${MARIA_CTEST_LABEL_I})

add_new_test(isort   
    LINKS maria_assert maria_utils ${MARIA_PRNG_IMPLEMENTATION}
    LABELS ${MARIA_CTEST_LABEL_UTILS} ${MARIA_CTEST_LABEL_I})

add_new_test(iselect_no_replacement
    LINKS maria_assert maria_utils ${MARIA_PRNG_IMPLEMENTATION}
    LABELS ${MARIA_CTEST_LABEL_UTILS} ${MARIA_CTEST_LABEL_I})

add_new_test(ipermute
    LINKS maria_assert maria_utils ${MARIA_PRNG_IMPLEMENTATION}
    LABELS ${MARIA_CTEST_LABEL_UTILS} ${MARIA_CTEST_LABEL_I})

add_new_test(srealloc
    LINKS maria_assert maria_kinds maria_constants maria_utils
    LABELS ${MARIA_CTEST_LABEL_UTILS} ${MARIA_CTEST_LABEL_S})

add_new_test(drealloc
    LINKS maria_assert maria_kinds maria_constants maria_utils
    LABELS ${MARIA_CTEST_LABEL_UTILS} ${MARIA_CTEST_LABEL_D})

