set(MARIA_CTEST_LABEL_LA_UTILS maria_la_utils)

add_new_test(sall_close
    LINKS maria_assert maria_kinds maria_constants maria_utils maria_la_utils
    LABELS ${MARIA_CTEST_LABEL_LA_UTILS} ${MARIA_CTEST_LABEL_S})

add_new_test(dall_close
    LINKS maria_assert maria_kinds maria_constants maria_utils maria_la_utils
    LABELS ${MARIA_CTEST_LABEL_LA_UTILS} ${MARIA_CTEST_LABEL_D})

add_new_test(sall_const
    LINKS maria_assert maria_kinds maria_constants maria_utils maria_la_utils
    LABELS ${MARIA_CTEST_LABEL_LA_UTILS} ${MARIA_CTEST_LABEL_S})

add_new_test(dall_const
    LINKS maria_assert maria_kinds maria_constants maria_utils maria_la_utils
    LABELS ${MARIA_CTEST_LABEL_LA_UTILS} ${MARIA_CTEST_LABEL_D})

add_new_test(sgeall_close
    LINKS maria_assert maria_kinds maria_constants maria_utils maria_la_core maria_la_utils
    LABELS ${MARIA_CTEST_LABEL_LA_UTILS} ${MARIA_CTEST_LABEL_S})

add_new_test(dgeall_close
    LINKS maria_assert maria_kinds maria_constants maria_utils maria_la_core maria_la_utils
    LABELS ${MARIA_CTEST_LABEL_LA_UTILS} ${MARIA_CTEST_LABEL_D})

add_new_test(sgeall_const
    LINKS maria_assert maria_kinds maria_constants maria_utils maria_la_core maria_la_utils
    LABELS ${MARIA_CTEST_LABEL_LA_UTILS} ${MARIA_CTEST_LABEL_S})

add_new_test(dgeall_const
    LINKS maria_assert maria_kinds maria_constants maria_utils maria_la_core maria_la_utils
    LABELS ${MARIA_CTEST_LABEL_LA_UTILS} ${MARIA_CTEST_LABEL_D})

add_new_test(srandort
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_core maria_la_utils
    LABELS ${MARIA_CTEST_LABEL_LA_UTILS} ${MARIA_CTEST_LABEL_S})

add_new_test(drandort
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_core maria_la_utils
    LABELS ${MARIA_CTEST_LABEL_LA_UTILS} ${MARIA_CTEST_LABEL_D})

add_new_test(srandsvd
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_core maria_la_utils
    LABELS ${MARIA_CTEST_LABEL_LA_UTILS} ${MARIA_CTEST_LABEL_S})

add_new_test(drandsvd
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_core maria_la_utils
    LABELS ${MARIA_CTEST_LABEL_LA_UTILS} ${MARIA_CTEST_LABEL_D})

