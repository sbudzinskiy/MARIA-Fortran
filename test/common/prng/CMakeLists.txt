set(MARIA_CTEST_LABEL_PRNG maria_prng)

add_new_test(snormal
    LINKS maria_assert maria_constants maria_comparison maria_kinds ${MARIA_PRNG_IMPLEMENTATION} maria_utils
    LABELS ${MARIA_CTEST_LABEL_PRNG} ${MARIA_CTEST_LABEL_S})

add_new_test(dnormal
    LINKS maria_assert maria_constants maria_comparison maria_kinds ${MARIA_PRNG_IMPLEMENTATION} maria_utils
    LABELS ${MARIA_CTEST_LABEL_PRNG} ${MARIA_CTEST_LABEL_D})

add_new_test(iuniform
    LINKS maria_assert ${MARIA_PRNG_IMPLEMENTATION}
    LABELS ${MARIA_CTEST_LABEL_PRNG} ${MARIA_CTEST_LABEL_I})

add_new_test(suniform
    LINKS maria_assert maria_constants maria_comparison maria_kinds ${MARIA_PRNG_IMPLEMENTATION} maria_utils
    LABELS ${MARIA_CTEST_LABEL_PRNG} ${MARIA_CTEST_LABEL_S})

add_new_test(duniform
    LINKS maria_assert maria_constants maria_comparison maria_kinds ${MARIA_PRNG_IMPLEMENTATION} maria_utils
    LABELS ${MARIA_CTEST_LABEL_PRNG} ${MARIA_CTEST_LABEL_D})
