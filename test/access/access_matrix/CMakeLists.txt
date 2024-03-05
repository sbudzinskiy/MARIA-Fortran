set(MARIA_CTEST_LABEL_ACCESS_MATRIX maria_access_matrix)

add_new_test(smatval2slc
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_comparison maria_access_matrix
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_S})

add_new_test(dmatval2slc
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_comparison maria_access_matrix
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_D})

add_new_test(smatval2full
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_comparison maria_access_matrix
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_S})

add_new_test(dmatval2full
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_comparison maria_access_matrix
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_D})

