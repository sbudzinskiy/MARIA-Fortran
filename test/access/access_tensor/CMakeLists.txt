set(MARIA_CTEST_LABEL_ACCESS_TENSOR maria_access_tensor)

add_new_test(mi2i
    LINKS maria_assert maria_kinds maria_constants maria_comparison maria_access_tensor ${MARIA_PRNG_IMPLEMENTATION}
    LABELS ${MARIA_CTEST_LABEL_ACCESS_TENSOR})

add_new_test(i2mi
    LINKS maria_assert maria_kinds maria_constants maria_comparison maria_access_tensor ${MARIA_PRNG_IMPLEMENTATION}
    LABELS ${MARIA_CTEST_LABEL_ACCESS_TENSOR})

add_new_test(stenval2fib
    LINKS maria_assert maria_kinds maria_constants maria_comparison maria_access_tensor ${MARIA_PRNG_IMPLEMENTATION}
    LABELS ${MARIA_CTEST_LABEL_ACCESS_TENSOR} ${MARIA_CTEST_LABEL_S})

add_new_test(dtenval2fib
    LINKS maria_assert maria_kinds maria_constants maria_comparison maria_access_tensor ${MARIA_PRNG_IMPLEMENTATION}
    LABELS ${MARIA_CTEST_LABEL_ACCESS_TENSOR} ${MARIA_CTEST_LABEL_D})

add_new_test(stenval2full
    LINKS maria_assert maria_kinds maria_constants maria_comparison maria_access_tensor ${MARIA_PRNG_IMPLEMENTATION}
    LABELS ${MARIA_CTEST_LABEL_ACCESS_TENSOR} ${MARIA_CTEST_LABEL_S})

add_new_test(dtenval2full
    LINKS maria_assert maria_kinds maria_constants maria_comparison maria_access_tensor ${MARIA_PRNG_IMPLEMENTATION}
    LABELS ${MARIA_CTEST_LABEL_ACCESS_TENSOR} ${MARIA_CTEST_LABEL_D})

