set(MARIA_CTEST_LABEL_LR_TSVD maria_lr_tsvd)

add_new_test(schop
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_comparison maria_utils maria_lr_tsvd
    LABELS ${MARIA_CTEST_LABEL_LR_TSVD} ${MARIA_CTEST_LABEL_S})

add_new_test(dchop
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_comparison maria_utils maria_lr_tsvd
    LABELS ${MARIA_CTEST_LABEL_LR_TSVD} ${MARIA_CTEST_LABEL_D})

add_new_test(sgersvd1
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_tsvd maria_access_matrix
    LABELS ${MARIA_CTEST_LABEL_LR_TSVD} ${MARIA_CTEST_LABEL_S})

add_new_test(dgersvd1
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_tsvd maria_access_matrix
    LABELS ${MARIA_CTEST_LABEL_LR_TSVD} ${MARIA_CTEST_LABEL_D})

add_new_test(sgersvd2
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_tsvd maria_access_matrix
    LABELS ${MARIA_CTEST_LABEL_LR_TSVD} ${MARIA_CTEST_LABEL_S})

add_new_test(dgersvd2
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_tsvd maria_access_matrix
    LABELS ${MARIA_CTEST_LABEL_LR_TSVD} ${MARIA_CTEST_LABEL_D})

