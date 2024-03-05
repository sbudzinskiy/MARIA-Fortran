set(MARIA_CTEST_LABEL_LR_CROSS maria_lr_cross)

add_new_test(smatcross_top
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_cross maria_access_matrix
    LABELS ${MARIA_CTEST_LABEL_LR_CROSS} ${MARIA_CTEST_LABEL_S})

add_new_test(dmatcross_top
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_cross maria_access_matrix
    LABELS ${MARIA_CTEST_LABEL_LR_CROSS} ${MARIA_CTEST_LABEL_D})

add_new_test(smatcross
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_cross maria_access_matrix
    LABELS ${MARIA_CTEST_LABEL_LR_CROSS} ${MARIA_CTEST_LABEL_S})

add_new_test(dmatcross
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_cross maria_access_matrix
    LABELS ${MARIA_CTEST_LABEL_LR_CROSS} ${MARIA_CTEST_LABEL_D})

add_new_test(smatcross_aca
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_cross maria_access_matrix
    LABELS ${MARIA_CTEST_LABEL_LR_CROSS} ${MARIA_CTEST_LABEL_S})

add_new_test(dmatcross_aca
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_cross maria_access_matrix
    LABELS ${MARIA_CTEST_LABEL_LR_CROSS} ${MARIA_CTEST_LABEL_D})
