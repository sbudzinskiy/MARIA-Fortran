set(MARIA_CTEST_LABEL_LR_LA maria_lr_la)

add_new_test(slrval
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_utils maria_la_core maria_lr_la
    LABELS ${MARIA_CTEST_LABEL_LR_LA} ${MARIA_CTEST_LABEL_S})

add_new_test(dlrval
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_utils maria_la_core maria_lr_la
    LABELS ${MARIA_CTEST_LABEL_LR_LA} ${MARIA_CTEST_LABEL_D})

add_new_test(slr2full
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_utils maria_la_utils maria_la_core maria_lr_la
    LABELS ${MARIA_CTEST_LABEL_LR_LA} ${MARIA_CTEST_LABEL_S})

add_new_test(dlr2full
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_utils maria_la_utils maria_la_core maria_lr_la
    LABELS ${MARIA_CTEST_LABEL_LR_LA} ${MARIA_CTEST_LABEL_D})

add_new_test(slrdotf
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_utils maria_la_core maria_lr_la
    LABELS ${MARIA_CTEST_LABEL_LR_LA} ${MARIA_CTEST_LABEL_S})

add_new_test(dlrdotf
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_utils maria_la_core maria_lr_la
    LABELS ${MARIA_CTEST_LABEL_LR_LA} ${MARIA_CTEST_LABEL_D})

add_new_test(slrnrmf
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_utils maria_la_core maria_lr_la
    LABELS ${MARIA_CTEST_LABEL_LR_LA} ${MARIA_CTEST_LABEL_S})

add_new_test(dlrnrmf
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_utils maria_la_core maria_lr_la
    LABELS ${MARIA_CTEST_LABEL_LR_LA} ${MARIA_CTEST_LABEL_D})

add_new_test(slrnrmf_diff
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_utils maria_la_core maria_lr_la
    LABELS ${MARIA_CTEST_LABEL_LR_LA} ${MARIA_CTEST_LABEL_S})

add_new_test(dlrnrmf_diff
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_utils maria_la_core maria_lr_la
    LABELS ${MARIA_CTEST_LABEL_LR_LA} ${MARIA_CTEST_LABEL_D})

add_new_test(lrort_rank
    LINKS maria_assert maria_lr_la
    LABELS ${MARIA_CTEST_LABEL_LR_LA} ${MARIA_CTEST_LABEL_I})

add_new_test(slrort
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_la
    LABELS ${MARIA_CTEST_LABEL_LR_LA} ${MARIA_CTEST_LABEL_S})

add_new_test(dlrort
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_la
    LABELS ${MARIA_CTEST_LABEL_LR_LA} ${MARIA_CTEST_LABEL_D})

add_new_test(slrsvd_ort
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_la
    LABELS ${MARIA_CTEST_LABEL_LR_LA} ${MARIA_CTEST_LABEL_S})

add_new_test(dlrsvd_ort
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_la
    LABELS ${MARIA_CTEST_LABEL_LR_LA} ${MARIA_CTEST_LABEL_D})

