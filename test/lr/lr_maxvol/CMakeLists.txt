set(MARIA_CTEST_LABEL_LR_MAXVOL maria_lr_maxvol)

add_new_test(sgevolume
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_tsvd maria_access_matrix maria_lr_maxvol
    LABELS ${MARIA_CTEST_LABEL_LR_MAXVOL} ${MARIA_CTEST_LABEL_S})

add_new_test(dgevolume
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_tsvd maria_access_matrix maria_lr_maxvol
    LABELS ${MARIA_CTEST_LABEL_LR_MAXVOL} ${MARIA_CTEST_LABEL_D})

add_new_test(sgemaxvol_swap_rows
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_tsvd maria_access_matrix maria_lr_maxvol
    LABELS ${MARIA_CTEST_LABEL_LR_MAXVOL} ${MARIA_CTEST_LABEL_S})

add_new_test(dgemaxvol_swap_rows
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_tsvd maria_access_matrix maria_lr_maxvol
    LABELS ${MARIA_CTEST_LABEL_LR_MAXVOL} ${MARIA_CTEST_LABEL_D})

add_new_test(sgemaxvol
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_tsvd maria_access_matrix maria_lr_cross maria_lr_maxvol
    LABELS ${MARIA_CTEST_LABEL_LR_MAXVOL} ${MARIA_CTEST_LABEL_S})

add_new_test(dgemaxvol
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_tsvd maria_access_matrix maria_lr_cross maria_lr_maxvol
    LABELS ${MARIA_CTEST_LABEL_LR_MAXVOL} ${MARIA_CTEST_LABEL_D})

add_new_test(sgemaxvol_rect_add_rows
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_tsvd maria_access_matrix maria_lr_maxvol
    LABELS ${MARIA_CTEST_LABEL_LR_MAXVOL} ${MARIA_CTEST_LABEL_S})

add_new_test(dgemaxvol_rect_add_rows
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_tsvd maria_access_matrix maria_lr_maxvol
    LABELS ${MARIA_CTEST_LABEL_LR_MAXVOL} ${MARIA_CTEST_LABEL_D})

add_new_test(sgemaxvol_rect_swap_rows
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_tsvd maria_access_matrix maria_lr_maxvol
    LABELS ${MARIA_CTEST_LABEL_LR_MAXVOL} ${MARIA_CTEST_LABEL_S})

add_new_test(dgemaxvol_rect_swap_rows
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_tsvd maria_access_matrix maria_lr_maxvol
    LABELS ${MARIA_CTEST_LABEL_LR_MAXVOL} ${MARIA_CTEST_LABEL_D})

add_new_test(sgemaxvol_rect_swap_cols
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_tsvd maria_access_matrix maria_lr_maxvol
    LABELS ${MARIA_CTEST_LABEL_LR_MAXVOL} ${MARIA_CTEST_LABEL_S})

add_new_test(dgemaxvol_rect_swap_cols
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_tsvd maria_access_matrix maria_lr_maxvol
    LABELS ${MARIA_CTEST_LABEL_LR_MAXVOL} ${MARIA_CTEST_LABEL_D})

add_new_test(sgemaxvol_proj
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_tsvd maria_access_matrix maria_lr_cross maria_lr_maxvol
    LABELS ${MARIA_CTEST_LABEL_LR_MAXVOL} ${MARIA_CTEST_LABEL_S})

add_new_test(dgemaxvol_proj
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_tsvd maria_access_matrix maria_lr_cross maria_lr_maxvol
    LABELS ${MARIA_CTEST_LABEL_LR_MAXVOL} ${MARIA_CTEST_LABEL_D})

