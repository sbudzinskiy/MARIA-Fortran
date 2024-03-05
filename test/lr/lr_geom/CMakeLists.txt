set(MARIA_CTEST_LABEL_LR_GEOM maria_lr_geom)

add_new_test(slrproj_tangent
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_geom
    LABELS ${MARIA_CTEST_LABEL_LR_GEOM} ${MARIA_CTEST_LABEL_S})

add_new_test(dlrproj_tangent
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_geom
    LABELS ${MARIA_CTEST_LABEL_LR_GEOM} ${MARIA_CTEST_LABEL_D})

add_new_test(slrdotf_tangent
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_geom
    LABELS ${MARIA_CTEST_LABEL_LR_GEOM} ${MARIA_CTEST_LABEL_S})

add_new_test(dlrdotf_tangent
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_geom
    LABELS ${MARIA_CTEST_LABEL_LR_GEOM} ${MARIA_CTEST_LABEL_D})

add_new_test(slrretr_tangent
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_geom
    LABELS ${MARIA_CTEST_LABEL_LR_GEOM} ${MARIA_CTEST_LABEL_S})

add_new_test(dlrretr_tangent
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_geom
    LABELS ${MARIA_CTEST_LABEL_LR_GEOM} ${MARIA_CTEST_LABEL_D})

