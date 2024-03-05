set(MARIA_CTEST_LABEL_TT_TSVD maria_tt_tsvd)

add_new_test(sttsvd
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_lr_tsvd maria_access_tensor maria_tt_utils maria_tt_tsvd
    LABELS ${MARIA_CTEST_LABEL_TT_TSVD} ${MARIA_CTEST_LABEL_S})

add_new_test(dttsvd
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_lr_tsvd maria_access_tensor maria_tt_utils maria_tt_tsvd
    LABELS ${MARIA_CTEST_LABEL_TT_TSVD} ${MARIA_CTEST_LABEL_D})

