set(MARIA_CTEST_LABEL_TT_UTILS maria_tt_utils)

add_new_test(ttrank
    LINKS maria_assert maria_tt_utils ${MARIA_PRNG_IMPLEMENTATION} 
    LABELS ${MARIA_CTEST_LABEL_TT_UTILS})

add_new_test(stt2full
    LINKS maria_assert maria_kinds maria_constants maria_la_utils maria_access_tensor maria_tt_utils
    LABELS ${MARIA_CTEST_LABEL_TT_UTILS} ${MARIA_CTEST_LABEL_S})

add_new_test(dtt2full
    LINKS maria_assert maria_kinds maria_constants maria_la_utils maria_access_tensor maria_tt_utils
    LABELS ${MARIA_CTEST_LABEL_TT_UTILS} ${MARIA_CTEST_LABEL_D})

#add_new_test(sgevolume
#    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_tsvd maria_access_matrix maria_lr_maxvol
#    LABELS ${MARIA_CTEST_LABEL_LR_MAXVOL} ${MARIA_CTEST_LABEL_S})
#
#add_new_test(dgevolume
#    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core maria_lr_tsvd maria_access_matrix maria_lr_maxvol
#    LABELS ${MARIA_CTEST_LABEL_LR_MAXVOL} ${MARIA_CTEST_LABEL_D})
