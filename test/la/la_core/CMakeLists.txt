set(MARIA_CTEST_LABEL_LA_CORE maria_la_core)

add_new_test(sgesdd_q
    LINKS maria_assert maria_kinds maria_constants maria_la_core
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_S})

add_new_test(dgesdd_q
    LINKS maria_assert maria_kinds maria_constants maria_la_core
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_D})

add_new_test(sroundup_lwork
    LINKS maria_assert maria_kinds maria_la_core
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_S})

add_new_test(droundup_lwork
    LINKS maria_assert maria_kinds maria_la_core
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_D})

add_new_test(sgescal
    LINKS maria_assert maria_kinds maria_constants maria_la_utils maria_la_core
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_S})

add_new_test(dgescal
    LINKS maria_assert maria_kinds maria_constants maria_la_utils maria_la_core
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_D})

add_new_test(sgedotf
    LINKS maria_assert maria_kinds maria_constants maria_utils maria_la_core
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_S})

add_new_test(dgedotf
    LINKS maria_assert maria_kinds maria_constants maria_utils maria_la_core
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_D})

add_new_test(sgenrmf
    LINKS maria_assert maria_kinds maria_constants maria_utils ${MARIA_PRNG_IMPLEMENTATION} maria_la_core
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_S})

add_new_test(dgenrmf
    LINKS maria_assert maria_kinds maria_constants maria_utils ${MARIA_PRNG_IMPLEMENTATION} maria_la_core
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_D})

add_new_test(sgenrmc
    LINKS maria_assert maria_kinds maria_constants maria_utils ${MARIA_PRNG_IMPLEMENTATION} maria_la_core
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_S})

add_new_test(dgenrmc
    LINKS maria_assert maria_kinds maria_constants maria_utils ${MARIA_PRNG_IMPLEMENTATION} maria_la_core
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_D})

add_new_test(sdgmm
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_S})

add_new_test(ddgmm
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_D})

add_new_test(sgepiv
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_utils maria_la_utils maria_la_core
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_S})

add_new_test(dgepiv
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_utils maria_la_utils maria_la_core
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_D})

add_new_test(sorcangles
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_S})

add_new_test(dorcangles
    LINKS maria_assert maria_kinds maria_constants ${MARIA_PRNG_IMPLEMENTATION} maria_la_utils maria_la_core
    LABELS ${MARIA_CTEST_LABEL_LA_CORE} ${MARIA_CTEST_LABEL_D})

