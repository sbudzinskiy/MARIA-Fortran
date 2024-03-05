set(MARIA_CTEST_LABEL_ARGCHECK maria_argcheck)

add_new_test(iarg_is_bad
    LINKS maria_assert maria_argcheck
    LABELS ${MARIA_CTEST_LABEL_ARGCHECK} ${MARIA_CTEST_LABEL_I})

add_new_test(sarg_is_bad
    LINKS maria_assert maria_constants maria_argcheck
    LABELS ${MARIA_CTEST_LABEL_ARGCHECK} ${MARIA_CTEST_LABEL_S})

add_new_test(darg_is_bad
    LINKS maria_assert maria_constants maria_argcheck
    LABELS ${MARIA_CTEST_LABEL_ARGCHECK} ${MARIA_CTEST_LABEL_D})
