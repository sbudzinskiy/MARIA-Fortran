set(MARIA_CTEST_LABEL_COMPARISON maria_comparison)

add_new_test(ssafe_less
    LINKS maria_assert maria_constants maria_comparison
    LABELS ${MARIA_CTEST_LABEL_COMPARISON} ${MARIA_CTEST_LABEL_S})

add_new_test(dsafe_less
    LINKS maria_assert maria_constants maria_comparison
    LABELS ${MARIA_CTEST_LABEL_COMPARISON} ${MARIA_CTEST_LABEL_D})

add_new_test(ssafe_leq
    LINKS maria_assert maria_constants maria_comparison
    LABELS ${MARIA_CTEST_LABEL_COMPARISON} ${MARIA_CTEST_LABEL_S})

add_new_test(dsafe_leq
    LINKS maria_assert maria_constants maria_comparison
    LABELS ${MARIA_CTEST_LABEL_COMPARISON} ${MARIA_CTEST_LABEL_D})

add_new_test(ssafe_eq
    LINKS maria_assert maria_constants maria_comparison
    LABELS ${MARIA_CTEST_LABEL_COMPARISON} ${MARIA_CTEST_LABEL_S})

add_new_test(dsafe_eq
    LINKS maria_assert maria_constants maria_comparison
    LABELS ${MARIA_CTEST_LABEL_COMPARISON} ${MARIA_CTEST_LABEL_D})
