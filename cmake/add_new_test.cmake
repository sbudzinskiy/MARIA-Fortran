function(add_new_test NAME)
    cmake_parse_arguments(TEST "" "" "LINKS;LABELS" ${ARGN})

    set_source_files_properties(${NAME}_test.f90
        PROPERTIES Fortran_PREPROCESS ON)

    add_executable(${NAME}_test
        ${NAME}_test.f90)
    target_link_libraries(${NAME}_test ${TEST_LINKS})
    add_test(
        NAME ${NAME}
        COMMAND $<TARGET_FILE:${NAME}_test> ${MARIA_TEST_SEED})
    set_property(
        TEST ${NAME}
        PROPERTY LABELS ${TEST_LABELS})
endfunction()
