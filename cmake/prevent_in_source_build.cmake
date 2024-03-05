function(prevent_in_source_build)
    get_filename_component(SRCDIR ${PROJECT_SOURCE_DIR} REALPATH)
    get_filename_component(BINDIR ${PROJECT_BINARY_DIR} REALPATH)

    if(${BINDIR} STREQUAL ${SRCDIR})
        message(FATAL_ERROR "In-source build is prohibited. Use: cmake -B build")
    endif()
endfunction()

prevent_in_source_build()
