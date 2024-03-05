add_library(maria_prng
    STATIC prng_mod.f90)
target_link_libraries(maria_prng
    PRIVATE maria_kinds)

if(MARIA_USE_MKL)
    add_library(mkl_vsl_type
        STATIC ${MKL_INCLUDE_DIRS}/mkl_vsl.f90)

    add_library(maria_prng_mkl
        STATIC prng_mkl_mod.f90 ./sub/prng_mkl_sub.f90)
    target_link_libraries(maria_prng_mkl
        PRIVATE ${MKL_LINK} mkl_vsl_type maria_kinds maria_constants maria_argcheck maria_reports maria_prng)
    target_compile_definitions(maria_prng_mkl
        INTERFACE MARIA_MKL)
else()
    add_library(maria_prng_builtin
        STATIC prng_builtin_mod.f90 ./sub/prng_builtin_sub.f90)
    target_link_libraries(maria_prng_builtin
        PRIVATE maria_kinds maria_constants maria_argcheck maria_reports maria_prng)
endif()
