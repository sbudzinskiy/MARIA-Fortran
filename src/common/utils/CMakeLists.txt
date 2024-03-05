add_library(maria_utils
    STATIC utils_mod.f90 ./sub/utils_sub.f90)
target_link_libraries(maria_utils
    PRIVATE maria_kinds maria_constants maria_argcheck maria_reports maria_prng)
