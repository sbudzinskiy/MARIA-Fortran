add_library(maria_tt_utils
    STATIC tt_utils_mod.f90 ./sub/tt_utils_sub.f90)
target_link_libraries(maria_tt_utils
    PRIVATE maria_access_matrix maria_kinds maria_constants maria_argcheck maria_reports maria_la_core)
