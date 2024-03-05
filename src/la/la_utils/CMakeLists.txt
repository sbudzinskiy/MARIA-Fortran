add_library(maria_la_utils
    STATIC la_utils_mod.f90 ./sub/la_utils_sub.f90)
target_link_libraries(maria_la_utils
    PRIVATE maria_kinds maria_constants maria_argcheck maria_reports maria_utils maria_la_core)
