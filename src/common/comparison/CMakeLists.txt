add_library(maria_comparison
    STATIC comparison_mod.f90 ./sub/comparison_sub.f90)
target_link_libraries(maria_comparison
    PRIVATE maria_kinds maria_constants)
