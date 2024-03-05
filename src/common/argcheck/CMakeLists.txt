add_library(maria_argcheck
    STATIC argcheck_mod.f90 ./sub/argcheck_sub.f90)
target_link_libraries(maria_argcheck
    PRIVATE maria_kinds maria_comparison maria_reports)
