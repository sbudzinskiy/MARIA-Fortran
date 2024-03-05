add_library(maria_lr_tsvd
    STATIC lr_tsvd_mod.f90 ./sub/lr_tsvd_sub.f90)
target_link_libraries(maria_lr_tsvd
    PRIVATE maria_kinds maria_constants maria_comparison maria_argcheck maria_reports maria_la_core)
