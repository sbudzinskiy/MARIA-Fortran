add_library(maria_tt_tsvd
    STATIC tt_tsvd_mod.f90 ./sub/tt_tsvd_sub.f90)
target_link_libraries(maria_tt_tsvd
    PRIVATE maria_access_matrix maria_kinds maria_constants maria_argcheck maria_reports maria_la_core maria_lr_tsvd)
