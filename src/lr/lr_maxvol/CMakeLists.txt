add_library(maria_lr_maxvol
    STATIC lr_maxvol_mod.f90 ./sub/lr_maxvol_sub.f90)
target_link_libraries(maria_lr_maxvol
    PRIVATE maria_access_matrix maria_kinds maria_constants maria_argcheck maria_reports maria_la_core)
