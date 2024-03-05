add_library(maria_lr_cross
    STATIC lr_cross_mod.f90 ./sub/lr_cross_sub.f90)
target_link_libraries(maria_lr_cross
    PRIVATE maria_kinds maria_constants maria_argcheck maria_reports maria_la_core maria_lr_la)
