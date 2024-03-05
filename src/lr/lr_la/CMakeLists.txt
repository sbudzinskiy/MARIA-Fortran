add_library(maria_lr_la
    STATIC lr_la_mod.f90 ./sub/lr_la_sub.f90)
target_link_libraries(maria_lr_la
    PRIVATE maria_kinds maria_constants maria_argcheck maria_reports maria_la_core)
