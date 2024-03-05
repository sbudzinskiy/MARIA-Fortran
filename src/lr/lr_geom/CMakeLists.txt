add_library(maria_lr_geom
    STATIC lr_geom_mod.f90 ./sub/lr_geom_sub.f90)
target_link_libraries(maria_lr_geom
    PRIVATE maria_kinds maria_constants maria_comparison maria_argcheck maria_reports maria_la_core)
