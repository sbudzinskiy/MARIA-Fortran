add_library(maria_access_matrix
    STATIC access_matrix_mod.f90 ./sub/access_matrix_sub.f90)
target_link_libraries(maria_access_matrix
    PRIVATE maria_kinds maria_constants maria_argcheck maria_reports)
