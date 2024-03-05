add_library(maria_access_tensor
    STATIC access_tensor_mod.f90 ./sub/access_tensor_sub.f90)
target_link_libraries(maria_access_tensor
    PRIVATE maria_kinds maria_constants maria_argcheck maria_reports)
