add_library(maria_constants
    STATIC constants_mod.f90)
target_link_libraries(maria_constants
    PRIVATE maria_kinds)
