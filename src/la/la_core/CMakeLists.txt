add_library(maria_la_core
    STATIC la_core_mod.f90 ./sub/la_core_sub.f90)
target_link_libraries(maria_la_core
    PRIVATE maria_kinds maria_constants maria_argcheck maria_reports maria_utils
            maria_blas1 maria_blas2 maria_blas3 maria_lapack)

add_library(maria_blas1
    STATIC blas1_mod.f90)
target_link_libraries(maria_blas1
    PRIVATE maria_kinds ${BLAS_LINK})

add_library(maria_blas2
    STATIC blas2_mod.f90)
target_link_libraries(maria_blas2
    PRIVATE maria_kinds ${BLAS_LINK})

add_library(maria_blas3
    STATIC blas3_mod.f90)
target_link_libraries(maria_blas3
    PRIVATE maria_kinds ${BLAS_LINK})

add_library(maria_lapack
    STATIC lapack_mod.f90)
target_link_libraries(maria_lapack
    PRIVATE maria_kinds ${LAPACK_LINK})
