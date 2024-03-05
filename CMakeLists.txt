# -- CMake configuration -------------------------------------------------------
cmake_minimum_required(VERSION 3.24)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# -- Project description -------------------------------------------------------
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/VERSION PROJECT_VERSION)
project(MARIA 
    LANGUAGES Fortran 
    VERSION ${PROJECT_VERSION})

# -- Options -------------------------------------------------------------------
option(MARIA_USE_MKL "Use MKL?" ON)
option(MARIA_MKL_SEQ "Use sequential version of MKL?" ON)
option(MARIA_USE_TESTS "Build unit tests?" ON)

# -- Languages -----------------------------------------------------------------
enable_language(Fortran)
if(MARIA_USE_MKL)
    message(STATUS "Enabling CXX to build with MKL...")
    enable_language(CXX)
endif()

# -- Set up build and additional directories -----------------------------------
include(prevent_in_source_build)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_Fortran_MODULE_DIRECTORY})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# -- Set up header files with macros -------------------------------------------
set(MACROS_HEADER_FILES
    ${PROJECT_SOURCE_DIR}/src/common/assert/maria_assert_mod.h)
file(COPY ${MACROS_HEADER_FILES} 
    DESTINATION ${PROJECT_BINARY_DIR}/include)
include_directories(${PROJECT_BINARY_DIR}/include)

# -- Set up linear algebra core ------------------------------------------------
if(MARIA_USE_MKL)
    if(MARIA_MKL_SEQ)
        set(BLA_VENDOR Intel10_64lp_seq)
    else()
        set(BLA_VENDOR Intel10_64lp)
    endif()

    if(DEFINED ENV{MKLROOT})
        find_package(MKL REQUIRED)
        #        add_compile_definitions(MARIA_MKL)
    else()
        message(WARNING "Intel MKL libraries must be configured beforehand:"
                        "\nbash> source $(PATH_TO_MKL)/bin/mklvars.sh intel64"
                        "\nor"
                        "\nbash> source $(PATH_TO_ONEAPI)/setvars.sh intel64")
        message(FATAL_ERROR "Intel MKL was not found.")
    endif()

    set(BLAS_LINK ${MKL_LINK})
    set(LAPACK_LINK ${MKL_LINK})
else()
    set(BLA_VENDOR Generic)
    find_package(BLAS)
    find_package(LAPACK)
    set(BLAS_LINK BLAS::BLAS)
    set(LAPACK_LINK LAPACK::LAPACK)
endif()

# -- Source files --------------------------------------------------------------
add_subdirectory(src)

# -- Set up unit tests ---------------------------------------------------------
if(MARIA_USE_TESTS)
    include(CTest)
    enable_testing()
    include(add_new_test)
    if(BUILD_TESTING)
        add_subdirectory(test)
    endif()
endif()
