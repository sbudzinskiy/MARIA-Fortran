# **MARIA**: **MA**trix and tenso**R** **I**nterpolation and **A**pproximation

## Purpose

**MARIA** is an open-source Fortran library for working with low-rank matrices and tensors:

- truncated SVD
- truncated SVD with sketching
- adaptive cross approximation with rook pivoting
- cross approximation based on the maximum-volume principle
- cross approximation based on the maximum-projective-volume principle
- elements of Riemannian geometry of fixed-rank matrices
- tensor-train SVD

**Note:** **MARIA** will not be developed further as a Fortran library, it will be ported to Julia

## Style

**MARIA** is inspired by [LAPACK](https://netlib.org/lapack/) and provides

- single- and double-precision versions of every subroutine
- maximum control over memory allocation

## Technical features

- written in modern Fortran 2008 
- source code structured with modules and submodules
- automated build with [CMake](https://cmake.org/) and 6 configuration presets
    - debug or release mode
    - [gfortran](https://gcc.gnu.org/wiki/GFortran) or [ifort](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html) compiler
    - [OpenBLAS](https://www.openblas.net/) or [Intel MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html)
- automated testing with [CTest](https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html)
- source code documentation with [FORD](https://github.com/Fortran-FOSS-Programmers/ford)

## How to build and test

The build and testing processes for **MARIA** are automated with CMake and CTest, respectively.

1. Navigate to the root folder of **MARIA** (`maria`, by default)
    ```console
    foo@bar:~$ cd <path_to_maria>/maria/
    ```
2. Choose one of the 6 available configuration presets

    | Preset  | Compiler| Mode | LAPACK |
    |---------|---------|------|--------|
    | `gnu-debug-lapack`  | `gfortran`  | Debug | ObenBLAS |
    | `intel-debug-mkl`  | `ifort`  | Debug | Intel MKL (parallel) |
    | `intel-debug-mklseq`  | `ifort`  | Debug | Intel MKL (sequential) |
    | `gnu-release-lapack`  | `gfortran`  | Release | ObenBLAS |
    | `intel-release-mkl`  | `ifort`  | Release | Intel MKL (parallel) |
    | `intel-release-mklseq`  | `ifort`  | Release | Intel MKL (sequential) |

    **Note:** when using one of the Intel presets, make sure that `ifort` and Intel MKL are available by running
    ```console
    foo@bar:~$ source /opt/intel/oneapi/setvars.sh intel64
    ```
    or
    ```console
    foo@bar:~$ source /opt/intel/mkl/mklvars.sh intel64
    ```
    In addition, it may be required to change the `-mkl` compiler flag in `CMakePresetes.json` to `-qmkl`.

3. Build the code by typing
    ```console
    foo@bar:~<path_to_maria>/maria$ cmake -S . --preset <preset_name>
    foo@bar:~<path_to_maria>/maria$ cmake --build --preset <preset_name>
    ```
    This will create a folder `<path_to_maria>/maria/build/<preset_name>` and put the corresponding binaries there.

4. You can check that the code is built correctly by running the unit tests:
    ```console
    foo@bar:~<path_to_maria>/maria$ ctest --preset <preset_name>
    ```

## Documentation

You can generate the source code documentation with [FORD](https://github.com/Fortran-FOSS-Programmers/ford) and view it as an HTML page. FORD can be installed with
```console
pip install ford
```
Then, go to the documentation folder and run FORD:
```console
foo@bar:~<path_to_maria>/maria$ cd ./doc 
foo@bar:~<path_to_maria>/maria/doc$ ford ./ford.md
```
This will create a `ford_html` folder with `index.html` inside.
