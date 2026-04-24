cmake_minimum_required(VERSION 3.24)

project(FastFortranBoltzmann
    VERSION 0.1.0
    LANGUAGES Fortran
)

# only support for Intel Fortran compiler
if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    message(FATAL_ERROR
    "This project requires Intel ifx, configured with -DCMAKE_Fortran_COMPILER=ifx"
  )
endif()

# only support for modern Fortran
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_EXTENSIONS OFF)

# executables in build/bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(FFB
    source/initialization.f90
    source/main.f90
)

# coarrays enabled and debug-friendly flags
target_compile_options(FFB PRIVATE
  -coarray
  $<$<CONFIG:Debug>:-O0>
  $<$<CONFIG:Debug>:-g>
  $<$<CONFIG:Debug>:-traceback>
  $<$<CONFIG:Debug>:-warn all>
  $<$<CONFIG:Debug>:-check all>
  $<$<CONFIG:Debug>:-standard-semantics>
  $<$<CONFIG:Release>:-O3>
)

target_link_options(FFB PRIVATE
  -coarray
)
