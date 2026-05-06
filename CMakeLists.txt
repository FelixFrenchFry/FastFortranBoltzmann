cmake_minimum_required(VERSION 3.24)

project(FastFortranBoltzmann
    VERSION 0.1.0
    LANGUAGES Fortran
)

# only support for Intel Fortran compiler
if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    message(FATAL_ERROR
    "Usage of the 'ifx' compiler is required, configure with -DCMAKE_Fortran_COMPILER=ifx"
  )
endif()

# only support for modern Fortran
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_EXTENSIONS OFF)

# compiler flags for different features
option(FFB_FP64 "Usage of real64 for floating point values enabled" OFF)
option(FFB_DENSITY_CHECKS "Run-time density checks enabled" OFF)
option(FFB_BOUNDARY_CHECKS "Run-time boundary condition checks enabled" OFF)

# executables in build/bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(FFB
    source/domain.f90
    source/export.f90
    source/initialization.f90
    source/main.f90
    source/settings.f90
    source/simulation.f90
    source/sim/couette_flow.f90
    source/sim/poiseuille_flow.f90
    source/sim/shear_wave.f90
    source/sim/sliding_lid.f90
)

# debug flags for ifx
target_compile_options(FFB PRIVATE
    -coarray
    -fpp
    -warn all
    -standard-semantics
    $<$<CONFIG:Debug>:-O0>
    $<$<CONFIG:Debug>:-g>
    $<$<CONFIG:Debug>:-traceback>
    $<$<CONFIG:Debug>:-check all>
)

if(FFB_FP64)
    target_compile_definitions(FFB PRIVATE FFB_FP64)
endif()
if(FFB_DENSITY_CHECKS)
    target_compile_definitions(FFB PRIVATE FFB_DENSITY_CHECKS)
endif()
if(FFB_BOUNDARY_CHECKS)
    target_compile_definitions(FFB PRIVATE FFB_BOUNDARY_CHECKS)
endif()

target_link_options(FFB PRIVATE
    -coarray
)
