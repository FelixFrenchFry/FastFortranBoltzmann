cmake_minimum_required(VERSION 3.24)

project(FastFortranBoltzmann
    VERSION 0.1.0
    LANGUAGES Fortran
)

find_package(Git QUIET)
find_package(MPI REQUIRED COMPONENTS Fortran)

if(GIT_FOUND)
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" rev-parse --short HEAD
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_VARIABLE FFB_GIT_COMMIT
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
endif()

if(NOT FFB_GIT_COMMIT)
    set(FFB_GIT_COMMIT "unknown")
endif()

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

# build feature options
option(FFB_FP64 "Usage of real64 for floating point values enabled" OFF)
option(FFB_DENSITY_CHECKS "Run-time density checks enabled" OFF)
set(FFB_EXPORT_NUM "" CACHE STRING "Export run directory name under output/")
set(FFB_SETTINGS_DEFINITIONS "" CACHE STRING "Fortran preprocessor overrides for source/settings.f90")

# executables in build/bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

set(FFB_SHARED_SOURCES
    source/domain.f90
    source/export.f90
    source/hardware_info.f90
    source/reporting.f90
    source/settings.f90
)

set(FFB_XYD_SHARED_SOURCES
    source/layouts/xyd/exchange.f90
    source/layouts/xyd/initialization.f90
    source/layouts/xyd/main.f90
    source/layouts/xyd/modes/couette_flow.f90
    source/layouts/xyd/modes/poiseuille_flow.f90
    source/layouts/xyd/modes/shear_wave.f90
    source/layouts/xyd/modes/sliding_lid.f90
)

set(FFB_XYD_SOURCES
    ${FFB_XYD_SHARED_SOURCES}
    source/layouts/xyd/simulation.f90
)

set(FFB_EXPERIMENT_SHARED_SOURCES source/experiments/exchange.f90 source/experiments/initialization.f90 source/experiments/main.f90)
string(TOUPPER "${CMAKE_BUILD_TYPE}" FFB_BUILD_TYPE_UPPER)
set(FFB_ACTIVE_FORTRAN_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${FFB_BUILD_TYPE_UPPER}}")
string(STRIP "${FFB_ACTIVE_FORTRAN_FLAGS}" FFB_ACTIVE_FORTRAN_FLAGS)

function(configure_ffb_target target_name dist_func_layout)
    set_target_properties(${target_name} PROPERTIES
        Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod/${target_name}"
    )

    # debug flags for ifx
    target_compile_options(${target_name} PRIVATE
        -coarray=distributed
        -fpp
        -warn all
        -standard-semantics
        $<$<CONFIG:Debug>:-O0>
        $<$<CONFIG:Debug>:-g>
        $<$<CONFIG:Debug>:-traceback>
        $<$<CONFIG:Debug>:-check all>
    )

    target_compile_definitions(${target_name} PRIVATE
        "FFB_COMPILER_ID=\"${CMAKE_Fortran_COMPILER_ID}\""
        "FFB_COMPILER_VERSION=\"${CMAKE_Fortran_COMPILER_VERSION}\""
        "FFB_FORTRAN_FLAGS=\"${FFB_ACTIVE_FORTRAN_FLAGS}\""
        "FFB_GIT_COMMIT=\"${FFB_GIT_COMMIT}\""
        "FFB_DIST_FUNC_LAYOUT=\"${dist_func_layout}\""
        ${FFB_SETTINGS_DEFINITIONS}
    )

    if(NOT FFB_EXPORT_NUM STREQUAL "")
        target_compile_definitions(${target_name} PRIVATE
            "FFB_EXPORT_NUM=\"${FFB_EXPORT_NUM}\""
        )
    endif()

    if(FFB_FP64)
        target_compile_definitions(${target_name} PRIVATE FFB_FP64)
    endif()
    if(FFB_DENSITY_CHECKS)
        target_compile_definitions(${target_name} PRIVATE FFB_DENSITY_CHECKS)
    endif()

    target_link_options(${target_name} PRIVATE
        -coarray=distributed
    )
endfunction()

add_executable(FFB ${FFB_SHARED_SOURCES} ${FFB_XYD_SOURCES})
configure_ffb_target(FFB xyd)
target_link_libraries(FFB PRIVATE MPI::MPI_Fortran)
add_executable(FFB_CHECKPOINT_0 ${FFB_SHARED_SOURCES} ${FFB_EXPERIMENT_SHARED_SOURCES} source/experiments/checkpoint_0/simulation.f90)
configure_ffb_target(FFB_CHECKPOINT_0 dxy)
target_link_libraries(FFB_CHECKPOINT_0 PRIVATE MPI::MPI_Fortran)
add_executable(FFB_CHECKPOINT_1 ${FFB_SHARED_SOURCES} ${FFB_EXPERIMENT_SHARED_SOURCES} source/experiments/checkpoint_1/simulation.f90)
configure_ffb_target(FFB_CHECKPOINT_1 dxy)
target_link_libraries(FFB_CHECKPOINT_1 PRIVATE MPI::MPI_Fortran)
add_executable(FFB_CHECKPOINT_2 ${FFB_SHARED_SOURCES} ${FFB_EXPERIMENT_SHARED_SOURCES} source/experiments/checkpoint_2/simulation.f90)
configure_ffb_target(FFB_CHECKPOINT_2 dxy)
target_link_libraries(FFB_CHECKPOINT_2 PRIVATE MPI::MPI_Fortran)
