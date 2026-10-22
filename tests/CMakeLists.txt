option(FFB_TEST_EXTENDED "Include four-image decomposition tests" OFF)

# compile the production modules with isolated module files
set(FFB_TEST_SOURCES ${FFB_SHARED_SOURCES} ${FFB_XYD_SOURCES})
list(REMOVE_ITEM FFB_TEST_SOURCES source/layouts/xyd/main.f90)
list(TRANSFORM FFB_TEST_SOURCES PREPEND "${PROJECT_SOURCE_DIR}/")
set(FFB_TEST_ROOT "${CMAKE_CURRENT_BINARY_DIR}/runs")
set(FFB_TEST_ENVIRONMENT
    "I_MPI_PIN=0"
    "I_MPI_FABRICS=shm"
    "I_MPI_HYDRA_BOOTSTRAP=fork"
)
set(FFB_TEST_ENVIRONMENT_MODIFICATIONS
    "FOR_COARRAY_CONFIG_FILE=unset:"
    "I_MPI_HYDRA_HOST_FILE=unset:"
    "I_MPI_PIN_DOMAIN=unset:"
    "I_MPI_PIN_ORDER=unset:"
    "I_MPI_PIN_PROCESSOR_LIST=unset:"
)

function(add_ffb_test_core target_name sim_mode coarray_mode)
    # local overrides keep tests independent of benchmark configuration
    set(FFB_EXPORT_NUM "")
    set(FFB_SETTINGS_DEFINITIONS
        FFB_USE_CMAKE_SETTINGS
        FFB_N_X=64 FFB_N_Y=48 FFB_N_STEPS=80 FFB_SIM_MODE=${sim_mode}
        FFB_RHO_0=1.0_FP FFB_OMEGA=1.1_FP FFB_U_MAX=0.01_FP FFB_N_SIN=1.0_FP
        FFB_U_WALL=0.02_FP FFB_U_LID=0.02_FP
        FFB_RHO_IN=1.00005_FP FFB_RHO_OUT=0.99995_FP
        FFB_USE_UNROLLED_KERNELS=.true.
        FFB_EXPORT_MACROS=.true. FFB_EXPORT_ENDPOINT_STATES=.true.
        FFB_EXPORT_INTERVAL=20 "FFB_EXPORT_NUM=\"test\""
        FFB_INTERACTIVE_PROGRESS=.false. FFB_PROGRESS_INTERVAL=80
    )
    add_library(${target_name} OBJECT EXCLUDE_FROM_ALL ${FFB_TEST_SOURCES})
    configure_ffb_target(${target_name} xyd ${coarray_mode})
    target_include_directories(${target_name} PUBLIC "${CMAKE_BINARY_DIR}/mod/${target_name}")
    target_link_libraries(${target_name} PUBLIC MPI::MPI_Fortran)
endfunction()

function(add_ffb_test_solver target_name sim_mode)
    add_ffb_test_core(${target_name}_core ${sim_mode} distributed)
    add_executable(${target_name} EXCLUDE_FROM_ALL "${PROJECT_SOURCE_DIR}/source/layouts/xyd/main.f90")
    configure_ffb_target(${target_name} xyd)
    target_link_libraries(${target_name} PRIVATE ${target_name}_core)
    add_dependencies(ffb_tests ${target_name})
endfunction()

function(configure_ffb_test name images ix iy labels)
    set_tests_properties(${name} PROPERTIES
        LABELS "${labels}"
        PROCESSORS ${images}
        TIMEOUT 40
        FAIL_REGULAR_EXPRESSION "forrtl: warning [(]406[)]"
        ENVIRONMENT "${FFB_TEST_ENVIRONMENT};FOR_COARRAY_NUM_IMAGES=${images};I_X=${ix};I_Y=${iy}"
        ENVIRONMENT_MODIFICATION "${FFB_TEST_ENVIRONMENT_MODIFICATIONS}"
    )
endfunction()

function(add_ffb_run name target_name run_dir images ix iy fixture labels)
    file(MAKE_DIRECTORY "${run_dir}")
    add_test(NAME ${name}
        COMMAND "${CMAKE_COMMAND}"
            "-DFFB_TEST_EXE=$<TARGET_FILE:${target_name}>"
            "-DFFB_TEST_ROOT=${FFB_TEST_ROOT}"
            "-DFFB_TEST_RUN_DIR=${run_dir}"
            ${ARGN}
            -P "${CMAKE_CURRENT_SOURCE_DIR}/run_solver.cmake"
    )
    configure_ffb_test(${name} ${images} ${ix} ${iy} "${labels}")
    if(NOT fixture STREQUAL "")
        set_tests_properties(${name} PROPERTIES FIXTURES_SETUP "${fixture}")
    endif()
endfunction()

function(add_ffb_test_program name source_file)
    set(target_name FFB_TEST_${name})
    add_executable(${target_name} EXCLUDE_FROM_ALL support/assertions.f90 ${source_file})
    configure_ffb_target(${target_name} xyd single)
    target_link_libraries(${target_name} PRIVATE FFB_TEST_serial_core)
    add_dependencies(ffb_tests ${target_name})
endfunction()

add_custom_target(ffb_tests)
set(FFB_TEST_MODES shear_wave couette_flow poiseuille_flow sliding_lid)
set(FFB_SIM_MODE 0)
foreach(mode IN LISTS FFB_TEST_MODES)
    math(EXPR FFB_SIM_MODE "${FFB_SIM_MODE} + 1")
    add_ffb_test_solver(FFB_TEST_${mode} ${FFB_SIM_MODE})
endforeach()

# validation helpers run in one process so diagnostics do not depend on MPI shutdown
add_ffb_test_core(FFB_TEST_serial_core 1 single)

foreach(name collision boundaries domain)
    add_ffb_test_program(${name} unit/test_${name}.f90)
    add_test(NAME unit.${name} COMMAND FFB_TEST_${name})
    configure_ffb_test(unit.${name} 1 1 1 "unit")
endforeach()
add_ffb_test_program(fields integration/test_fields.f90)

set(FFB_SIM_MODE 0)
foreach(mode IN LISTS FFB_TEST_MODES)
    math(EXPR FFB_SIM_MODE "${FFB_SIM_MODE} + 1")
    set(serial_dir "${FFB_TEST_ROOT}/${mode}/serial")
    add_ffb_run(run.${mode}.serial FFB_TEST_${mode} "${serial_dir}" 1 1 1 "${mode}_serial" "smoke")

    set(check_name smoke.${mode})
    set(check_labels "smoke")
    if(mode STREQUAL "shear_wave")
        set(check_name numerical.shear_wave)
        set(check_labels "smoke;numerical")
    endif()
    add_test(NAME ${check_name} COMMAND FFB_TEST_fields ${FFB_SIM_MODE} 1 1)
    configure_ffb_test(${check_name} 1 1 1 "${check_labels}")
    set_tests_properties(${check_name} PROPERTIES
        WORKING_DIRECTORY "${serial_dir}"
        FIXTURES_REQUIRED "${mode}_serial"
    )

    # each comparison requires both the serial and distributed solver output
    set(decompositions x y)
    if(FFB_TEST_EXTENDED)
        list(APPEND decompositions xy)
    endif()
    foreach(decomposition IN LISTS decompositions)
        set(ix 1)
        set(iy 1)
        if(decomposition MATCHES "x")
            set(ix 2)
        endif()
        if(decomposition MATCHES "y")
            set(iy 2)
        endif()
        math(EXPR images "${ix} * ${iy}")
        set(run_dir "${FFB_TEST_ROOT}/${mode}/${decomposition}")
        set(labels "distributed")
        if(decomposition STREQUAL "xy")
            set(labels "distributed;extended")
        endif()
        add_ffb_run(run.${mode}.${decomposition} FFB_TEST_${mode} "${run_dir}"
            ${images} ${ix} ${iy} "${mode}_${decomposition}" "${labels}")
        set(check_name distributed.${mode}.${decomposition})
        add_test(NAME ${check_name}
            COMMAND FFB_TEST_fields ${FFB_SIM_MODE} ${ix} ${iy} "${serial_dir}")
        configure_ffb_test(${check_name} 1 1 1 "${labels}")
        set_tests_properties(${check_name} PROPERTIES
            WORKING_DIRECTORY "${run_dir}"
            FIXTURES_REQUIRED "${mode}_serial;${mode}_${decomposition}"
        )
    endforeach()
endforeach()

add_ffb_run(configuration.positive FFB_TEST_domain "${FFB_TEST_ROOT}/invalid/positive" 1 0 1 "" "configuration"
    "-DFFB_EXPECTED_ERROR=I_X and I_Y must be positive")
add_ffb_run(configuration.product FFB_TEST_domain "${FFB_TEST_ROOT}/invalid/product" 1 2 1 "" "configuration"
    "-DFFB_EXPECTED_ERROR=I_X * I_Y must match coarray images")

# the solver must reject a grid that does not match the real image count
add_ffb_run(configuration.solver FFB_TEST_shear_wave "${FFB_TEST_ROOT}/invalid/solver" 2 1 1 "" "configuration;distributed"
    "-DFFB_EXPECTED_ERROR=I_X * I_Y must match coarray images")
