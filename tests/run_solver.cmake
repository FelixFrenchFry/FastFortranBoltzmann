cmake_path(GET FFB_TEST_EXE FILENAME test_program)

# only remove outputs inside the build's dedicated test run directory
file(REAL_PATH "${FFB_TEST_ROOT}" test_root)
file(REAL_PATH "${FFB_TEST_RUN_DIR}" run_dir)
cmake_path(IS_PREFIX test_root "${run_dir}" NORMALIZE inside_test_root)
if(NOT inside_test_root OR run_dir STREQUAL test_root)
    message(FATAL_ERROR "Test run directory must be below ${test_root}")
endif()
file(REMOVE_RECURSE "${run_dir}/output/test")

execute_process(
    COMMAND "${FFB_TEST_EXE}"
    WORKING_DIRECTORY "${run_dir}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE error
    TIMEOUT 30
)
file(WRITE "${run_dir}/solver.log" "exit status: ${result}\n${output}\n${error}")

# a timeout or launcher signal must not count as an expected domain rejection
if(NOT result MATCHES "^[0-9]+$")
    message(FATAL_ERROR "${test_program} did not finish: ${result}\n${output}\n${error}")
endif()
if("${output}\n${error}" MATCHES "forrtl: warning [(]406[)]")
    message(FATAL_ERROR "${test_program} emitted warning 406 with exit status ${result}\n${output}\n${error}")
endif()
if(DEFINED FFB_EXPECTED_ERROR)
    string(FIND "${output}\n${error}" "${FFB_EXPECTED_ERROR}" diagnostic_position)
    if(result EQUAL 0 OR diagnostic_position EQUAL -1)
        message(FATAL_ERROR "${test_program} exited with ${result}, expected a nonzero exit and '${FFB_EXPECTED_ERROR}'\n${output}\n${error}")
    endif()
elseif(NOT result EQUAL 0)
    message(FATAL_ERROR "${test_program} exited with ${result}\n${output}\n${error}")
endif()
