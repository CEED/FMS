# fms_add_cpp_test(
#   TEST      - The test name & name of .cpp file
#   [SOURCES] - Any addition source files
#   [ARGS]    - Command line arguments to invoke the test with
#   [DEPENDS] - Other CMake targets that this test needs to link to
#   [PROPERTIES]  - Test properties (https://cmake.org/cmake/help/v3.10/manual/cmake-properties.7.html#test-properties)
# )
# Defines a target and a test, both have the name specified in the TEST argument.
function(fms_add_cpp_test)
  set(options)
  set(singleValueArgs TEST)
  set(multiValueArgs ARGS SOURCES DEPENDS PROPERTIES)
  cmake_parse_arguments(
    arg
    "${options}"
    "${singleValueArgs}"
    "${multiValueArgs}" ${ARGN})

  # Create the executable and link it to anything the test depends on
  add_executable(${arg_TEST} ${arg_TEST}.cpp ${arg_SOURCES})
  target_link_libraries(${arg_TEST} PRIVATE fms_testing_common fms gtest_main ${arg_DEPENDS})

  # Add the actual test to "make test" with optional command line arguments
  add_test(
    NAME ${arg_TEST}
    COMMAND ${arg_TEST} ${arg_ARGS})

  # Add any test properties that were given to the function
  if("${arg_PROPERTIES}")
    set_tests_properties(${arg_TEST}
      PROPERTIES ${arg_PROPERTIES})
  endif()

  # Make sure this test gets the Fms label, have todo it like this incase someone gave us other labels too
  get_test_property(${arg_TEST} LABELS ${arg_TEST}_LABELS)
  if("${${arg_TEST}_LABELS}" STREQUAL "NOTFOUND")
    set_tests_properties(${arg_TEST}
      PROPERTIES
        LABELS Fms)
  else()
    set_tests_properties(${arg_TEST}
      PROPERTIES
        LABELS "Fms;${${arg_TEST}_LABELS}")
  endif()
endfunction()
# ------------------------------------------------------------------------

# Add FmsBaselines
add_library(fms_testing_common STATIC fmstesting.hpp fmstesting.cpp)
target_link_libraries(fms_testing_common PUBLIC fms)
target_include_directories(fms_testing_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

fms_add_cpp_test(
  TEST t_fms
  PROPERTIES
    TIMEOUT 1
)

fms_add_cpp_test(
  TEST t_fmsio
)