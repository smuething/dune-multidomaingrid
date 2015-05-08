# This target will be used to build all tests
add_custom_target(build_tests)

# Function to add a new test
# This function
# - creates the executable target and registers it as a dependecy
#   for build_tests
# - optionally sets compiler flags
# - registers the test
function(add_mdgrid_test)
  include(CMakeParseArguments)
  set(OPTIONS)
  set(SINGLEARGS NAME)
  set(MULTIARGS SOURCES COMPILE_DEFINITIONS COMPILE_OPTIONS COMMAND)
  cmake_parse_arguments(MDTEST "${OPTIONS}" "${SINGLEARGS}" "${MULTIARGS}" ${ARGN})

  if(MDTEST_UNPARSED_ARGUMENTS)
    message(WARNING "Unrecognized arguments ('${MDTEST_UNPARSED_ARGUMENTS}') for add_mdgrid_test()!")
  endif()

  # by default, a test is built from a file with the same name and a ".cc" suffix
  if("${MDTEST_SOURCES}" STREQUAL "")
    set(MDTEST_SOURCES "${MDTEST_NAME}.cc")
  endif()

  add_executable(
    "${MDTEST_NAME}"
    EXCLUDE_FROM_ALL
    ${MDTEST_SOURCES}
    )

  add_dependencies(build_tests "${MDTEST_NAME}")

  target_compile_definitions(
    "${MDTEST_NAME}"
    PUBLIC
    ${MDTEST_COMPILE_DEFINITIONS}
    )

  target_compile_options(
    "${MDTEST_NAME}"
    PUBLIC
    ${MDTEST_COMPILE_OPTIONS}
    )

  if("${MDTEST_COMMAND}" STREQUAL "")
    set(MDTEST_COMMAND "${MDTEST_NAME}")
  endif()

  # by default, a test is run by simply invoking the built executable
  add_test(
    NAME "${MDTEST_NAME}"
    COMMAND ${MDTEST_COMMAND}
    )
endfunction()

add_executable(iterateallinterfaces EXCLUDE_FROM_ALL iterateallinterfaces.cc)
add_test(NAME iterateallinterfaces
  COMMAND iterateallinterfaces "3")
add_executable(multidomain-leveliterator-bug EXCLUDE_FROM_ALL multidomain-leveliterator-bug.cc)
add_executable(profilerepartitioning EXCLUDE_FROM_ALL profilerepartitioning.cc)
set_target_properties(profilerepartitioning
  PROPERTIES COMPILE_FLAGS "-DGRIDDIM=2 -DYASPGRID")
add_executable(testadaptation EXCLUDE_FROM_ALL testadaptation.cc)
add_executable(testintersectionconversion EXCLUDE_FROM_ALL testintersectionconversion.cc)
add_executable(testintersectiongeometrytypes EXCLUDE_FROM_ALL testintersectiongeometrytypes.cc)
add_executable(testlargedomainnumbers EXCLUDE_FROM_ALL testlargedomainnumbers.cc)
set_target_properties(testlargedomainnumbers
  PROPERTIES COMPILE_FLAGS "-DGRIDDIM=2 -DYASPGRID")
add_executable(testparallel EXCLUDE_FROM_ALL testparallel.cc)
add_test(NAME testparallel
  COMMAND testparallel "5" "1")
add_executable(testpartitioning EXCLUDE_FROM_ALL testpartitioning.cc)
add_executable(test-wrapper EXCLUDE_FROM_ALL test-wrapper.cc)
