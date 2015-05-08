set(tests
  multidomain-leveliterator-bug
  testadaptation
  testintersectionconversion
  testintersectiongeometrytypes
  testlargedomainnumbers
  testpartitioning
  test-wrapper
)

set(tests_with_arguments
  iterateallinterfaces
  testparallel
)

set(all_tests "${tests};${tests_with_arguments}")

add_executable(iterateallinterfaces EXCLUDE_FROM_ALL iterateallinterfaces.cc)
add_test(NAME iterateallinterfaces
  COMMAND iterateallinterfaces "3")
add_executable(multidomain-leveliterator-bug EXCLUDE_FROM_ALL multidomain-leveliterator-bug.cc)
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

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${all_tests})

# install the header as done for the auto-tools
install(FILES output.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/test)

foreach(_test ${tests})
  add_test(${_test} ${_test})
endforeach()

foreach(_test ${all_tests})
  target_link_libraries(${_test} ${DUNE_LIBS})
endforeach()
