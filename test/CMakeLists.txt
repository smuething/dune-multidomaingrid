set(TESTS
  iterateallinterfaces
  multidomain-leveliterator-bug
  profilerepartitioning
  testadaptation
  testintersectionconversion
  testintersectiongeometrytypes
  testlargedomainnumbers
  testparallel
  testpartitioning
  test_wrapper)

add_executable(iterateallinterfaces EXCLUDE_FROM_ALL iterateallinterfaces.cc)
add_executable(multidomain-leveliterator-bug EXCLUDE_FROM_ALL multidomain-leveliterator-bug.cc)
add_executable(profilerepartitioning EXCLUDE_FROM_ALL profilerepartitioning.cc)
add_executable(testadaptation EXCLUDE_FROM_ALL testadaptation.cc)
add_executable(testintersectionconversion EXCLUDE_FROM_ALL testintersectionconversion.cc)
add_executable(testintersectiongeometrytypes EXCLUDE_FROM_ALL testintersectiongeometrytypes.cc)
add_executable(testlargedomainnumbers EXCLUDE_FROM_ALL testlargedomainnumbers.cc)
add_executable(testparallel EXCLUDE_FROM_ALL testparallel.cc)
add_executable(testpartitioning EXCLUDE_FROM_ALL testpartitioning.cc)
add_executable(test_wrapper EXCLUDE_FROM_ALL test-wrapper.cc)

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS})

# install the header as done for the auto-tools
install(FILES output.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/test)

foreach(_test ${TESTS})
  add_test(${_test} ${_test})
  target_link_libraries(${_test} dunecommon dunegeometry dunegrid)
endforeach(_test ${TESTS})
