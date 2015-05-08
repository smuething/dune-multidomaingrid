testdir = $(includedir)/test

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt

check_PROGRAMS = \
  test-wrapper-multidomaingrid \
  test-wrapper-subdomaingrid \
  testpartitioning \
  testlargedomainnumbers \
  testadaptation \
  profilerepartitioning \
  iterateallinterfaces \
  testparallel

test_wrapper_multidomaingrid_SOURCES = test-wrapper.cc
test_wrapper_multidomaingrid_CXXFLAGS = \
  $(AM_CXXFLAGS) \
  $(ALL_PKG_CPPFLAGS) \
  -DCHECK_DIMENSION=2 \
  -DCHECK_GRIDSIZE=16 \
  -DCHECK_MULTIDOMAINGRID
test_wrapper_multidomaingrid_LDFLAGS = $(ALL_PKG_LDFLAGS)
test_wrapper_multidomaingrid_LIBS = $(ALL_PKG_LIBS)

test_wrapper_subdomaingrid_SOURCES = test-wrapper.cc
test_wrapper_subdomaingrid_CXXFLAGS = \
  $(AM_CXXFLAGS) \
  $(ALL_PKG_CPPFLAGS) \
  -DCHECK_DIMENSION=2 \
  -DCHECK_GRIDSIZE=16 \
  -DCHECK_SUBDOMAINGRID
test_wrapper_subdomaingrid_LDFLAGS = $(ALL_PKG_LDFLAGS)
test_wrapper_subdomaingrid_LIBS = $(ALL_PKG_LIBS)

testpartitioning_SOURCES = testpartitioning.cc
testpartitioning_CXXFLAGS = $(AM_CXXFLAGS) $(ALL_PKG_CPPFLAGS)
testpartitioning_LDFLAGS = $(ALL_PKG_LDFLAGS)
testpartitioning_LIBS = $(ALL_PKG_LIBS)

testlargedomainnumbers_SOURCES = testlargedomainnumbers.cc
testlargedomainnumbers_CXXFLAGS = $(AM_CXXFLAGS) -DGRIDDIM=2 -DYASPGRID $(ALL_PKG_CPPFLAGS)
testlargedomainnumbers_LDFLAGS = $(ALL_PKG_LDFLAGS)
testlargedomainnumbers_LIBS = $(ALL_PKG_LIBS)

testadaptation_SOURCES = testadaptation.cc
testadaptation_CXXFLAGS = $(AM_CXXFLAGS) $(ALL_PKG_CPPFLAGS)
testadaptation_LDFLAGS = $(ALL_PKG_LDFLAGS)
testadaptation_LIBS = $(ALL_PKG_LIBS)

profilerepartitioning_SOURCES = profilerepartitioning.cc
profilerepartitioning_CXXFLAGS = $(AM_CXXFLAGS) -DGRIDDIM=2 -DYASPGRID $(ALL_PKG_CPPFLAGS)
profilerepartitioning_LDFLAGS = $(ALL_PKG_LDFLAGS)
profilerepartitioning_LIBS = $(ALL_PKG_LIBS)

iterateallinterfaces_SOURCES = iterateallinterfaces.cc

testparallel_SOURCES = testparallel.cc
testparallel_CXXFLAGS = $(AM_CXXFLAGS) $(ALL_PKG_CPPFLAGS)
testparallel_LDFLAGS = $(ALL_PKG_LDFLAGS)
testparallel_LIBS = $(ALL_PKG_LIBS)

check_PROGRAMS += multidomain-leveliterator-bug
multidomain_leveliterator_bug_SOURCES = multidomain-leveliterator-bug.cc

check_PROGRAMS += testintersectionconversion
testintersectionconversion_SOURCES = testintersectionconversion.cc

check_PROGRAMS += testintersectiongeometrytypes
testintersectiongeometrytypes_SOURCES = testintersectiongeometrytypes.cc
testintersectiongeometrytypes_CXXFLAGS = $(AM_CXXFLAGS) $(ALL_PKG_CPPFLAGS)
testintersectiongeometrytypes_LDFLAGS = $(AM_LDFLAGS) $(ALL_PKG_LDFLAGS)
testintersectiongeometrytypes_LIBS = $(AM_LIBS) $(ALL_PKG_LIBS)
