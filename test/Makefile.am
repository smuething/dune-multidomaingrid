testdir = $(includedir)/test

include $(top_srcdir)/am/global-rules

EXTRA_DIST = CMakeLists.txt

check_PROGRAMS = test-wrapper testpartitioning testlargedomainnumbers testadaptation profilerepartitioning iterateallinterfaces testparallel

test_wrapper_SOURCES = test-wrapper.cc
test_wrapper_CXXFLAGS = $(AM_CXXFLAGS) $(ALL_PKG_CPPFLAGS) $(BOOST_CPPFLAGS)
test_wrapper_LDFLAGS = $(ALL_PKG_LDFLAGS)
test_wrapper_LIBS = $(ALL_PKG_LIBS)

testpartitioning_SOURCES = testpartitioning.cc
testpartitioning_CXXFLAGS = $(AM_CXXFLAGS) $(ALL_PKG_CPPFLAGS) $(BOOST_CPPFLAGS)
testpartitioning_LDFLAGS = $(ALL_PKG_LDFLAGS)
testpartitioning_LIBS = $(ALL_PKG_LIBS)

testlargedomainnumbers_SOURCES = testlargedomainnumbers.cc
testlargedomainnumbers_CXXFLAGS = $(AM_CXXFLAGS) -DGRIDDIM=2 -DYASPGRID $(ALL_PKG_CPPFLAGS) $(BOOST_CPPFLAGS)
testlargedomainnumbers_LDFLAGS = $(ALL_PKG_LDFLAGS)
testlargedomainnumbers_LIBS = $(ALL_PKG_LIBS)

testadaptation_SOURCES = testadaptation.cc
testadaptation_CXXFLAGS = $(AM_CXXFLAGS) $(ALL_PKG_CPPFLAGS) $(BOOST_CPPFLAGS)
testadaptation_LDFLAGS = $(ALL_PKG_LDFLAGS) -lalugrid
testadaptation_LIBS = $(ALL_PKG_LIBS) -lalugrid

profilerepartitioning_SOURCES = profilerepartitioning.cc
profilerepartitioning_CXXFLAGS = $(AM_CXXFLAGS) -DGRIDDIM=2 -DYASPGRID $(ALL_PKG_CPPFLAGS) $(BOOST_CPPFLAGS)
profilerepartitioning_LDFLAGS = $(ALL_PKG_LDFLAGS)
profilerepartitioning_LIBS = $(ALL_PKG_LIBS)

iterateallinterfaces_SOURCES = iterateallinterfaces.cc

testparallel_SOURCES = testparallel.cc
testparallel_CXXFLAGS = $(AM_CXXFLAGS) $(ALL_PKG_CPPFLAGS) $(BOOST_CPPFLAGS)
testparallel_LDFLAGS = $(ALL_PKG_LDFLAGS)
testparallel_LIBS = $(ALL_PKG_LIBS)

check_PROGRAMS += multidomain-leveliterator-bug
multidomain_leveliterator_bug_SOURCES = multidomain-leveliterator-bug.cc

check_PROGRAMS += testintersectionconversion
testintersectionconversion_SOURCES = testintersectionconversion.cc

check_PROGRAMS += testintersectiongeometrytypes
testintersectiongeometrytypes_SOURCES = testintersectiongeometrytypes.cc
testintersectiongeometrytypes_CXXFLAGS = $(AM_CXXFLAGS) $(ALL_PKG_CPPFLAGS) $(BOOST_CPPFLAGS)
testintersectiongeometrytypes_LDFLAGS = $(AM_LDFLAGS) $(ALL_PKG_LDFLAGS)
testintersectiongeometrytypes_LIBS = $(AM_LIBS) $(ALL_PKG_LIBS)
