dunedir = $(includedir)/dune

include $(top_srcdir)/am/global-rules

#noinst_PROGRAMS = testindexset testwrapper

#testsubdomainset_SOURCES = testsubdomainset.cc

#testindexset_SOURCES = testindexset.cc
#testindexset_CXXFLAGS = $(AM_CXXFLAGS) -DGRIDDIM=2 -DYASPGRID $(ALL_PKG_CPPFLAGS) $(BOOST_CPPFLAGS)
#testindexset_LDFLAGS = $(ALL_PKG_LDFLAGS)
#testindexset_LIBS = $(ALL_PKG_LIBS)

#testwrapper_SOURCES = testwrapper.cc
#testwrapper_CXXFLAGS = $(AM_CXXFLAGS) -DGRIDDIM=2 -DYASPGRID $(ALL_PKG_CPPFLAGS) $(BOOST_CPPFLAGS)
#testwrapper_LDFLAGS = $(ALL_PKG_LDFLAGS)
#testwrapper_LIBS = $(ALL_PKG_LIBS)

check_PROGRAMS = test-wrapper testpartitioning testlargedomainnumbers

test_wrapper_SOURCES = test-wrapper.cc
test_wrapper_CXXFLAGS = $(AM_CXXFLAGS) -DGRIDDIM=2 -DYASPGRID $(ALL_PKG_CPPFLAGS) $(BOOST_CPPFLAGS)
test_wrapper_LDFLAGS = $(ALL_PKG_LDFLAGS)
test_wrapper_LIBS = $(ALL_PKG_LIBS)

testpartitioning_SOURCES = testpartitioning.cc
testpartitioning_CXXFLAGS = $(AM_CXXFLAGS) -DGRIDDIM=2 -DYASPGRID $(ALL_PKG_CPPFLAGS) $(BOOST_CPPFLAGS)
testpartitioning_LDFLAGS = $(ALL_PKG_LDFLAGS)
testpartitioning_LIBS = $(ALL_PKG_LIBS)

testlargedomainnumbers_SOURCES = testlargedomainnumbers.cc
testlargedomainnumbers_CXXFLAGS = $(AM_CXXFLAGS) -DGRIDDIM=2 -DYASPGRID $(ALL_PKG_CPPFLAGS) $(BOOST_CPPFLAGS)
testlargedomainnumbers_LDFLAGS = $(ALL_PKG_LDFLAGS)
testlargedomainnumbers_LIBS = $(ALL_PKG_LIBS)
