dunedir = $(includedir)/dune

include $(top_srcdir)/am/global-rules

noinst_PROGRAMS = testsubdomainset testindexset

testsubdomainset_SOURCES = testsubdomainset.cc
testindexset_SOURCES = testindexset.cc
testindexset_CXXFLAGS = $(AM_CXXFLAGS) -DGRIDDIM=2 -DYASPGRID $(ALL_PKG_CPPFLAGS) $(BOOST_CPPFLAGS)
testindexset_LDFLAGS = $(ALL_PKG_LDFLAGS)
testindexset_LIBS = $(ALL_PKG_LIBS)