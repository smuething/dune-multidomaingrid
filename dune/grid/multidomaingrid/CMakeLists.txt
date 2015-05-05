add_subdirectory(subdomaingrid)

#install headers
install(FILES
  allsubdomaininterfacesiterator.hh
  arraybasedset.hh
  entity.hh
  entitypointer.hh
  geometry.hh
  gridview.hh
  hierarchiciterator.hh
  hostgridaccessor.hh
  idsets.hh
  indexsets.hh
  intersection.hh
  intersectioniterator.hh
  iterator.hh
  localgeometry.hh
  mdgridtraits.hh
  multidomaingrid.hh
  multidomainmcmgmapper.hh
  singlevalueset.hh
  subdomaininterfaceiterator.hh
  subdomainset.hh
  subdomaintosubdomaininterfaceiterator.hh
  utility.hh
  DESTINATION include/dune/grid/multidomaingrid)
