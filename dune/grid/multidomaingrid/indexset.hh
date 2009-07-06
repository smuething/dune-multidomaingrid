#ifndef DUNE_MULTIDOMAINGRID_INDEXSET_HH
#define DUNE_MULTIDOMAINGRID_INDEXSET_HH



class IndexSet {

  static const std::size_t maxSubDomains = SubDomainSet::maxSize;

  struct MapEntry {
    SubDomainSet domains;
    IndexType index;
  };

  void update() {
    HostIndexSet his = ...;
    _indexMap.resize(his.size(0));
    _multiIndexMap.clear();
    const std::size_t numGeomTypes = his.geomTypes().size();
    std::array<std::map<GeometryType,IndexType>,maxSubDomains> indices();
    HostEntityIterator end = ...;
    for (HostEntityIterator it  = ...; it != end; ++it) {
      IndexType hostIndex = his.index(*it);
      MapEntry& e = _indexMap[hostIndex];
      switch (domains.state()) {
      case SubDomainSet::simpleSet:
	e.index = indices[*set.begin()]++;
	break;
      case SubDomainSet::multipleSet:
	_multiIndexMap.push_back(MultiIndexContainer());
	MultiIndexContainer& mic = _multiIndexMap.back();
	e.index = _multiIndexMap.size();
	for (SubDomainSet::Iterator it = e.domains.begin(); it != e.domains.end(); ++it) {
	  mic[*it] = indices[*it]++;
	}
      }
    }
  }

private:
  std::vector<MapEntry> _indexMap;
  std::vector<MultiIndexContainer> _multiIndexMap;

};


#endif // DUNE_MULTIDOMAINGRID_INDEXSET_HH
