#ifndef DUNE_MULTIDOMAINGRID_INDEXSET_HH
#define DUNE_MULTIDOMAINGRID_INDEXSET_HH

#include <unordered_map>
#include <vector>
#include <array>

namespace Dune {

namespace multidomaingrid {

struct GeometryTypeHash {

  std::size_t operator()(GeometryType gt) const {
    std::size_t hash = gt.dim() * 509;
    return gt.dim() < 2 ? hash : hash + static_cast<std::size_t>(gt.basicType());
  }

};

template<typename HostGridView, typename SubDomainSetType>
class IndexSet {

public:

  typedef SubDomainSetType SubDomainSet;

  typedef typename HostGridView::IndexSet HostIndexSet;
  typedef typename HostGridView::template Codim<0>::Iterator HostEntityIterator;
  typedef typename HostGridView::template Codim<0>::Entity HostEntity;
  typedef typename HostIndexSet::IndexType IndexType;
  typedef typename SubDomainSet::DomainType DomainType;
  typedef typename HostGridView::Grid::ctype ctype;
  static const std::size_t maxSubDomains = SubDomainSet::maxSize;
  static const std::size_t dimension = HostGridView::dimension;

  struct MapEntry {
    SubDomainSet domains;
    IndexType index;
  };

  typedef std::unordered_map<GeometryType,std::vector<MapEntry>,GeometryTypeHash> CellMap;
  typedef std::unordered_map<GeometryType,std::vector<IndexType>,GeometryTypeHash> SubIndexMap;
  typedef std::unordered_map<GeometryType,std::array<IndexType,maxSubDomains>,GeometryTypeHash> SizeMap;
  typedef std::array<IndexType,maxSubDomains> MultiIndexContainer;

  void update(bool full) {
    const HostIndexSet& his = _gridView.indexSet();
    const std::vector<GeometryType>& hgtList = his.geomTypes(0);
    if (full) {
      _indexMap.swap(CellMap());
      _indexMap.rehash(hgtList.size());
      _sizeMap.swap(SizeMap());
      _sizeMap.rehash(hgtList.size());
    }
    for (std::vector<GeometryType>::const_iterator it = hgtList.begin(); it != hgtList.end(); ++it) {
      _indexMap[*it].resize(his.size(*it));
      _sizeMap[*it].assign(0);
    }
    _multiIndexMap.clear();
    HostEntityIterator end = _gridView.template end<0>();
    for (HostEntityIterator it  = _gridView.template begin<0>(); it != end; ++it) {
      const HostEntity& he = *it;
      const GeometryType hgt = he.type();
      IndexType hostIndex = his.index(he);
      MapEntry& e = _indexMap[hgt][hostIndex];
      switch (e.domains.state()) {
      case SubDomainSet::simpleSet:
	e.index = _sizeMap[hgt][*e.domains.begin()]++;
	break;
      case SubDomainSet::multipleSet:
	e.index = _multiCellMap.size();
	_multiCellMap.push_back(MultiIndexContainer());
	MultiIndexContainer& mic = _multiIndexMap.back();
	for (typename SubDomainSet::Iterator sdit = e.domains.begin(); sdit != e.domains.end(); ++sdit) {
	  mic[*sdit] = _sizeMap[hgt][*sdit]++;
	}
      }
    }
  }

  template<int codim>
  void updateSubIndices(const HostEntity& e, const GeometryType gt, const SubDomainSet& domains, const HostIndexSet& his) {
    int count = GenericReferenceElements<ctype,dim-codim>::general(gt).size(codim);
    for (int i = 0; i < count; ++i) {
      IndexType hostIndex = his.subIndex(e,i,codim);
    }
  }

  SubDomainSet& subDomainSet(const HostEntity& e) {
    return _indexMap[e.type()][_gridView.indexSet().index(e)].domains;
  }

  IndexType indexForSubDomain(DomainType subDomain, const HostEntity& e) {
    GeometryType gt = e.type();
    IndexType hostIndex = _gridView.indexSet().index(e);
    const MapEntry& me = _indexMap[gt][hostIndex];
    assert(me.domains.contains(subDomain));
    if (me.domains.simple()) {
      return me.index;
    } else {
      return _multiIndexMap[me.index][subDomain];
    }
  }

  explicit IndexSet(HostGridView gv) :
    _gridView(gv)
  {}

private:
  CellMap _indexMap;
  std::array<SizeMap,dimension+1> _sizeMap;
  std::vector<MultiIndexContainer> _multiIndexMap;
  HostGridView _gridView;

};

} // namespace multidomaingrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_INDEXSET_HH
