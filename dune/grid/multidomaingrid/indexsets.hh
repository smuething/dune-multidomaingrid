#ifndef DUNE_MULTIDOMAINGRID_INDEXSETS_HH
#define DUNE_MULTIDOMAINGRID_INDEXSETS_HH

#include <unordered_map>
#include <vector>
#include <array>
#include <algorithm>
#include <type_traits>
#include <tuple>
#include <boost/scoped_ptr.hpp>
#include <boost/bind.hpp>
#include <dune/grid/multidomaingrid/utility.hh>

namespace Dune {

namespace mdgrid {

template<typename HostGrid>
class MultiDomainGrid;


template<typename GridImp, typename HostGridViewType>
class IndexSetWrapper :
    public Dune::IndexSet<GridImp,IndexSetWrapper<GridImp,HostGridViewType>,
                          typename HostGridViewType::IndexSet::IndexType>
{

  template<typename, typename>
  friend class IndexSetWrapper;

  template<typename HostGrid>
  friend class MultiDomainGrid;

  typedef IndexSetWrapper<GridImp,HostGridViewType> ThisType;

  typedef typename remove_const<GridImp>::type::HostGridType HostGrid;
  typedef HostGridViewType HostGridView;
  typedef typename HostGridView::IndexSet HostIndexSet;
  typedef typename remove_const<GridImp>::type::ctype ctype;

public:

  typedef typename remove_const<GridImp>::type::SubDomainSet SubDomainSet;
  typedef typename SubDomainSet::DomainType SubDomainType;
  typedef typename HostIndexSet::IndexType IndexType;
  static const int dimension = remove_const<GridImp>::type::dimension;
  typedef typename SubDomainSet::DomainType DomainType;
  static const std::size_t maxSubDomains = SubDomainSet::maxSize;

private:

  typedef typename HostGridView::template Codim<0>::Iterator HostEntityIterator;
  typedef typename HostGridView::template Codim<0>::Entity HostEntity;
  typedef typename HostGridView::template Codim<0>::EntityPointer HostEntityPointer;
  typedef typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity Codim0Entity;

  struct MapEntry {
    SubDomainSet domains;
    IndexType index;
  };

  typedef std::array<IndexType,maxSubDomains> SizeContainer;
  typedef std::unordered_map<GeometryType,std::vector<MapEntry>,util::GeometryTypeHash> IndexMap;
  typedef std::unordered_map<GeometryType,SizeContainer,util::GeometryTypeHash> SizeMap;
  typedef std::array<IndexType,maxSubDomains> MultiIndexContainer;

  typedef std::vector<boost::shared_ptr<IndexSetWrapper<GridImp, typename HostGridView::Grid::LevelGridView> > > LevelIndexSets;

public:

  template<int codim>
  IndexType index(const typename remove_const<GridImp>::type::Traits::template Codim<codim>::Entity& e) const {
    return _hostGridView.indexSet().index(_grid.hostEntity(e));
  }

  template<typename Entity>
  IndexType index(const Entity& e) const {
    return _hostGridView.indexSet().index(_grid.hostEntity(e));
  }

  template<int codim>
  IndexType subIndex(const Codim0Entity& e, int i) const {
    return _hostGridView.indexSet().subIndex(_grid.hostEntity(e),i,codim);
  }

  IndexType subIndex(const Codim0Entity& e, int i, unsigned int codim) const {
    return _hostGridView.indexSet().subIndex(_grid.hostEntity(e),i,codim);
  }

  const std::vector<GeometryType>& geomTypes(int codim) const {
    return _hostGridView.indexSet().geomTypes(codim);
  }

  IndexType size(GeometryType type) const {
    return _hostGridView.indexSet().size(type);
  }

  IndexType size(int codim) const {
    return _hostGridView.indexSet().size(codim);
  }

  template<typename EntityType>
  bool contains(const EntityType& e) const {
    return _hostGridView.indexSet().contains(_grid.hostEntity(e));
  }

  template<typename EntityType>
  const SubDomainSet& subDomains(const EntityType& e) const {
    return subDomains<EntityType::codimension>(e);
  }

  template<int cc>
  const SubDomainSet& subDomains(const typename remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& e) const {
    return indexMap(cc).find(e.type())->second[_hostGridView.indexSet().index(_grid.hostEntity(e))].domains;
  }

  template<class EntityType>
  IndexType index(DomainType subDomain, const EntityType& e) const {
    return index<EntityType::codimension>(subDomain,e);
  }

  template<int cc>
  IndexType index(DomainType subDomain, const typename remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& e) const {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    const MapEntry& me = indexMap(cc).find(gt)->second[hostIndex];
    assert(me.domains.contains(subDomain));
    if (me.domains.simple()) {
      return me.index;
    } else {
      return _multiIndexMap[me.index][subDomain];
    }
  }

  template<int cc>
  IndexType indexForSubDomain(DomainType subDomain, const typename remove_const<GridImp>::type::HostGridType::Traits::template Codim<cc>::Entity& he) const {
    const GeometryType gt = he.type();
    const IndexType hostIndex = _hostGridView.indexSet().index(he);
    const MapEntry& me = indexMap(cc).find(gt)->second[hostIndex];
    assert(me.domains.contains(subDomain));
    if (me.domains.simple()) {
      return me.index;
    } else {
      return _multiIndexMap[me.index][subDomain];
    }
  }

  IndexType subIndexForSubDomain(DomainType subDomain, const typename remove_const<GridImp>::type::HostGridType::Traits::template Codim<0>::Entity& he, int i, int codim) const {
    const GeometryType gt = GenericReferenceElements<ctype,dimension>::general(he.type()).type(i,codim);
    const IndexType hostIndex = _hostGridView.indexSet().subIndex(he,i,codim);
    const MapEntry& me = indexMap(codim).find(gt)->second[hostIndex];
    assert(me.domains.contains(subDomain));
    if (me.domains.simple()) {
      return me.index;
    } else {
      return _multiIndexMap[me.index][subDomain];
    }
  }

  const std::vector<GeometryType>& geomTypesForSubDomain(DomainType subDomain, int codim) const {
    return geomTypes(codim);
  }

  IndexType sizeForSubDomain(DomainType subDomain, GeometryType type) const {
    return sizeMap(dimension-type.dim()).find(type)->second[subDomain];
  }

  IndexType sizeForSubDomain(DomainType subDomain, int codim) const {
    return _codimSizes[codim][subDomain];
  }

  template<typename EntityType>
  bool containsForSubDomain(DomainType subDomain, const EntityType& he) const {
    const GeometryType gt = he.type();
    const IndexType hostIndex = _hostGridView.indexSet().index(he);
    const MapEntry& me = indexMap(EntityType::codimension).find(gt)->second[hostIndex];
    return me.domains.contains(subDomain);
  }

private:

  const GridImp& _grid;
  HostGridView _hostGridView;
  std::array<boost::scoped_ptr<IndexMap>,dimension+1> _indexMap;
  std::array<boost::scoped_ptr<SizeMap>,dimension+1> _sizeMap;
  std::array<std::array<IndexType,maxSubDomains>,dimension+1> _codimSizes;
  std::vector<MultiIndexContainer> _multiIndexMap;

  void swap(ThisType& rhs) {
    assert(&_grid == &rhs._grid);
    std::swap(_indexMap,rhs._indexMap);
    std::swap(_sizeMap,rhs._sizeMap);
    std::swap(_codimSizes,rhs._codimSizes);
    std::swap(_multiIndexMap,rhs._multiIndexMap);
  }

  void addToSubDomain(SubDomainType subDomain, const Codim0Entity& e) {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    indexMap(0)[gt][hostIndex].domains.add(subDomain);
  }

  void removeFromSubDomain(SubDomainType subDomain, const Codim0Entity& e) {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    indexMap(0)[gt][hostIndex].domains.remove(subDomain);
  }

  void assignToSubDomain(SubDomainType subDomain, const Codim0Entity& e) {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    indexMap(0)[gt][hostIndex].domains.set(subDomain);
  }

  IndexSetWrapper(const GridImp& grid, HostGridView hostGridView) :
    _grid(grid),
    _hostGridView(hostGridView)
  {
    for (int codim = 0; codim <= dimension; ++codim) {
      _indexMap[codim].reset(new IndexMap());
      _sizeMap[codim].reset(new SizeMap());
    }
  }

  explicit IndexSetWrapper(const ThisType& rhs) :
    _grid(rhs._grid),
    _hostGridView(rhs._hostGridView),
    _indexMap(),
    _sizeMap(),
    _codimSizes(rhs._codimSizes),
    _multiIndexMap(rhs._multiIndexMap)
  {
    for (int codim = 0; codim <= dimension; ++codim) {
      _indexMap[codim].reset(new IndexMap(rhs.indexMap(codim)));
      _sizeMap[codim].reset(new SizeMap(rhs.sizeMap(codim)));
    }
  }

  IndexMap& indexMap(std::size_t codim) {
    return *(_indexMap[codim]);
  }

  SizeMap& sizeMap(std::size_t codim) {
    return *(_sizeMap[codim]);
  }

  const IndexMap& indexMap(std::size_t codim) const {
    return *(_indexMap[codim]);
  }

  const SizeMap& sizeMap(std::size_t codim) const {
    return *(_sizeMap[codim]);
  }

  void reset(bool full) {
    const HostIndexSet& his = _hostGridView.indexSet();
    for (int codim = 0; codim <= dimension; ++codim) {
      if (full) {
        indexMap(codim).swap(IndexMap());
        indexMap(codim).rehash(his.geomTypes(codim).size());
        sizeMap(codim).swap(SizeMap());
        sizeMap(codim).rehash(his.geomTypes(codim).size());
      }
      for (std::vector<GeometryType>::const_iterator it = his.geomTypes(codim).begin(); it != his.geomTypes(codim).end(); ++it) {
        if (full) {
          indexMap(codim)[*it].resize(his.size(*it));
        }
	sizeMap(codim)[*it].assign(0);
      }
    }
    _multiIndexMap.clear();
  }

  void update(LevelIndexSets& levelIndexSets, bool full) {
    const HostIndexSet& his = _hostGridView.indexSet();
    //reset(full);
    for (typename LevelIndexSets::iterator it = levelIndexSets.begin(); it != levelIndexSets.end(); ++it) {
      (*it)->reset(full);
    }
    HostEntityIterator end = _hostGridView.template end<0>();
    IndexMap& im = indexMap(0);
    SizeMap& sm = sizeMap(0);
    for (HostEntityIterator it  = _hostGridView.template begin<0>(); it != end; ++it) {
      const HostEntity& he = *it;
      const GeometryType hgt = he.type();
      IndexType hostIndex = his.index(he);
      MapEntry& me = im[hgt][hostIndex];
      levelIndexSets[he.level()]->indexMap(0)[hgt][levelIndexSets[he.level()]->_hostGridView.indexSet().index(he)].domains.add(me.domains);
      markAncestors(levelIndexSets,HostEntityPointer(he),me.domains);
      updateMapEntry(me,sm[hgt]);
      markSubIndices(he,me.domains,his,GenericReferenceElements<ctype,dimension>::general(hgt),1);
    }
    updateSubIndices(1);
    updatePerCodimSizes();
    for(typename LevelIndexSets::iterator it = levelIndexSets.begin(); it != levelIndexSets.end(); ++it) {
      (*it)->updateLevelIndexSet();
    }
  }


  void updateLevelIndexSet() {
    const HostIndexSet& his = _hostGridView.indexSet();
    HostEntityIterator end = _hostGridView.template end<0>();
    IndexMap& im = indexMap(0);
    SizeMap& sm = sizeMap(0);
    for (HostEntityIterator it  = _hostGridView.template begin<0>(); it != end; ++it) {
      const HostEntity& he = *it;
      const GeometryType hgt = he.type();
      IndexType hostIndex = his.index(he);
      MapEntry& me = im[hgt][hostIndex];
      updateMapEntry(me,sm[hgt]);
      markSubIndices(he,me.domains,his,GenericReferenceElements<ctype,dimension>::general(hgt),1);
    }
    updateSubIndices(1);
    updatePerCodimSizes();
  }

  void updateMapEntry(MapEntry& me, SizeContainer& sizes) {
    switch (me.domains.state()) {
    case SubDomainSet::simpleSet:
      me.index = sizes[*me.domains.begin()]++;
      break;
    case SubDomainSet::multipleSet:
      me.index = _multiIndexMap.size();
      _multiIndexMap.push_back(MultiIndexContainer());
      MultiIndexContainer& mic = _multiIndexMap.back();
      for (typename SubDomainSet::Iterator it = me.domains.begin(); it != me.domains.end(); ++it) {
	mic[*it] = sizes[*it]++;
      }
    }
  }

  void markAncestors(LevelIndexSets& levelIndexSets, HostEntityPointer he, const SubDomainSet& domains) {
    while (he->level() > 0) {
      he = he->father();
      SubDomainSet& fatherDomains = levelIndexSets[he->level()]->indexMap(0)[he->type()][levelIndexSets[he->level()]->_hostGridView.indexSet().index(*he)].domains;
      if (fatherDomains.contains(domains))
        break;
      fatherDomains.add(domains);
    }
  }

  void markSubIndices(const HostEntity& e, const SubDomainSet& domains, const HostIndexSet& his, const GenericReferenceElement<ctype,dimension>& refEl, int codim) {
    const int size = refEl.size(codim);
    IndexMap& im = indexMap(codim);
    for (int i = 0; i < size; ++i) {
      IndexType hostIndex = his.subIndex(e,i,codim);
      GeometryType gt = refEl.type(i,codim);
      im[gt][hostIndex].domains.add(domains);
    }
    if (codim < dimension) {
      markSubIndices(e,domains,his,refEl,codim+1);
    }
  }

  void updateSubIndices(int codim) {
    const typename IndexMap::iterator end = indexMap(codim).end();
    for (typename IndexMap::iterator it = indexMap(codim).begin(); it != end; ++it) {
      const GeometryType gt = it->first;
      std::vector<MapEntry>& indices = it->second;
      SizeContainer& sizes = sizeMap(codim)[gt];
      std::for_each(indices.begin(),indices.end(),boost::bind(&ThisType::updateMapEntry,this,_1,boost::ref(sizes)));
    }
    if (codim < dimension) {
      updateSubIndices(codim+1);
    }
  }

  void updatePerCodimSizes() {
    for (int codim = 0; codim <= dimension; ++codim) {
      _codimSizes[codim].assign(0);
      std::for_each(util::value_iterator(sizeMap(codim).begin()),
		    util::value_iterator(sizeMap(codim).end()),
		    util::collect_elementwise<std::plus<IndexType> >(_codimSizes[codim]));
    }
  }

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_INDEXSETS_HH
