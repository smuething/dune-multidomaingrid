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

  template<typename HostGrid>
  friend class MultiDomainGrid;

  typedef IndexSetWrapper<GridImp,HostGridViewType> ThisType;

  typedef typename remove_const<GridImp>::type::HostGridType HostGrid;
  typedef HostGridViewType HostGridView;
  typedef typename HostGridView::IndexSet HostIndexSet;
  typedef typename remove_const<GridImp>::type::ctype ctype;

public:

  typedef typename remove_const<GridImp>::type::SubDomainSet SubDomainSet;
  typedef typename HostIndexSet::IndexType IndexType;
  static const int dimension = remove_const<GridImp>::type::dimension;
  typedef typename SubDomainSet::DomainType DomainType;
  static const std::size_t maxSubDomains = SubDomainSet::maxSize;

private:

  typedef typename HostGridView::template Codim<0>::Iterator HostEntityIterator;
  typedef typename HostGridView::template Codim<0>::Entity HostEntity;
  typedef typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity Codim0Entity;

  struct MapEntry {
    SubDomainSet domains;
    IndexType index;
  };

  typedef std::array<IndexType,maxSubDomains> SizeContainer;
  typedef std::unordered_map<GeometryType,std::vector<MapEntry>,util::GeometryTypeHash> IndexMap;
  typedef std::unordered_map<GeometryType,SizeContainer,util::GeometryTypeHash> SizeMap;
  typedef std::array<IndexType,maxSubDomains> MultiIndexContainer;

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

private:

  const GridImp& _grid;
  HostGridView _hostGridView;
  std::array<boost::scoped_ptr<IndexMap>,dimension+1> _indexMap;
  std::array<boost::scoped_ptr<SizeMap>,dimension+1> _sizeMap;
  std::array<std::array<IndexType,maxSubDomains>,dimension+1> _codimSizes;
  std::vector<MultiIndexContainer> _multiIndexMap;


  IndexSetWrapper(const GridImp& grid, HostGridView hostGridView) :
    _grid(grid),
    _hostGridView(hostGridView)
  {
    for (int codim = 0; codim <= dimension; ++codim) {
      _indexMap[codim].reset(new IndexMap());
      _sizeMap[codim].reset(new SizeMap());
    }
  }

  IndexMap& indexMap(std::size_t codim) {
    return *(_indexMap[codim]);
  }

  SizeMap& sizeMap(std::size_t codim) {
    return *(_sizeMap[codim]);
  }

  void update(bool full) {
    const HostIndexSet& his = _hostGridView.indexSet();
    if (full) {
      for (int codim = 0; codim <= dimension; ++codim) {
	indexMap(codim).swap(IndexMap());
	indexMap(codim).rehash(his.geomTypes(codim).size());
	sizeMap(codim).swap(SizeMap());
	sizeMap(codim).rehash(his.geomTypes(codim).size());
      }
    }
    for (int codim = 0; codim <= dimension; ++codim) {
      for (std::vector<GeometryType>::const_iterator it = his.geomTypes(codim).begin(); it != his.geomTypes(codim).end(); ++it) {
	indexMap(codim)[*it].resize(his.size(*it));
	sizeMap(codim)[*it].assign(0);
      }
    }
    _multiIndexMap.clear();
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

  void markSubIndices(const HostEntity& e, const SubDomainSet& domains, const HostIndexSet& his, const GenericReferenceElement<ctype,dimension>& refEl, int codim) {
    const int size = refEl.size(codim);
    IndexMap& im = indexMap(codim);
    for (int i = 0; i < size; ++i) {
      IndexType hostIndex = his.subIndex(e,i,codim);
      GeometryType gt = refEl.type(i,codim);
      im[gt][hostIndex].domains.addAll(domains);
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
