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
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/push_front.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/algorithm/iteration.hpp>
#include <boost/fusion/algorithm/transformation/zip.hpp>

#include <dune/grid/multidomaingrid/utility.hh>

namespace Dune {

namespace mdgrid {

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

template<typename HostGrid, typename MDGridTraits>
class MultiDomainGrid;

namespace detail {

template<template<int> class buildForCodim, int codim>
struct buildMap {

  typedef typename mpl::push_back<typename buildMap<buildForCodim,codim-1>::type,
                                  typename buildForCodim<codim>::type
                                  >::type type;

};

template<template<int> class buildForCodim>
struct buildMap<buildForCodim,0> {

  typedef fusion::vector<typename buildForCodim<0>::type> type;

};

template<bool doDispatch, typename Impl, typename resulttype>
struct invokeIf {

  typedef resulttype result_type;

  template<int codim>
  result_type invoke() {
    _impl.invoke<codim>();
  }

  Impl& _impl;

  invokeIf(Impl& impl) :
    _impl(impl)
  {}

};

template<typename Impl, typename resulttype>
struct invokeIf<false,Impl, resulttype> {

  typedef resulttype result_type;

  template<int codim>
  result_type invoke() {
    DUNE_THROW(GridError,"codim not supported");
  }

  Impl& _impl;

  invokeIf(Impl& impl) :
    _impl(impl)
  {}

};

template<typename Impl, typename resulttype, typename MDGridTraits, int codim>
struct dispatchToCodim : public dispatchToCodim<Impl,resulttype,MDGridTraits,codim-1> {

  typedef resulttype result_type;

  result_type dispatch(int cc) {
    if (cc == codim)
      return invokeIf<MDGridTraits::template Codim<codim>::supported,Impl,result_type>(static_cast<Impl&>(*this)).template invoke<codim>();
    return static_cast<dispatchToCodim<Impl,result_type,MDGridTraits,codim-1>&>(*this).dispatch(cc);
  }

};

template<typename Impl, typename resulttype, typename MDGridTraits>
struct dispatchToCodim<Impl,resulttype,MDGridTraits,0> {

  typedef resulttype result_type;

  result_type dispatch(int cc) {
    if (cc == 0)
      return invokeIf<MDGridTraits::template Codim<0>::supported,Impl,result_type>(static_cast<Impl&>(*this)).template invoke<0>();
    DUNE_THROW(GridError,"invalid codimension specified");
  }

};

template<bool doApply, typename Impl>
struct applyIf {

  template<int codim, typename T>
  void apply(T& t) const {
    _impl.template apply<codim>(t);
  }

  Impl& _impl;

  applyIf(Impl& impl) :
    _impl(impl)
  {}

};

template<typename Impl>
struct applyIf<false,Impl> {

  template<int codim, typename T>
  void apply(T& t) const {
  }

  Impl& _impl;

  applyIf(Impl& impl) :
    _impl(impl)
  {}

};

template<typename T>
T& rw(const T& t) {
  return const_cast<T&>(t);
}

}


template<typename GridImp, typename HostGridViewType>
class IndexSetWrapper :
    public Dune::IndexSet<GridImp,IndexSetWrapper<GridImp,HostGridViewType>,
                          typename HostGridViewType::IndexSet::IndexType>
{

  template<typename, typename>
  friend class IndexSetWrapper;

  template<typename, typename>
  friend class MultiDomainGrid;

  typedef IndexSetWrapper<GridImp,HostGridViewType> ThisType;

  typedef typename remove_const<GridImp>::type::HostGridType HostGrid;
  typedef HostGridViewType HostGridView;
  typedef typename HostGridView::IndexSet HostIndexSet;
  typedef typename remove_const<GridImp>::type::ctype ctype;

public:

  typedef typename remove_const<GridImp>::type::MDGridTraits MDGridTraits;
  typedef typename MDGridTraits::template Codim<0>::SubDomainSet SubDomainSet;
  typedef typename MDGridTraits::SubDomainType SubDomainType;
  typedef typename HostIndexSet::IndexType IndexType;
  static const int dimension = remove_const<GridImp>::type::dimension;
  typedef SubDomainType DomainType;
  static const std::size_t maxSubDomains = SubDomainSet::maxSize;

private:

  typedef typename HostGridView::template Codim<0>::Iterator HostEntityIterator;
  typedef typename HostGridView::template Codim<0>::Entity HostEntity;
  typedef typename HostGridView::template Codim<0>::EntityPointer HostEntityPointer;
  typedef typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity Codim0Entity;

  template<int cc>
  struct MapEntry {

    typedef typename MDGridTraits::template Codim<cc>::SubDomainSet SubDomainSet;

    SubDomainSet domains;
    IndexType index;
  };

  struct NotSupported {};

  template<int codim>
  struct buildIndexMapForCodim {

    static const bool supported = remove_const<GridImp>::type::MDGridTraits::template Codim<codim>::supported;

    dune_static_assert(codim > 0 || supported, "index mapping of codimension 0 must be supported!");

    typedef typename SelectType<supported,
                                std::map<GeometryType,std::vector<MapEntry<codim> > >,
                                NotSupported
                                >::Type type;
  };

  template<int codim>
  struct buildSizeMapForCodim {

    static const bool supported = remove_const<GridImp>::type::MDGridTraits::template Codim<codim>::supported;

    typedef typename SelectType<supported,
                                std::map<GeometryType,typename remove_const<GridImp>::type::MDGridTraits::template Codim<codim>::SizeContainer>,
                                NotSupported
                                >::Type type;
  };

  template<int codim>
  struct buildCodimSizeMapForCodim {

    static const bool supported = remove_const<GridImp>::type::MDGridTraits::template Codim<codim>::supported;

    typedef typename SelectType<supported,
                                typename remove_const<GridImp>::type::MDGridTraits::template Codim<codim>::SizeContainer,
                                NotSupported
                                >::Type type;
  };


  template<int codim>
  struct buildMultiIndexMapForCodim {

    static const bool supported = remove_const<GridImp>::type::MDGridTraits::template Codim<codim>::supported;

    typedef typename SelectType<supported,
                                std::vector<typename remove_const<GridImp>::type::MDGridTraits::template Codim<codim>::MultiIndexContainer>,
                                NotSupported
                                >::Type type;
  };

  typedef std::array<IndexType,maxSubDomains> SizeContainer;
  typedef typename detail::buildMap<buildIndexMapForCodim,dimension>::type IndexMap;
  typedef typename detail::buildMap<buildSizeMapForCodim,dimension>::type SizeMap;
  typedef typename detail::buildMap<buildCodimSizeMapForCodim,dimension>::type CodimSizeMap;
  typedef typename detail::buildMap<buildMultiIndexMapForCodim,dimension>::type MultiIndexMap;

  typedef std::vector<boost::shared_ptr<IndexSetWrapper<GridImp, typename HostGridView::Grid::LevelGridView> > > LevelIndexSets;

  template<typename Impl,typename result_type>
  struct dispatchToCodim : public detail::dispatchToCodim<Impl,result_type,MDGridTraits,dimension> {};

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
  const typename MapEntry<EntityType::codimension>::SubDomainSet& subDomains(const EntityType& e) const {
    return subDomains<EntityType::codimension>(e);
  }

  template<int cc>
  const typename MapEntry<cc>::SubDomainSet& subDomains(const typename remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& e) const {
    return indexMap<cc>().find(e.type())->second[_hostGridView.indexSet().index(_grid.hostEntity(e))].domains;
  }

  template<class EntityType>
  IndexType index(DomainType subDomain, const EntityType& e) const {
    return index<EntityType::codimension>(subDomain,e);
  }

  template<int cc>
  IndexType index(DomainType subDomain, const typename remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& e) const {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    const MapEntry<cc>& me = indexMap<cc>().at(gt).at(hostIndex);
    assert(me.domains.contains(subDomain));
    if (me.domains.simple()) {
      return me.index;
    } else {
      return multiIndexMap<cc>()[me.index].at(subDomain);
    }
  }

  template<int cc>
  IndexType indexForSubDomain(DomainType subDomain, const typename remove_const<GridImp>::type::HostGridType::Traits::template Codim<cc>::Entity& he) const {
    const GeometryType gt = he.type();
    const IndexType hostIndex = _hostGridView.indexSet().index(he);
    const MapEntry<cc>& me = indexMap<cc>().find(gt)->second[hostIndex];
    assert(me.domains.contains(subDomain));
    if (me.domains.simple()) {
      return me.index;
    } else {
      return multiIndexMap<cc>()[me.index].at(subDomain);
    }
  }

  struct getSubIndexForSubDomain : public dispatchToCodim<getSubIndexForSubDomain,IndexType> {

    template<int codim>
    IndexType invoke() const {
      const MapEntry<codim>& me = _indexSet.indexMap<codim>().at(_gt)[_hostIndex];
      assert(me.domains.contains(_subDomain));
      if (me.domains.simple()) {
        return me.index;
      } else {
        return _indexSet.multiIndexMap<codim>()[me.index][_subDomain];
      }
    }

    DomainType _subDomain;
    GeometryType _gt;
    IndexType _hostIndex;
    const ThisType& _indexSet;

    getSubIndexForSubDomain(DomainType subDomain, GeometryType gt, IndexType hostIndex, const ThisType& indexSet) :
      _subDomain(subDomain),
      _gt(gt),
      _hostIndex(hostIndex),
      _indexSet(indexSet)
    {}

  };

  IndexType subIndexForSubDomain(DomainType subDomain, const typename remove_const<GridImp>::type::HostGridType::Traits::template Codim<0>::Entity& he, int i, int codim) const {
    return getSubIndexForSubDomain(subDomain,
                                   GenericReferenceElements<ctype,dimension>::general(he.type()).type(i,codim),
                                   _hostGridView.indexSet().subIndex(he,i,codim),
                                   *this).dispatch(codim);
  }

  const std::vector<GeometryType>& geomTypesForSubDomain(DomainType subDomain, int codim) const {
    return geomTypes(codim);
  }

  struct getGeometryTypeSizeForSubDomain : public dispatchToCodim<getGeometryTypeSizeForSubDomain,IndexType> {

    template<int codim>
    IndexType invoke() const {
      return _indexSet.sizeMap<codim>().find(_gt)->second[_subDomain];
    }

    DomainType _subDomain;
    GeometryType _gt;
    const ThisType& _indexSet;

    getGeometryTypeSizeForSubDomain(DomainType subDomain, GeometryType gt, const ThisType& indexSet) :
      _subDomain(subDomain),
      _gt(gt),
      _indexSet(indexSet)
    {}

  };

  IndexType sizeForSubDomain(DomainType subDomain, GeometryType type) const {
    return getGeometryTypeSizeForSubDomain(subDomain,type,*this).dispatch(dimension-type.dim());
  }

  struct getCodimSizeForSubDomain : public dispatchToCodim<getCodimSizeForSubDomain,IndexType> {

    template<int codim>
    IndexType invoke() const {
      return _indexSet.codimSizes<codim>()[_subDomain];
    }

    DomainType _subDomain;
    const ThisType& _indexSet;

    getCodimSizeForSubDomain(DomainType subDomain, const ThisType& indexSet) :
      _subDomain(subDomain),
      _indexSet(indexSet)
    {}

  };

  IndexType sizeForSubDomain(DomainType subDomain, int codim) const {
    return getCodimSizeForSubDomain(subDomain,*this).dispatch(codim);
  }

  template<typename EntityType>
  bool containsForSubDomain(DomainType subDomain, const EntityType& he) const {
    const GeometryType gt = he.type();
    const IndexType hostIndex = _hostGridView.indexSet().index(he);
    const MapEntry<EntityType::codimension>& me = indexMap<EntityType::codimension>().find(gt)->second[hostIndex];
    return me.domains.contains(subDomain);
  }

  IndexType subIndex(DomainType subDomain, const typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i, int codim) const {
    return subIndexForSubDomain(subDomain,_grid.hostEntity(e),i,codim);
  }

  const std::vector<GeometryType>& geomTypes(DomainType subDomain, int codim) const {
    return geomTypes(codim);
  }

  IndexType size(DomainType subDomain, GeometryType type) const {
    return sizeForSubDomain(subDomain,type);
  }

  IndexType size(DomainType subDomain, int codim) const {
    return sizeForSubDomain(subDomain,codim);
  }

  template<typename EntityType>
  bool contains(DomainType subDomain, const EntityType& e) const {
    const GeometryType gt = e.type();
    const IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    const MapEntry<EntityType::codimension>& me = indexMap<EntityType::codimension>().find(gt)->second[hostIndex];
    return me.domains.contains(subDomain);
  }

private:

  const GridImp& _grid;
  HostGridView _hostGridView;
  IndexMap _indexMap;
  SizeMap _sizeMap;
  CodimSizeMap _codimSizes;
  MultiIndexMap _multiIndexMap;

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
    indexMap<0>().at(gt)[hostIndex].domains.add(subDomain);
  }

  void removeFromSubDomain(SubDomainType subDomain, const Codim0Entity& e) {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    indexMap<0>()[gt][hostIndex].domains.remove(subDomain);
  }

  void assignToSubDomain(SubDomainType subDomain, const Codim0Entity& e) {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    indexMap<0>()[gt][hostIndex].domains.set(subDomain);
  }

  IndexSetWrapper(const GridImp& grid, HostGridView hostGridView) :
    _grid(grid),
    _hostGridView(hostGridView)
  {}

  explicit IndexSetWrapper(const ThisType& rhs) :
    _grid(rhs._grid),
    _hostGridView(rhs._hostGridView),
    _indexMap(rhs._indexMap),
    _sizeMap(rhs._sizeMap),
    _codimSizes(rhs._codimSizes),
    _multiIndexMap(rhs._multiIndexMap)
  {}

  template<int cc>
  typename fusion::result_of::at_c<IndexMap,cc>::type indexMap() {
    return fusion::at_c<cc>(_indexMap);
  }

  template<int cc>
  typename fusion::result_of::at_c<SizeMap,cc>::type sizeMap() {
    return fusion::at_c<cc>(_sizeMap);
  }

  template<int cc>
  typename fusion::result_of::at_c<CodimSizeMap,cc>::type codimSizes() {
    return fusion::at_c<cc>(_codimSizes);
  }

  template<int cc>
  typename fusion::result_of::at_c<MultiIndexMap,cc>::type multiIndexMap() {
    return fusion::at_c<cc>(_multiIndexMap);
  }

  template<int cc>
  typename fusion::result_of::at_c<const IndexMap,cc>::type indexMap() const {
    return fusion::at_c<cc>(_indexMap);
  }

  template<int cc>
  typename fusion::result_of::at_c<const SizeMap,cc>::type sizeMap() const {
    return fusion::at_c<cc>(_sizeMap);
  }

  template<int cc>
  typename fusion::result_of::at_c<const CodimSizeMap,cc>::type codimSizes() const {
    return fusion::at_c<cc>(_codimSizes);
  }

  template<int cc>
  typename fusion::result_of::at_c<const MultiIndexMap,cc>::type multiIndexMap() const {
    return fusion::at_c<cc>(_multiIndexMap);
  }

  template<typename Sequence, typename Functor>
  void applyToCodims(Sequence seq, Functor func) const {
    fusion::for_each(fusion::zip(mpl::range_c<int,0,dimension+1>(),seq),func);
  }

  template<typename Sequence, typename Functor>
  void applyToCodims(Sequence seq, Functor func) {
    fusion::for_each(fusion::zip(mpl::range_c<int,0,dimension+1>(),seq),func);
  }

  template<typename Impl>
  struct applyToCodim {

    template<typename T>
    void operator()(T t) const {
      typedef mpl::int_<fusion::result_of::value_at_c<T,0>::type::value> codim;
      detail::applyIf<MDGridTraits::template Codim<codim::value>::supported,const Impl>(static_cast<const Impl&>(*this)).template apply<codim::value>(fusion::at_c<1>(t));
    }

  };

  struct resetPerCodim : public applyToCodim<resetPerCodim> {

    template<int codim, typename T>
    void apply(T& t) const {
      typedef typename std::remove_const<typename std::remove_reference<typename fusion::result_of::value_at_c<T,0>::type>::type>::type IndexMap;
      typedef typename std::remove_const<typename std::remove_reference<typename fusion::result_of::value_at_c<T,1>::type>::type>::type SizeMap;
      typedef typename std::remove_const<typename std::remove_reference<typename fusion::result_of::value_at_c<T,3>::type>::type>::type MultiIndexMap;
      IndexMap& indexMap = const_cast<IndexMap&>(fusion::at_c<0>(t));
      SizeMap& sizeMap = const_cast<SizeMap&>(fusion::at_c<1>(t));
      MultiIndexMap& multiIndexMap = const_cast<MultiIndexMap&>(fusion::at_c<3>(t));
      for (std::vector<GeometryType>::const_iterator it = _his.geomTypes(codim).begin(); it != _his.geomTypes(codim).end(); ++it) {
        if (_full) {
          // resize index map
          indexMap[*it].resize(_his.size(*it));
        }
        // reset SizeMap counter
        sizeMap[*it].assign(0);
      }
      // clear MultiIndexMap
      multiIndexMap.clear();
    }

    resetPerCodim(bool full, const HostIndexSet& his) :
      _full(full),
      _his(his)
    {}

    const bool _full;
    const HostIndexSet& _his;
  };

  void reset(bool full) {
    const HostIndexSet& his = _hostGridView.indexSet();
    if (full) {
      IndexMap him;
      SizeMap hsm;
      std::swap(_indexMap,him);
      std::swap(_sizeMap,hsm);
    }
    applyToCodims(fusion::zip(_indexMap,_sizeMap,_codimSizes,_multiIndexMap),resetPerCodim(full,his));
  }

  struct updatePerCodimSizes : public applyToCodim<updatePerCodimSizes> {

    template<int codim, typename T>
    void apply(T& t) const {
      // reset size for this codim to zero
      typedef typename std::remove_const<typename std::remove_reference<typename fusion::result_of::value_at_c<T,0>::type>::type>::type PerCodimSize;
      typedef typename std::remove_const<typename std::remove_reference<typename fusion::result_of::value_at_c<T,1>::type>::type>::type SizeMap;
      PerCodimSize& perCodimSize = const_cast<PerCodimSize&>(fusion::at_c<0>(t));
      SizeMap& sizeMap = const_cast<SizeMap&>(fusion::at_c<1>(t));
      perCodimSize.assign(0);
      // collect per-geometrytype sizes into codim size structure
      std::for_each(util::value_iterator(sizeMap.begin()),
		    util::value_iterator(sizeMap.end()),
		    util::collect_elementwise<std::plus<IndexType> >(perCodimSize));
    }

  };

  void update(LevelIndexSets& levelIndexSets, bool full) {
    const HostIndexSet& his = _hostGridView.indexSet();
    //reset(full);
    for (typename LevelIndexSets::iterator it = levelIndexSets.begin(); it != levelIndexSets.end(); ++it) {
      (*it)->reset(full);
    }
    HostEntityIterator end = _hostGridView.template end<0>();
    typename fusion::result_of::at_c<IndexMap,0>::type im = indexMap<0>();
    typename fusion::result_of::at_c<SizeMap,0>::type sm = sizeMap<0>();
    for (HostEntityIterator it  = _hostGridView.template begin<0>(); it != end; ++it) {
      const HostEntity& he = *it;
      const GeometryType hgt = he.type();
      IndexType hostIndex = his.index(he);
      MapEntry<0>& me = im[hgt][hostIndex];
      levelIndexSets[he.level()]->indexMap<0>()[hgt][levelIndexSets[he.level()]->_hostGridView.indexSet().index(he)].domains.addAll(me.domains);
      markAncestors(levelIndexSets,HostEntityPointer(he),me.domains);
      updateMapEntry(me,sm[hgt],multiIndexMap<0>());
      applyToCodims(fusion::zip(_indexMap,_sizeMap),markSubIndices(he,me.domains,his,GenericReferenceElements<ctype,dimension>::general(hgt)));
    }
    applyToCodims(fusion::zip(_indexMap,_sizeMap,_multiIndexMap),updateSubIndices(*this));
    applyToCodims(fusion::zip(_codimSizes,_sizeMap),updatePerCodimSizes());
    for(typename LevelIndexSets::iterator it = levelIndexSets.begin(); it != levelIndexSets.end(); ++it) {
      (*it)->updateLevelIndexSet();
    }
  }


  void updateLevelIndexSet() {
    const HostIndexSet& his = _hostGridView.indexSet();
    HostEntityIterator end = _hostGridView.template end<0>();
    typename fusion::result_of::at_c<IndexMap,0>::type im = indexMap<0>();
    typename fusion::result_of::at_c<SizeMap,0>::type sm = sizeMap<0>();
    for (HostEntityIterator it  = _hostGridView.template begin<0>(); it != end; ++it) {
      const HostEntity& he = *it;
      const GeometryType hgt = he.type();
      IndexType hostIndex = his.index(he);
      MapEntry<0>& me = im[hgt][hostIndex];
      updateMapEntry(me,sm[hgt],multiIndexMap<0>());
      applyToCodims(fusion::zip(_indexMap,_sizeMap),markSubIndices(he,me.domains,his,GenericReferenceElements<ctype,dimension>::general(hgt)));
    }
    applyToCodims(fusion::zip(_indexMap,_sizeMap,_multiIndexMap),updateSubIndices(*this));
    applyToCodims(fusion::zip(_codimSizes,_sizeMap),updatePerCodimSizes());
  }

  template<int codim, typename SizeContainer, typename MultiIndexContainer>
  void updateMapEntry(MapEntry<codim>& me, SizeContainer& sizes, std::vector<MultiIndexContainer>& multiIndexMap) {
    switch (me.domains.state()) {
    case SubDomainSet::simpleSet:
      me.index = sizes[*me.domains.begin()]++;
      break;
    case SubDomainSet::multipleSet:
      me.index = multiIndexMap.size();
      multiIndexMap.push_back(MultiIndexContainer());
      MultiIndexContainer& mic = multiIndexMap.back();
      for (typename SubDomainSet::Iterator it = me.domains.begin(); it != me.domains.end(); ++it) {
	mic[*it] = sizes[*it]++;
      }
    }
  }

  template<typename SubDomainSet>
  void markAncestors(LevelIndexSets& levelIndexSets, HostEntityPointer he, const SubDomainSet& domains) {
    while (he->level() > 0) {
      he = he->father();
      SubDomainSet& fatherDomains = levelIndexSets[he->level()]->indexMap<0>()[he->type()][levelIndexSets[he->level()]->_hostGridView.indexSet().index(*he)].domains;
      if (fatherDomains.containsAll(domains))
        break;
      fatherDomains.addAll(domains);
    }
  }

  struct markSubIndices : public applyToCodim<const markSubIndices> {

    template<int codim, typename T>
    void apply(T& t) const {
      if (codim == 0)
        return;
      typedef typename std::remove_const<typename std::remove_reference<typename fusion::result_of::value_at_c<T,0>::type>::type>::type IndexSet;
      IndexSet& indexSet = const_cast<IndexSet&>(fusion::at_c<0>(t));
      const int size = _refEl.size(codim);
      for (int i = 0; i < size; ++i) {
        IndexType hostIndex = _his.subIndex(_he,i,codim);
        GeometryType gt = _refEl.type(i,codim);
        indexSet.at(gt)[hostIndex].domains.addAll(_domains);
      }
    }

    typedef typename MapEntry<0>::SubDomainSet& DomainSet;

    const HostEntity& _he;
    const DomainSet& _domains;
    const HostIndexSet& _his;
    const GenericReferenceElement<ctype,dimension>& _refEl;

    markSubIndices(const HostEntity& he, const DomainSet& domains, const HostIndexSet& his, const GenericReferenceElement<ctype,dimension>& refEl) :
      _he(he),
      _domains(domains),
      _his(his),
      _refEl(refEl)
    {}

  };

  struct updateSubIndices : public applyToCodim<const updateSubIndices> {

    template<int codim, typename T>
    void apply(T& t) const {
      if (codim == 0)
        return;
      typedef typename std::remove_const<typename std::remove_reference<typename fusion::result_of::value_at_c<T,0>::type>::type>::type IndexMap;
      IndexMap& indexMap = const_cast<IndexMap&>(fusion::at_c<0>(t));
      typedef typename std::remove_const<typename std::remove_reference<typename fusion::result_of::value_at_c<T,1>::type>::type>::type SizeMap;
      SizeMap& sizeMap = const_cast<SizeMap&>(fusion::at_c<1>(t));
      typedef typename std::remove_const<typename std::remove_reference<typename fusion::result_of::value_at_c<T,2>::type>::type>::type MultiIndexMap;
      MultiIndexMap& multiIndexMap = const_cast<MultiIndexMap&>(fusion::at_c<2>(t));
      const typename IndexMap::iterator end = indexMap.end();
      for (typename IndexMap::iterator it = indexMap.begin(); it != end; ++it) {
        const GeometryType gt = it->first;
        typedef typename remove_const<typename IndexMap::mapped_type>::type::iterator Iterator;
        const Iterator end2 = it->second.end();
        for (Iterator it2 = it->second.begin(); it2 != end2; ++it2)
          _indexSet.updateMapEntry(*it2,sizeMap[gt],multiIndexMap);
      }
    }

    ThisType& _indexSet;

    updateSubIndices(ThisType& indexSet) :
      _indexSet(indexSet)
    {}
  };

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_INDEXSETS_HH
