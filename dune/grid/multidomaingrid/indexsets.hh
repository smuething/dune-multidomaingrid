#ifndef DUNE_MULTIDOMAINGRID_INDEXSETS_HH
#define DUNE_MULTIDOMAINGRID_INDEXSETS_HH

#include <map>
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
#include <boost/swap.hpp>

#include <dune/grid/multidomaingrid/utility.hh>

namespace Dune {

namespace mdgrid {

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

template<typename HostGrid, typename MDGridTraits>
class MultiDomainGrid;

//! @cond DEV_DOC

//! \internal
namespace detail {

//! \internal template meta program for assembling the per-codim map vector
template<template<int> class StructForCodim, int codim>
struct buildMap {

  typedef typename mpl::push_back<typename buildMap<StructForCodim,codim-1>::type,
                                  StructForCodim<codim>
                                  >::type type;

};

//! \internal recursion limit for buildMap template meta program
template<template<int> class StructForCodim>
struct buildMap<StructForCodim,0> {

  typedef fusion::vector<StructForCodim<0> > type;

};

//! \internal Helper mechanism for dispatching a call to a possibly non-existent method
/**
 * This trick is necessary because some calls in the indexset (e.g. IndexSet::size(int codim) )
 * pass the codim as a run-time parameter. But as we do not necessarily support all codimensions,
 * we have to make sure that we only call the actual method if there is a data structure to
 * operate on. This is handled by the template parameter doDispatch: The specialisation for
 * doDispatch == false will simply throw an exception.
 */
template<bool doDispatch, typename Impl, typename resulttype>
struct invokeIf {

  typedef resulttype result_type;

  template<int codim>
  result_type invoke() {
    return _impl.invoke<codim>();
  }

  Impl& _impl;

  invokeIf(Impl& impl) :
    _impl(impl)
  {}

};

//! \internal
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

//! \internal Template meta program for dispatching a method to the correctly parameterised template method base on a run-time parameter
template<typename Impl, typename resulttype, typename MDGridTraits, int codim, bool protect = true>
struct dispatchToCodim : public dispatchToCodim<Impl,resulttype,MDGridTraits,codim-1,protect> {

  typedef resulttype result_type;

  result_type dispatch(int cc) {
    if (cc == codim)
      return invokeIf<MDGridTraits::template Codim<codim>::supported,Impl,result_type>(static_cast<Impl&>(*this)).template invoke<codim>();
    return static_cast<dispatchToCodim<Impl,result_type,MDGridTraits,codim-1,protect>&>(*this).dispatch(cc);
  }

};

//! \internal Recursion limit for dispatchToCodim
template<typename Impl, typename resulttype, typename MDGridTraits>
struct dispatchToCodim<Impl,resulttype,MDGridTraits,0,true> {

  typedef resulttype result_type;

  result_type dispatch(int cc) {
    if (cc == 0)
      return invokeIf<MDGridTraits::template Codim<0>::supported,Impl,result_type>(static_cast<Impl&>(*this)).template invoke<0>();
    DUNE_THROW(GridError,"invalid codimension specified");
  }

};

//! \internal Template meta program for dispatching a method to the correctly parameterised template method base on a run-time parameter
template<typename Impl, typename resulttype, typename MDGridTraits, int codim>
struct dispatchToCodim<Impl,resulttype,MDGridTraits,codim,false> : public dispatchToCodim<Impl,resulttype,MDGridTraits,codim-1,false> {

  typedef resulttype result_type;

  result_type dispatch(int cc) {
    if (cc == codim)
      return static_cast<Impl&>(*this).template invoke<codim>();
    return static_cast<dispatchToCodim<Impl,result_type,MDGridTraits,codim-1,false>&>(*this).dispatch(cc);
  }

};

//! \internal Recursion limit for dispatchToCodim
template<typename Impl, typename resulttype, typename MDGridTraits>
struct dispatchToCodim<Impl,resulttype,MDGridTraits,0,false> {

  typedef resulttype result_type;

  result_type dispatch(int cc) {
    if (cc == 0)
      return static_cast<Impl&>(*this).template invoke<0>();
    DUNE_THROW(GridError,"invalid codimension specified");
  }

};

/**
 * \internal Helper structure for protecting calls at non-supported codims, similar to invokeIf.
 *
 * The main difference is that failures are silently ignored and that the called method may not return anything.
 */
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

//! \internal
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

//! \internal Little helper function for casting away constness. Avoids having to spell out the typename.
template<typename T>
T& rw(const T& t) {
  return const_cast<T&>(t);
}

}

//! @endcond

/**
 * Index set for the MultiDomainGrid.
 */
template<typename GridImp, typename HostGridViewType>
class IndexSetWrapper :
    public Dune::IndexSet<GridImp,IndexSetWrapper<GridImp,HostGridViewType>,
                          typename HostGridViewType::IndexSet::IndexType>
{

  template<typename, typename>
  friend class IndexSetWrapper;

  template<typename, typename>
  friend class MultiDomainGrid;

  template<typename, typename>
  friend class subdomain::IndexSetWrapper;

  template<typename, typename, typename, typename>
  friend class SubDomainInterface;

  template<typename>
  friend class SubDomainToSubDomainController;

  template<typename>
  friend class AllInterfacesController;

  typedef IndexSetWrapper<GridImp,HostGridViewType> ThisType;

  typedef typename remove_const<GridImp>::type::HostGridType HostGrid;
  typedef HostGridViewType HostGridView;
  typedef typename HostGridView::IndexSet HostIndexSet;
  typedef typename remove_const<GridImp>::type::ctype ctype;

public:

  typedef typename remove_const<GridImp>::type::MDGridTraits MDGridTraits;
  typedef typename MDGridTraits::template Codim<0>::SubDomainSet SubDomainSet;
  typedef typename MDGridTraits::SubDomainIndexType SubDomainType DUNE_DEPRECATED;
  typedef typename MDGridTraits::SubDomainIndexType SubDomainIndexType;

  typedef typename HostIndexSet::IndexType IndexType;
  static const int dimension = remove_const<GridImp>::type::dimension;
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

  //! Placeholder struct for non-supported codimensions in data structures.
  struct NotSupported {};

  template<int codim>
  struct Containers {

    static const bool supported = remove_const<GridImp>::type::MDGridTraits::template Codim<codim>::supported;

    dune_static_assert(codim > 0 || supported, "index mapping of codimension 0 must be supported!");

    typedef typename SelectType<supported,
                                std::map<GeometryType,std::vector<MapEntry<codim> > >,
                                NotSupported
                                >::Type IndexMap;

    typedef typename SelectType<supported,
                                std::map<GeometryType,typename remove_const<GridImp>::type::MDGridTraits::template Codim<codim>::SizeContainer>,
                                NotSupported
                                >::Type SizeMap;

    typedef typename SelectType<supported,
                                typename remove_const<GridImp>::type::MDGridTraits::template Codim<codim>::SizeContainer,
                                NotSupported
                                >::Type CodimSizeMap;

    typedef typename SelectType<supported,
                                std::vector<typename remove_const<GridImp>::type::MDGridTraits::template Codim<codim>::MultiIndexContainer>,
                                NotSupported
                                >::Type MultiIndexMap;

    IndexMap indexMap;
    SizeMap sizeMap;
    CodimSizeMap codimSizeMap;
    MultiIndexMap multiIndexMap;

    // containers should not be assignable
    const Containers& operator=(const Containers& r);

    Containers(const Containers& r) :
      indexMap(r.indexMap),
      sizeMap(r.sizeMap),
      codimSizeMap(r.codimSizeMap),
      multiIndexMap(r.multiIndexMap) {
    }

    void swap(Containers& rhs) {
      boost::swap(indexMap,rhs.indexMap);
      boost::swap(sizeMap,rhs.sizeMap);
      boost::swap(codimSizeMap,rhs.codimSizeMap);
      boost::swap(multiIndexMap,rhs.multiIndexMap);
    }

    Containers() {}

  };

  typedef typename detail::buildMap<Containers,dimension>::type ContainerMap;

  typedef std::vector<boost::shared_ptr<IndexSetWrapper<GridImp, typename HostGridView::Grid::LevelGridView> > > LevelIndexSets;

  //! Convenience subclass of dispatchToCodim for automatically passing in the MDGridTraits and the dimension
  template<typename Impl,typename result_type, bool protect = true>
  struct dispatchToCodim : public detail::dispatchToCodim<Impl,result_type,MDGridTraits,dimension,protect> {};

public:

  //! Returns the index of the entity with codimension codim.
  template<int codim>
  IndexType index(const typename remove_const<GridImp>::type::Traits::template Codim<codim>::Entity& e) const {
    return _hostGridView.indexSet().index(_grid.hostEntity(e));
  }

  //! Returns the index of the entity.
  template<typename Entity>
  IndexType index(const Entity& e) const {
    return _hostGridView.indexSet().index(_grid.hostEntity(e));
  }

  //! Returns the subdindex of the i-th subentity of e with codimension codim.
  template<int codim>
  IndexType subIndex(const Codim0Entity& e, int i) const {
    return _hostGridView.indexSet().subIndex(_grid.hostEntity(e),i,codim);
  }

  //! Returns the subdindex of the i-th subentity of e with codimension codim.
  IndexType subIndex(const Codim0Entity& e, int i, unsigned int codim) const {
    IndexType r = _hostGridView.indexSet().subIndex(_grid.hostEntity(e),i,codim);
    return r;
  }

  //! Returns a list of all geometry types with codimension codim contained in the grid.
  const std::vector<GeometryType>& geomTypes(int codim) const {
    return _hostGridView.indexSet().geomTypes(codim);
  }

  //! Returns the number of entities with GeometryType type in the grid.
  IndexType size(GeometryType type) const {
    return _hostGridView.indexSet().size(type);
  }

  //! Returns the number of entities with codimension codim in the grid.
  IndexType size(int codim) const {
    return _hostGridView.indexSet().size(codim);
  }

  //! Returns true if the entity is contained in the grid.
  template<typename EntityType>
  bool contains(const EntityType& e) const {
    return _hostGridView.indexSet().contains(_grid.hostEntity(e));
  }

  //! Returns a constant reference to the SubDomainSet of the given entity.
  template<typename EntityType>
  const typename MapEntry<EntityType::codimension>::SubDomainSet& subDomains(const EntityType& e) const {
    return subDomainsForHostEntity(_grid.hostEntity(e));
  }

  //! Returns a constant reference to the SubDomainSet of the given entity with codimension cc.
  //! \tparam cc the codimension of the entity.
  template<int cc>
  const typename MapEntry<cc>::SubDomainSet& subDomains(const typename remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& e) const {
    return subDomainsForHostEntity<cc>(_grid.hostEntity(e));
  }

  //! Returns the index of the entity in a specific subdomain.
  template<class EntityType>
  IndexType index(SubDomainIndexType subDomain, const EntityType& e) const {
    return index<EntityType::codimension>(subDomain,e);
  }

  //! Returns the index of the entity with codimension cc in a specific subdomain.
  //! \tparam the codimension of the entity.
  template<int cc>
  IndexType index(SubDomainIndexType subDomain, const typename remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& e) const {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    const MapEntry<cc>& me = indexMap<cc>().at(gt).at(hostIndex);
    assert(me.domains.contains(subDomain));
    if (me.domains.simple()) {
      return me.index;
    } else {
      return multiIndexMap<cc>()[me.index].at(me.domains.domainOffset(subDomain));
    }
  }

private:

  //! Returns a constant reference to the SubDomainSet of the given host entity.
  template<typename HostEntity>
  const typename MapEntry<HostEntity::codimension>::SubDomainSet& subDomainsForHostEntity(const HostEntity& e) const {
    return subDomainsForHostEntity<HostEntity::codimension>(e);
  }

  //! Returns a constant reference to the SubDomainSet of the given entity with codimension cc.
  //! \tparam cc the codimension of the entity.
  template<int cc>
  const typename MapEntry<cc>::SubDomainSet& subDomainsForHostEntity(const typename remove_const<GridImp>::type::HostGridType::Traits::template Codim<cc>::Entity& he) const {
    return indexMap<cc>().find(he.type())->second[_hostGridView.indexSet().index(he)].domains;
  }

  template<int cc>
  IndexType indexForSubDomain(SubDomainIndexType subDomain, const typename remove_const<GridImp>::type::HostGridType::Traits::template Codim<cc>::Entity& he) const {
    const GeometryType gt = he.type();
    const IndexType hostIndex = _hostGridView.indexSet().index(he);
    const MapEntry<cc>& me = indexMap<cc>().find(gt)->second[hostIndex];
    assert(me.domains.contains(subDomain));
    if (me.domains.simple()) {
      return me.index;
    } else {
      return multiIndexMap<cc>()[me.index].at(me.domains.domainOffset(subDomain));
    }
  }

  //! functor template for retrieving a subindex.
  struct getSubIndexForSubDomain : public dispatchToCodim<getSubIndexForSubDomain,IndexType> {

    template<int codim>
    IndexType invoke() const {
      const MapEntry<codim>& me = _indexSet.indexMap<codim>().at(_gt)[_hostIndex];
      assert(me.domains.contains(_subDomain));
      if (me.domains.simple()) {
        return me.index;
      } else {
        return _indexSet.multiIndexMap<codim>()[me.index].at(me.domains.domainOffset(_subDomain));
      }
    }

    SubDomainIndexType _subDomain;
    GeometryType _gt;
    IndexType _hostIndex;
    const ThisType& _indexSet;

    getSubIndexForSubDomain(SubDomainIndexType subDomain, GeometryType gt, IndexType hostIndex, const ThisType& indexSet) :
      _subDomain(subDomain),
      _gt(gt),
      _hostIndex(hostIndex),
      _indexSet(indexSet)
    {}

  };

  IndexType subIndexForSubDomain(SubDomainIndexType subDomain, const typename remove_const<GridImp>::type::HostGridType::Traits::template Codim<0>::Entity& he, int i, int codim) const {
    return getSubIndexForSubDomain(subDomain,
                                   GenericReferenceElements<ctype,dimension>::general(he.type()).type(i,codim),
                                   _hostGridView.indexSet().subIndex(he,i,codim),
                                   *this).dispatch(codim);
  }

  const std::vector<GeometryType>& geomTypesForSubDomain(SubDomainIndexType subDomain, int codim) const {
    return geomTypes(codim);
  }

  struct getGeometryTypeSizeForSubDomain : public dispatchToCodim<getGeometryTypeSizeForSubDomain,IndexType> {

    template<int codim>
    IndexType invoke() const {
      return _indexSet.sizeMap<codim>().find(_gt)->second[_subDomain];
    }

    SubDomainIndexType _subDomain;
    GeometryType _gt;
    const ThisType& _indexSet;

    getGeometryTypeSizeForSubDomain(SubDomainIndexType subDomain, GeometryType gt, const ThisType& indexSet) :
      _subDomain(subDomain),
      _gt(gt),
      _indexSet(indexSet)
    {}

  };

  IndexType sizeForSubDomain(SubDomainIndexType subDomain, GeometryType type) const {
    return getGeometryTypeSizeForSubDomain(subDomain,type,*this).dispatch(dimension-type.dim());
  }

  struct getCodimSizeForSubDomain : public dispatchToCodim<getCodimSizeForSubDomain,IndexType> {

    template<int codim>
    IndexType invoke() const {
      return _indexSet.codimSizes<codim>()[_subDomain];
    }

    SubDomainIndexType _subDomain;
    const ThisType& _indexSet;

    getCodimSizeForSubDomain(SubDomainIndexType subDomain, const ThisType& indexSet) :
      _subDomain(subDomain),
      _indexSet(indexSet)
    {}

  };

  IndexType sizeForSubDomain(SubDomainIndexType subDomain, int codim) const {
    return getCodimSizeForSubDomain(subDomain,*this).dispatch(codim);
  }

  template<typename EntityType>
  bool containsForSubDomain(SubDomainIndexType subDomain, const EntityType& he) const {
    const GeometryType gt = he.type();
    const IndexType hostIndex = _hostGridView.indexSet().index(he);
    const MapEntry<EntityType::codimension>& me = indexMap<EntityType::codimension>().find(gt)->second[hostIndex];
    return me.domains.contains(subDomain);
  }

public:

  IndexType subIndex(SubDomainIndexType subDomain, const typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i, int codim) const {
    return subIndexForSubDomain(subDomain,_grid.hostEntity(e),i,codim);
  }

  const std::vector<GeometryType>& geomTypes(SubDomainIndexType subDomain, int codim) const {
    return geomTypes(codim);
  }

  IndexType size(SubDomainIndexType subDomain, GeometryType type) const {
    return sizeForSubDomain(subDomain,type);
  }

  IndexType size(SubDomainIndexType subDomain, int codim) const {
    return sizeForSubDomain(subDomain,codim);
  }

  //! Returns true if the entity is contained in a specific subdomain.
  template<typename EntityType>
  bool contains(SubDomainIndexType subDomain, const EntityType& e) const {
    const GeometryType gt = e.type();
    const IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    const MapEntry<EntityType::codimension>& me = indexMap<EntityType::codimension>().find(gt)->second[hostIndex];
    return me.domains.contains(subDomain);
  }

private:

  const GridImp& _grid;
  HostGridView _hostGridView;
  ContainerMap _containers;

  void swap(ThisType& rhs) {
    assert(&_grid == &rhs._grid);
    util::swap(_containers,rhs._containers);
  }

  void addToSubDomain(SubDomainIndexType subDomain, const Codim0Entity& e) {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    indexMap<0>().at(gt)[hostIndex].domains.add(subDomain);
  }

  void removeFromSubDomain(SubDomainIndexType subDomain, const Codim0Entity& e) {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    indexMap<0>()[gt][hostIndex].domains.remove(subDomain);
  }

  void removeFromAllSubDomains(const Codim0Entity& e) {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    indexMap<0>()[gt][hostIndex].domains.clear();
  }

  void assignToSubDomain(SubDomainIndexType subDomain, const Codim0Entity& e) {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    indexMap<0>()[gt][hostIndex].domains.set(subDomain);
  }

  void addToSubDomains(const typename MDGridTraits::template Codim<0>::SubDomainSet& subDomains, const Codim0Entity& e) {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    indexMap<0>().at(gt)[hostIndex].domains.addAll(subDomains);
  }

  IndexSetWrapper(const GridImp& grid, HostGridView hostGridView) :
    _grid(grid),
    _hostGridView(hostGridView)
  {}

  explicit IndexSetWrapper(const ThisType& rhs) :
    _grid(rhs._grid),
    _hostGridView(rhs._hostGridView),
    _containers(rhs._containers)
    {}

  //! Returns the index map for the given codimension.
  //! \tparam cc The requested codimension.
  template<int cc>
  typename Containers<cc>::IndexMap& indexMap() {
    return fusion::at_c<cc>(_containers).indexMap;
  }

  template<int cc>
  typename Containers<cc>::SizeMap& sizeMap() {
    return fusion::at_c<cc>(_containers).sizeMap;
  }

  template<int cc>
  typename Containers<cc>::CodimSizeMap& codimSizes() {
    return fusion::at_c<cc>(_containers).codimSizeMap;
  }

  template<int cc>
  typename Containers<cc>::MultiIndexMap& multiIndexMap() {
    return fusion::at_c<cc>(_containers).multiIndexMap;
  }

  template<int cc>
  const typename Containers<cc>::IndexMap& indexMap() const {
    return fusion::at_c<cc>(_containers).indexMap;
  }

  template<int cc>
  const typename Containers<cc>::SizeMap& sizeMap() const {
    return fusion::at_c<cc>(_containers).sizeMap;
  }

  template<int cc>
  const typename Containers<cc>::CodimSizeMap& codimSizes() const {
    return fusion::at_c<cc>(_containers).codimSizeMap;
  }

  template<int cc>
  const typename Containers<cc>::MultiIndexMap& multiIndexMap() const {
    return fusion::at_c<cc>(_containers).multiIndexMap;
  }

  template<typename Functor>
  void applyToCodims(Functor func) const {
    fusion::for_each(fusion::zip(mpl::range_c<int,0,dimension+1>(),_containers),func);
  }

  template<typename Functor>
  void applyToCodims(Functor func) {
    // we can't just use fusion::zip here as it will slap a const on the ContainerMap type,
    // making it impossible to modify it in the called functor. So we just create the
    // zip_view by hand... (idea lifted from fusion::swap()).
    typedef mpl::range_c<int,0,dimension+1> CodimIndex;
    typedef fusion::vector<CodimIndex&,ContainerMap&> References;
    CodimIndex codimIndex;
    fusion::for_each(fusion::zip_view<References>(References(codimIndex,_containers)),func);
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

    template<int codim>
    void apply(Containers<codim>& c) const {
      for (std::vector<GeometryType>::const_iterator it = _his.geomTypes(codim).begin(); it != _his.geomTypes(codim).end(); ++it) {
        if (_full) {
          // resize index map
          c.indexMap[*it].resize(_his.size(*it));
        }
        // reset SizeMap counter
        std::fill(c.sizeMap[*it].begin(),c.sizeMap[*it].end(),0);
      }
      // clear MultiIndexMap
      c.multiIndexMap.clear();
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
      ContainerMap cm = ContainerMap();
      util::swap(_containers,cm);
    }
    applyToCodims(resetPerCodim(full,his));
  }

  struct updatePerCodimSizes : public applyToCodim<updatePerCodimSizes> {

    template<int codim>
    void apply(Containers<codim>& c) const {
      // reset size for this codim to zero
      std::fill(c.codimSizeMap.begin(),c.codimSizeMap.end(),0);
      // collect per-geometrytype sizes into codim size structure
      std::for_each(util::value_iterator(c.sizeMap.begin()),
		    util::value_iterator(c.sizeMap.end()),
		    util::collect_elementwise<std::plus<IndexType> >(c.codimSizeMap));
    }

  };

  void update(LevelIndexSets& levelIndexSets, bool full) {
    const HostIndexSet& his = _hostGridView.indexSet();
    //reset(full);
    for (typename LevelIndexSets::iterator it = levelIndexSets.begin(); it != levelIndexSets.end(); ++it) {
      (*it)->reset(full);
    }
    HostEntityIterator end = _hostGridView.template end<0>();
    typename Containers<0>::IndexMap& im = indexMap<0>();
    typename Containers<0>::SizeMap& sm = sizeMap<0>();
    for (HostEntityIterator it  = _hostGridView.template begin<0>(); it != end; ++it) {
      const HostEntity& he = *it;
      const GeometryType hgt = he.type();
      IndexType hostIndex = his.index(he);
      MapEntry<0>& me = im[hgt][hostIndex];

      if (_grid.supportLevelIndexSets()) {
        levelIndexSets[he.level()]->indexMap<0>()[hgt][levelIndexSets[he.level()]->_hostGridView.indexSet().index(he)].domains.addAll(me.domains);
        markAncestors(levelIndexSets,HostEntityPointer(he),me.domains);
      }
      updateMapEntry(me,sm[hgt],multiIndexMap<0>());
      applyToCodims(markSubIndices(he,me.domains,his,GenericReferenceElements<ctype,dimension>::general(hgt)));
    }
    applyToCodims(updateSubIndices(*this));
    applyToCodims(updatePerCodimSizes());
    for(typename LevelIndexSets::iterator it = levelIndexSets.begin(); it != levelIndexSets.end(); ++it) {
      (*it)->updateLevelIndexSet();
    }
  }


  void updateLevelIndexSet() {
    const HostIndexSet& his = _hostGridView.indexSet();
    HostEntityIterator end = _hostGridView.template end<0>();
    typename Containers<0>::IndexMap& im = indexMap<0>();
    typename Containers<0>::SizeMap& sm = sizeMap<0>();
    for (HostEntityIterator it  = _hostGridView.template begin<0>(); it != end; ++it) {
      const HostEntity& he = *it;
      const GeometryType hgt = he.type();
      IndexType hostIndex = his.index(he);
      MapEntry<0>& me = im[hgt][hostIndex];
      updateMapEntry(me,sm[hgt],multiIndexMap<0>());
      applyToCodims(markSubIndices(he,me.domains,his,GenericReferenceElements<ctype,dimension>::general(hgt)));
    }
    applyToCodims(updateSubIndices(*this));
    applyToCodims(updatePerCodimSizes());
  }

  template<int codim, typename SizeContainer, typename MultiIndexContainer>
  void updateMapEntry(MapEntry<codim>& me, SizeContainer& sizes, std::vector<MultiIndexContainer>& multiIndexMap) {
    switch (me.domains.state()) {
    case SubDomainSet::emptySet:
      break;
    case SubDomainSet::simpleSet:
      me.index = sizes[*me.domains.begin()]++;
      break;
    case SubDomainSet::multipleSet:
      me.index = multiIndexMap.size();
      multiIndexMap.push_back(MultiIndexContainer());
      MultiIndexContainer& mic = multiIndexMap.back();
      for (typename SubDomainSet::Iterator it = me.domains.begin(); it != me.domains.end(); ++it) {
	mic[me.domains.domainOffset(*it)] = sizes[*it]++;
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

    template<int codim>
    void apply(Containers<codim>& c) const {
      if (codim == 0)
        return;
      const int size = _refEl.size(codim);
      for (int i = 0; i < size; ++i) {
        IndexType hostIndex = _his.subIndex(_he,i,codim);
        GeometryType gt = _refEl.type(i,codim);
        c.indexMap.at(gt)[hostIndex].domains.addAll(_domains);
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

    template<int codim>
    void apply(Containers<codim>& c) const {
      if (codim == 0)
        return;
      const typename Containers<codim>::IndexMap::iterator end = c.indexMap.end();
      for (typename Containers<codim>::IndexMap::iterator it = c.indexMap.begin(); it != end; ++it) {
        const GeometryType gt = it->first;
        typedef typename remove_const<typename Containers<codim>::IndexMap::mapped_type>::type::iterator Iterator;
        const Iterator end2 = it->second.end();
        for (Iterator it2 = it->second.begin(); it2 != end2; ++it2)
          _indexSet.updateMapEntry(*it2,c.sizeMap[gt],c.multiIndexMap);
      }
    }

    ThisType& _indexSet;

    updateSubIndices(ThisType& indexSet) :
      _indexSet(indexSet)
    {}
  };

  //! functor template for retrieving a subindex.
  struct getSupportsCodim : public dispatchToCodim<getSupportsCodim,bool,false> {

    template<int codim>
    bool invoke() const {
      return MDGridTraits::template Codim<codim>::supported && Dune::Capabilities::hasEntity<HostGrid,codim>::v;
    }

  };

  bool supportsCodim(int codim) const
  {
    return getSupportsCodim().dispatch(codim);
  }

  template<typename Impl>
  struct SubDomainSetDataHandleBase
    : public Dune::CommDataHandleIF<Impl,
                                    typename MapEntry<0>::SubDomainSet::DataHandle::DataType
                                    >
  {
    typedef typename MapEntry<0>::SubDomainSet SubDomainSet;
    typedef typename SubDomainSet::DataHandle DataHandle;

    bool fixedsize(int dim, int codim) const
    {
      return DataHandle::fixedsize(dim,codim);
    }

    template<typename Entity>
    std::size_t size(const Entity& e) const
    {
      return MapEntry<Entity::codimension>::SubDomainSet::DataHandle::size(_indexSet.subDomainsForHostEntity(e));
    }

    template<typename MessageBufferImp, typename Entity>
    void gather(MessageBufferImp& buf, const Entity& e) const
    {
      MapEntry<Entity::codimension>::SubDomainSet::DataHandle::gather(buf,_indexSet.subDomainsForHostEntity(e));
    }

    template<typename MessageBufferImp, typename Entity>
    void scatter(MessageBufferImp& buf, const Entity& e, std::size_t n)
    {
      MapEntry<Entity::codimension>::SubDomainSet::DataHandle::scatter(buf,_indexSet.subDomainsForHostEntity(e),n);
    }

    SubDomainSetDataHandleBase(ThisType& indexSet)
      : _indexSet(indexSet)
    {}

    ThisType& _indexSet;

  };

  struct SelectionDataHandle
    : public SubDomainSetDataHandleBase<SelectionDataHandle>
  {

    bool contains(int dim, int codim) const
    {
      return codim == 0;
    }

    SelectionDataHandle(ThisType& indexSet)
      : SubDomainSetDataHandleBase<SelectionDataHandle>(indexSet)
    {}

  };

  void communicateSubDomainSelection()
  {
    SelectionDataHandle dh(*this);
    _hostGridView.template communicate<SelectionDataHandle>(dh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
  }

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_INDEXSETS_HH
