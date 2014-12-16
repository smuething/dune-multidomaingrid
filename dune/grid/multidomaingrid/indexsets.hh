#ifndef DUNE_MULTIDOMAINGRID_INDEXSETS_HH
#define DUNE_MULTIDOMAINGRID_INDEXSETS_HH

#include <map>
#include <unordered_map>
#include <vector>
#include <array>
#include <algorithm>
#include <type_traits>
#include <tuple>

#include <dune/common/shared_ptr.hh>
#include <dune/grid/common/exceptions.hh>
#include <dune/grid/common/indexidset.hh>

#include <dune/grid/multidomaingrid/utility.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/indexsets.hh>

#include <dune/typetree/utility.hh>

namespace Dune {

namespace mdgrid {

template<typename HostGrid, typename MDGridTraits>
class MultiDomainGrid;

//! @cond DEV_DOC

//! \internal
namespace detail {

template<template<int> class StructForCodim, typename Codims>
struct _buildMap;

//! \internal template meta program for assembling the per-codim map vector
template<template<int> class StructForCodim, std::size_t... codim>
struct _buildMap<StructForCodim,TypeTree::index_pack<codim...> > {

  typedef std::tuple<StructForCodim<codim>... > type;

};

template<template<int> class StructForCodim, int dimension>
struct buildMap {

  typedef typename _buildMap<
    StructForCodim,
    decltype(TypeTree::index_range<dimension+1>())
    >::type type;

};

//! \internal Helper mechanism for dispatching a call to a possibly non-existent method
/**
 * This trick is necessary because some calls in the indexset (e.g. IndexSet::size(int codim) )
 * pass the codim as a run-time parameter. But as we do not necessarily support all codimensions,
 * we have to make sure that we only call the actual method if there is a data structure to
 * operate on. This is handled by the template parameter doDispatch: The specialisation for
 * doDispatch == false will simply throw an exception.
 */
template<bool doDispatch, typename Impl, typename resulttype,bool alternate_dispatch>
struct invokeIf;

template<typename Impl, typename resulttype, bool alternate_dispatch>
struct invokeIf<true,Impl,resulttype,alternate_dispatch> {

  typedef resulttype result_type;

  template<int codim>
  result_type invoke() {
    return _impl.template invoke<codim>();
  }

  Impl& _impl;

  invokeIf(Impl& impl) :
    _impl(impl)
  {}

};

//! \internal
template<typename Impl, typename resulttype>
struct invokeIf<false,Impl, resulttype,false> {

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

template<typename Impl, typename resulttype>
struct invokeIf<false,Impl, resulttype,true> {

  typedef resulttype result_type;

  template<int codim>
  result_type invoke() {
    return _impl.template invoke_unsupported<codim>();
  }

  Impl& _impl;

  invokeIf(Impl& impl) :
    _impl(impl)
  {}

};

//! \internal Template meta program for dispatching a method to the correctly parameterised template method base on a run-time parameter
template<typename Impl, typename resulttype, typename MDGridTraits, int codim, bool protect = true, bool alternate_dispatch = false>
struct dispatchToCodim
  : public dispatchToCodim<Impl,resulttype,MDGridTraits,codim-1,protect,alternate_dispatch> {

  typedef resulttype result_type;

  result_type dispatch(int cc) {
    if (cc == codim)
      return invokeIf<MDGridTraits::template Codim<codim>::supported,Impl,result_type,alternate_dispatch>(static_cast<Impl&>(*this)).template invoke<codim>();
    return static_cast<dispatchToCodim<Impl,result_type,MDGridTraits,codim-1,protect,alternate_dispatch>&>(*this).dispatch(cc);
  }

};

//! \internal Recursion limit for dispatchToCodim
template<typename Impl, typename resulttype, typename MDGridTraits, bool alternate_dispatch>
struct dispatchToCodim<Impl,resulttype,MDGridTraits,0,true,alternate_dispatch> {

  typedef resulttype result_type;

  result_type dispatch(int cc) {
    if (cc == 0)
      return invokeIf<MDGridTraits::template Codim<0>::supported,Impl,result_type,alternate_dispatch>(static_cast<Impl&>(*this)).template invoke<0>();
    DUNE_THROW(GridError,"invalid codimension specified");
  }

};

//! \internal Template meta program for dispatching a method to the correctly parameterised template method base on a run-time parameter
template<typename Impl, typename resulttype, typename MDGridTraits, int codim, bool alternate_dispatch>
struct dispatchToCodim<Impl,resulttype,MDGridTraits,codim,false,alternate_dispatch>
  : public dispatchToCodim<Impl,resulttype,MDGridTraits,codim-1,false,alternate_dispatch> {

  typedef resulttype result_type;

  result_type dispatch(int cc) {
    if (cc == codim)
      return static_cast<Impl&>(*this).template invoke<codim>();
    return static_cast<dispatchToCodim<Impl,result_type,MDGridTraits,codim-1,false,alternate_dispatch>&>(*this).dispatch(cc);
  }

};

//! \internal Recursion limit for dispatchToCodim
template<typename Impl, typename resulttype, typename MDGridTraits, bool alternate_dispatch>
struct dispatchToCodim<Impl,resulttype,MDGridTraits,0,false,alternate_dispatch> {

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
                          typename HostGridViewType::IndexSet::IndexType,
                          typename HostGridViewType::IndexSet::Types>
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

  typedef Dune::IndexSet<
    GridImp,
    IndexSetWrapper<
      GridImp,
      HostGridViewType
      >,
    typename HostGridViewType::IndexSet::IndexType,
    typename HostGridViewType::IndexSet::Types
    > BaseT;

public:

  typedef typename BaseT::Types Types;

  typedef typename remove_const<GridImp>::type::MDGridTraits MDGridTraits;
  typedef typename MDGridTraits::template Codim<0>::SubDomainSet SubDomainSet;
  typedef typename MDGridTraits::SubDomainIndex SubDomainIndex;


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

    static_assert((codim > 0 || supported), "index mapping of codimension 0 must be supported!");

    typedef typename conditional<supported,
                                 std::map<GeometryType,std::vector<MapEntry<codim> > >,
                                 NotSupported
                                 >::type IndexMap;

    typedef typename conditional<supported,
                                 std::map<GeometryType,typename remove_const<GridImp>::type::MDGridTraits::template Codim<codim>::SizeContainer>,
                                 NotSupported
                                 >::type SizeMap;

    typedef typename conditional<supported,
                                 typename remove_const<GridImp>::type::MDGridTraits::template Codim<codim>::SizeContainer,
                                 NotSupported
                                 >::type CodimSizeMap;

    typedef typename conditional<supported,
                                 std::vector<typename remove_const<GridImp>::type::MDGridTraits::template Codim<codim>::MultiIndexContainer>,
                                 NotSupported
                                 >::type MultiIndexMap;

    IndexMap indexMap;
    SizeMap sizeMap;
    CodimSizeMap codimSizeMap;
    MultiIndexMap multiIndexMap;

    // containers should not be assignable...
    Containers& operator=(const Containers&) = delete;

    // ...but must be movable
    Containers& operator=(Containers&&) = default;

    // make sure the compiler generates all necessary constructors...
    Containers(const Containers&) = default;
    Containers(Containers&&) = default;
    Containers() = default;

  };

  typedef typename detail::buildMap<Containers,dimension>::type ContainerMap;

  typedef std::vector<shared_ptr<IndexSetWrapper<GridImp, typename HostGridView::Grid::LevelGridView> > > LevelIndexSets;

  //! Convenience subclass of dispatchToCodim for automatically passing in the MDGridTraits and the dimension
  template<typename Impl,typename result_type, bool protect = true, bool alternate_dispatch = false>
  struct dispatchToCodim : public detail::dispatchToCodim<Impl,result_type,MDGridTraits,dimension,protect,alternate_dispatch> {};

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
  const std::vector<GeometryType>& DUNE_DEPRECATED_MSG("Use types() instead") geomTypes(int codim) const {
    return _hostGridView.indexSet().geomTypes(codim);
  }

    //! Returns a list of all geometry types with codimension codim contained in the grid.
  Types types(int codim) const {
    return _hostGridView.indexSet().types(codim);
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
  IndexType index(SubDomainIndex subDomain, const EntityType& e) const {
    return index<EntityType::codimension>(subDomain,e);
  }

  //! Returns the index of the entity with codimension cc in a specific subdomain.
  //! \tparam the codimension of the entity.
  template<int cc>
  IndexType index(SubDomainIndex subDomain, const typename remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& e) const {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    const MapEntry<cc>& me = indexMap<cc>().at(gt).at(hostIndex);
    assert(me.domains.contains(subDomain));
    if (me.domains.simple()) {
      return me.index;
    } else {
      return multiIndexMap<cc>()[me.index][me.domains.domainOffset(subDomain)];
    }
  }

private:

  //! Returns a mutable reference to the SubDomainSet of the given entity.
  template<typename EntityType>
  typename MapEntry<EntityType::codimension>::SubDomainSet& subDomains(const EntityType& e) {
    return subDomainsForHostEntity(_grid.hostEntity(e));
  }

  //! Returns a mutable reference to the SubDomainSet of the given entity with codimension cc.
  //! \tparam cc the codimension of the entity.
  template<int cc>
  typename MapEntry<cc>::SubDomainSet& subDomains(const typename remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& e) {
    return subDomainsForHostEntity<cc>(_grid.hostEntity(e));
  }

  //! Returns a mutable reference to the SubDomainSet of the given host entity.
  template<typename HostEntity>
  typename MapEntry<HostEntity::codimension>::SubDomainSet& subDomainsForHostEntity(const HostEntity& e) {
    return subDomainsForHostEntity<HostEntity::codimension>(e);
  }

  //! Returns a mutable reference to the SubDomainSet of the given entity with codimension cc.
  //! \tparam cc the codimension of the entity.
  template<int cc>
  typename MapEntry<cc>::SubDomainSet& subDomainsForHostEntity(const typename remove_const<GridImp>::type::HostGridType::Traits::template Codim<cc>::Entity& he) {
    return indexMap<cc>().find(he.type())->second[_hostGridView.indexSet().index(he)].domains;
  }

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
  IndexType indexForSubDomain(SubDomainIndex subDomain, const typename remove_const<GridImp>::type::HostGridType::Traits::template Codim<cc>::Entity& he) const {
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
        return _indexSet.multiIndexMap<codim>()[me.index][me.domains.domainOffset(_subDomain)];
      }
    }

    SubDomainIndex _subDomain;
    GeometryType _gt;
    IndexType _hostIndex;
    const ThisType& _indexSet;

    getSubIndexForSubDomain(SubDomainIndex subDomain, GeometryType gt, IndexType hostIndex, const ThisType& indexSet) :
      _subDomain(subDomain),
      _gt(gt),
      _hostIndex(hostIndex),
      _indexSet(indexSet)
    {}

  };

  IndexType subIndexForSubDomain(SubDomainIndex subDomain, const typename remove_const<GridImp>::type::HostGridType::Traits::template Codim<0>::Entity& he, int i, int codim) const {
    return getSubIndexForSubDomain(subDomain,
                                   ReferenceElements<ctype,dimension>::general(he.type()).type(i,codim),
                                   _hostGridView.indexSet().subIndex(he,i,codim),
                                   *this).dispatch(codim);
  }

  const std::vector<GeometryType>&
  DUNE_DEPRECATED_MSG("Use typesForSubDomain() instead")
  geomTypesForSubDomain(SubDomainIndex subDomain, int codim) const {
    return geomTypes(codim);
  }

  Types typesForSubDomain(SubDomainIndex subDomain, int codim) const {
    return types(codim);
  }

  struct getGeometryTypeSizeForSubDomain : public dispatchToCodim<getGeometryTypeSizeForSubDomain,IndexType,true,true> {

    template<int codim>
    IndexType invoke() const {
      return _indexSet.sizeMap<codim>().find(_gt)->second[_subDomain];
    }

    template<int codim>
    IndexType invoke_unsupported() const {
      return 0;
    }

    SubDomainIndex _subDomain;
    GeometryType _gt;
    const ThisType& _indexSet;

    getGeometryTypeSizeForSubDomain(SubDomainIndex subDomain, GeometryType gt, const ThisType& indexSet) :
      _subDomain(subDomain),
      _gt(gt),
      _indexSet(indexSet)
    {}

  };

  IndexType sizeForSubDomain(SubDomainIndex subDomain, GeometryType type) const {
    return getGeometryTypeSizeForSubDomain(subDomain,type,*this).dispatch(dimension-type.dim());
  }

  struct getCodimSizeForSubDomain : public dispatchToCodim<getCodimSizeForSubDomain,IndexType,true,true> {

    template<int codim>
    IndexType invoke() const {
      return _indexSet.codimSizes<codim>()[_subDomain];
    }

    template<int codim>
    IndexType invoke_unsupported() const {
      return 0;
    }

    SubDomainIndex _subDomain;
    const ThisType& _indexSet;

    getCodimSizeForSubDomain(SubDomainIndex subDomain, const ThisType& indexSet) :
      _subDomain(subDomain),
      _indexSet(indexSet)
    {}

  };

  IndexType sizeForSubDomain(SubDomainIndex subDomain, int codim) const {
    return getCodimSizeForSubDomain(subDomain,*this).dispatch(codim);
  }

  template<typename EntityType>
  bool containsForSubDomain(SubDomainIndex subDomain, const EntityType& he) const {
    const GeometryType gt = he.type();
    const IndexType hostIndex = _hostGridView.indexSet().index(he);
    const MapEntry<EntityType::codimension>& me = indexMap<EntityType::codimension>().find(gt)->second[hostIndex];
    return me.domains.contains(subDomain);
  }

public:

  IndexType subIndex(SubDomainIndex subDomain, const typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i, int codim) const {
    return subIndexForSubDomain(subDomain,_grid.hostEntity(e),i,codim);
  }

  const std::vector<GeometryType>& DUNE_DEPRECATED_MSG("Use types() instead") geomTypes(SubDomainIndex subDomain, int codim) const {
    return geomTypes(codim);
  }

  Types types(SubDomainIndex subDomain, int codim) const {
    return types(codim);
  }

  IndexType size(SubDomainIndex subDomain, GeometryType type) const {
    return sizeForSubDomain(subDomain,type);
  }

  IndexType size(SubDomainIndex subDomain, int codim) const {
    return sizeForSubDomain(subDomain,codim);
  }

  //! Returns true if the entity is contained in a specific subdomain.
  template<typename EntityType>
  bool contains(SubDomainIndex subDomain, const EntityType& e) const {
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
    std::swap(_containers,rhs._containers);
  }

  void addToSubDomain(SubDomainIndex subDomain, const Codim0Entity& e) {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    indexMap<0>().at(gt)[hostIndex].domains.add(subDomain);
  }

  void removeFromSubDomain(SubDomainIndex subDomain, const Codim0Entity& e) {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    indexMap<0>()[gt][hostIndex].domains.remove(subDomain);
  }

  void removeFromAllSubDomains(const Codim0Entity& e) {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    indexMap<0>()[gt][hostIndex].domains.clear();
  }

  void assignToSubDomain(SubDomainIndex subDomain, const Codim0Entity& e) {
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
    return std::get<cc>(_containers).indexMap;
  }

  template<int cc>
  typename Containers<cc>::SizeMap& sizeMap() {
    return std::get<cc>(_containers).sizeMap;
  }

  template<int cc>
  typename Containers<cc>::CodimSizeMap& codimSizes() {
    return std::get<cc>(_containers).codimSizeMap;
  }

  template<int cc>
  typename Containers<cc>::MultiIndexMap& multiIndexMap() {
    return std::get<cc>(_containers).multiIndexMap;
  }

  template<int cc>
  const typename Containers<cc>::IndexMap& indexMap() const {
    return std::get<cc>(_containers).indexMap;
  }

  template<int cc>
  const typename Containers<cc>::SizeMap& sizeMap() const {
    return std::get<cc>(_containers).sizeMap;
  }

  template<int cc>
  const typename Containers<cc>::CodimSizeMap& codimSizes() const {
    return std::get<cc>(_containers).codimSizeMap;
  }

  template<int cc>
  const typename Containers<cc>::MultiIndexMap& multiIndexMap() const {
    return std::get<cc>(_containers).multiIndexMap;
  }

  template<typename Functor>
  void applyToCodims(Functor func) const {
    TypeTree::apply_to_tuple(
      _containers,
      func,
      TypeTree::apply_to_tuple_policy::pass_index()
      );
  }

  template<typename Functor>
  void applyToCodims(Functor func) {
    TypeTree::apply_to_tuple(
      _containers,
      func,
      TypeTree::apply_to_tuple_policy::pass_index()
      );
  }

  template<typename Impl>
  struct applyToCodim {

    template<typename I, typename T>
    void operator()(I i, T&& t) const {
      const int codim = I::value;
      detail::applyIf<MDGridTraits::template Codim<codim>::supported,const Impl>(static_cast<const Impl&>(*this)).template apply<codim>(std::forward<T>(t));
    }

  };

  struct resetPerCodim : public applyToCodim<resetPerCodim> {

    template<int codim>
    void apply(Containers<codim>& c) const {
      // setup per-codim sizemap
      _traits.template setupSizeContainer<codim>(c.codimSizeMap);
      for (auto gt : _his.types(codim)) {
        if (_full) {
          // resize index map
          c.indexMap[gt].resize(_his.size(gt));
        } else if (codim > 0) {
          // clear out marked state for codim > 0 (we cannot keep the old
          // state for subentities, as doing so will leave stale entries if
          // elements are removed from a subdomain
          for (auto& mapEntry : c.indexMap[gt])
            mapEntry.domains.clear();
        }
        // setup / reset SizeMap counter
        _traits.template setupSizeContainer<codim>(c.sizeMap[gt]);
        std::fill(c.sizeMap[gt].begin(),c.sizeMap[gt].end(),0);
      }
      // clear MultiIndexMap
      c.multiIndexMap.clear();
    }

    resetPerCodim(bool full, const HostIndexSet& his, const MDGridTraits& traits) :
      _full(full),
      _his(his),
      _traits(traits)
    {}

    const bool _full;
    const HostIndexSet& _his;
    const MDGridTraits& _traits;
  };

  void reset(bool full) {
    const HostIndexSet& his = _hostGridView.indexSet();
    if (full) {
      ContainerMap cm = ContainerMap();
      std::swap(_containers,cm);
    }
    applyToCodims(resetPerCodim(full,his,_grid._traits));
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

    this->communicateSubDomainSelection();

    HostEntityIterator end = _hostGridView.template end<0>();
    typename Containers<0>::IndexMap& im = indexMap<0>();
    typename Containers<0>::SizeMap& sm = sizeMap<0>();
    for (HostEntityIterator it  = _hostGridView.template begin<0>(); it != end; ++it) {
      const HostEntity& he = *it;
      const GeometryType hgt = he.type();
      IndexType hostIndex = his.index(he);
      MapEntry<0>& me = im[hgt][hostIndex];

      if (_grid.supportLevelIndexSets()) {
        levelIndexSets[he.level()]->template indexMap<0>()[hgt][levelIndexSets[he.level()]->_hostGridView.indexSet().index(he)].domains.addAll(me.domains);
        markAncestors(levelIndexSets,HostEntityPointer(he),me.domains);
      }
      updateMapEntry(me,sm[hgt],multiIndexMap<0>());
      applyToCodims(markSubIndices(he,me.domains,his,ReferenceElements<ctype,dimension>::general(hgt)));
    }

    propagateBorderEntitySubDomains();

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

    communicateSubDomainSelection();

    for (HostEntityIterator it  = _hostGridView.template begin<0>(); it != end; ++it) {
      const HostEntity& he = *it;
      const GeometryType hgt = he.type();
      IndexType hostIndex = his.index(he);
      MapEntry<0>& me = im[hgt][hostIndex];
      updateMapEntry(me,sm[hgt],multiIndexMap<0>());
      applyToCodims(markSubIndices(he,me.domains,his,ReferenceElements<ctype,dimension>::general(hgt)));
    }

    propagateBorderEntitySubDomains();

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
      SubDomainSet& fatherDomains = levelIndexSets[he->level()]->template indexMap<0>()[he->type()][levelIndexSets[he->level()]->_hostGridView.indexSet().index(*he)].domains;
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
    DomainSet& _domains;
    const HostIndexSet& _his;
    const ReferenceElement<ctype,dimension>& _refEl;

    markSubIndices(const HostEntity& he, DomainSet& domains, const HostIndexSet& his, const ReferenceElement<ctype,dimension>& refEl) :
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

  struct BorderPropagationDataHandle
    : public SubDomainSetDataHandleBase<BorderPropagationDataHandle>
  {

    bool contains(int dim, int codim) const
    {
      return codim > 0 && this->_indexSet.supportsCodim(codim);
    }

    BorderPropagationDataHandle(ThisType& indexSet)
      : SubDomainSetDataHandleBase<BorderPropagationDataHandle>(indexSet)
    {}

  };


  void communicateSubDomainSelection()
  {
    SelectionDataHandle dh(*this);
    _hostGridView.template communicate<SelectionDataHandle>(dh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
  }

  void propagateBorderEntitySubDomains()
  {
    BorderPropagationDataHandle dh(*this);
    _hostGridView.communicate(dh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
  }

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_INDEXSETS_HH
