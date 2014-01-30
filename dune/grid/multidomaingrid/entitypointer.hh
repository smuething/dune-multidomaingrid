#ifndef DUNE_MULTIDOMAINGRID_ENTITYPOINTER_HH
#define DUNE_MULTIDOMAINGRID_ENTITYPOINTER_HH

#include <dune/grid/common/gridenums.hh>

#include "entity.hh"
#include "subdomaingrid/entity.hh"

namespace Dune {

namespace mdgrid {

template<PartitionIteratorType>
struct pitype_holder
{};

template<int codim, typename GridImp>
class EntityPointerWrapper
{

  static const int dim = GridImp::dimension;

  template<int, int, typename>
  friend class subdomain::EntityWrapper;

  struct Invalid {};
  struct InvalidHierarchic {};

  template<PartitionIteratorType pitype>
  struct HostLeafIterator
  {
    typedef typename GridImp::HostGridType::template Codim<codim>::template Partition<pitype>::LeafIterator type;
  };

  template<PartitionIteratorType pitype>
  struct HostLevelIterator
  {
    typedef typename GridImp::HostGridType::template Codim<codim>::template Partition<pitype>::LevelIterator _type;

    typedef typename conditional<
      !is_same<
        typename HostLeafIterator<pitype>::type,
        _type
        >::value,
      _type,
      Invalid
      >::type type;
  };

public:

  typedef EntityPointerWrapper EntityPointerImp;

  static const int codimension = codim;

  typedef typename GridImp::Traits::template Codim<codim>::Entity Entity;
  typedef EntityPointerWrapper<codim,GridImp> Base;

  typedef typename GridImp::HostGridType::Traits::template Codim<codim>::EntityPointer HostEntityPointer;
  typedef typename GridImp::HostGridType::Traits::HierarchicIterator HostHierarchicIterator;

  typedef typename GridImp::HostGridType::Traits::template Codim<codim>::Entity HostEntity;

  explicit EntityPointerWrapper(const HostEntityPointer& hostEntityPointer) :
    _entityWrapper(hostEntityPointer)
  {}

  template<PartitionIteratorType pitype>
  EntityPointerWrapper(pitype_holder<pitype> pth, const typename HostLeafIterator<pitype>::type& hostEntityPointer)
    : _entityWrapper(hostEntityPointer)
  {}

  template<PartitionIteratorType pitype>
  EntityPointerWrapper(pitype_holder<pitype> pth, const typename HostLevelIterator<pitype>::type& hostEntityPointer)
    : _entityWrapper(hostEntityPointer)
  {}

  explicit EntityPointerWrapper(typename SelectType<codim==0,const HostHierarchicIterator&,InvalidHierarchic>::Type hostEntityPointer) :
    _entityWrapper(hostEntityPointer)
  {}

  explicit EntityPointerWrapper(const HostEntity& hostEntity) :
    _entityWrapper(hostEntity)
  {}

  explicit EntityPointerWrapper(const EntityWrapper<codim,dim,GridImp>& entity) :
    _entityWrapper(entity._hostEntityPointer)
  {}

  bool equals(const EntityPointerWrapper<codim,GridImp>& rhs) const {
    return _entityWrapper.hostEntityPointer() == rhs._entityWrapper.hostEntityPointer();
  }

  Entity& dereference() const {
    return _entityWrapper;
  }

  void compactify() {
    _entityWrapper.compactify();
  }

  int level() const {
    return _entityWrapper.level();
  }

  const HostEntityPointer& hostEntityPointer() const {
    return _entityWrapper.hostEntityPointer();
  }

protected:

  mutable MakeableEntityWrapper<codim,dim,GridImp> _entityWrapper;

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_ENTITYPOINTER_HH
