#ifndef DUNE_MULTIDOMAINGRID_ENTITYPOINTER_HH
#define DUNE_MULTIDOMAINGRID_ENTITYPOINTER_HH

#include <dune/grid/common/gridenums.hh>

#include "entity.hh"
#include "subdomaingrid/entity.hh"

namespace Dune {

namespace mdgrid {

template<typename>
class HierarchicIteratorWrapper;

template<typename,int,PartitionIteratorType,typename>
class IteratorWrapper;

template<int codim, typename GridImp>
class EntityPointerWrapper
{

  static const int dim = GridImp::dimension;

  template<int, int, typename>
  friend class subdomain::EntityWrapper;

  template<typename,int,PartitionIteratorType,typename>
  friend class IteratorWrapper;

  template<typename>
  friend class HierarchicIteratorWrapper;

public:

  static const int codimension = codim;

  using Entity            = typename GridImp::template Codim<codim>::Entity;
  using EntityWrapper     = Dune::mdgrid::EntityWrapper<codim,dim,GridImp>;
  using HostEntityPointer = typename GridImp::HostGridType::Traits::template Codim<codim>::EntityPointer;
  using HostEntity        = typename GridImp::HostGridType::Traits::template Codim<codim>::Entity;

  EntityPointerWrapper() = default;

  explicit EntityPointerWrapper(const Entity& e)
    : _hostEntityPointer(GridImp::getRealImplementation(e).hostEntity())
  {}

  explicit EntityPointerWrapper(const HostEntityPointer& hostEntityPointer)
    : _hostEntityPointer(hostEntityPointer)
  {}

  explicit EntityPointerWrapper(const HostEntity& hostEntity)
    : _hostEntityPointer(hostEntity)
  {}

  bool equals(const EntityPointerWrapper& rhs) const {
    return _hostEntityPointer == rhs._hostEntityPointer;
  }

  Entity dereference() const {
    return {EntityWrapper(*_hostEntityPointer)};
  }

  int level() const {
    return _hostEntityPointer.level();
  }

  const HostEntityPointer& hostEntityPointer() const {
    return _hostEntityPointer;
  }

protected:

  HostEntityPointer _hostEntityPointer;

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_ENTITYPOINTER_HH
