#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_ENTITYPOINTER_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_ENTITYPOINTER_HH

#include <dune/grid/common/gridenums.hh>

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<int, int, typename>
class EntityWrapper;

template<int codim, typename GridImp>
class EntityPointerWrapper
{

  static const int dim = GridImp::dimension;

  template<int, int, typename>
  friend class subdomain::EntityWrapper;

public:

  static const int codimension = codim;

  using Entity                   = typename GridImp::template Codim<codim>::Entity;
  using EntityWrapper            = Dune::mdgrid::subdomain::EntityWrapper<codim,dim,GridImp>;
  using MultiDomainEntityPointer = typename GridImp::MultiDomainGrid::Traits::template Codim<codim>::EntityPointer;
  using MultiDomainEntity        = typename GridImp::MultiDomainGrid::Traits::template Codim<codim>::Entity;

  EntityPointerWrapper()
    : _grid(nullptr)
  {}

  explicit EntityPointerWrapper(const Entity& e)
    : _grid(&GridImp::getRealImplementation(e).grid())
    , _multiDomainEntityPointer(GridImp::getRealImplementation(e).multiDomainEntity())
  {}

  EntityPointerWrapper(const GridImp* grid, const MultiDomainEntityPointer& multiDomainEntityPointer)
    : _grid(grid)
    , _multiDomainEntityPointer(multiDomainEntityPointer)
  {}

  explicit EntityPointerWrapper(const GridImp* grid, const MultiDomainEntity& multiDomainEntity)
    : _grid(grid)
    , _multiDomainEntityPointer(multiDomainEntity)
  {}

  bool equals(const EntityPointerWrapper& rhs) const {
    return _grid == rhs._grid && _multiDomainEntityPointer == rhs._multiDomainEntityPointer;
  }

  Entity dereference() const {
    return {EntityWrapper(_grid,*_multiDomainEntityPointer)};
  }

  int level() const {
    return _multiDomainEntityPointer.level();
  }

private:

  const GridImp* _grid;
  MultiDomainEntityPointer _multiDomainEntityPointer;

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_ENTITYPOINTER_HH
