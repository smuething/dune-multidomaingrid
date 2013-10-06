#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_HIERARCHICITERATOR_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_HIERARCHICITERATOR_HH

#include <dune/grid/multidomaingrid/entity.hh>

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<typename GridImp>
class HierarchicIteratorWrapper :
    public EntityPointerWrapper<0,GridImp>
{

  template< int cd, class Grid, class IteratorImp >
  friend class ::Dune::EntityIterator;

  template<int, int, typename>
  friend class MakeableEntityWrapper;

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename, typename>
  friend class ::Dune::EntityPointer;

  typedef typename GridImp::MDGridType::Traits::template Codim<0>::Entity::HierarchicIterator MultiDomainHierarchicIterator;
  typedef typename GridImp::Traits::LevelIndexSet IndexSet;

  explicit HierarchicIteratorWrapper(const GridImp& grid, const MultiDomainHierarchicIterator& multiDomainIterator, const MultiDomainHierarchicIterator& multiDomainEnd) :
    EntityPointerWrapper<0,GridImp>(grid,multiDomainIterator),
    _grid(grid),
    _multiDomainIterator(multiDomainIterator),
    _multiDomainEnd(multiDomainEnd)
  {
    incrementToNextValidPosition();
  }

  void incrementToNextValidPosition() {
    while(_multiDomainIterator != _multiDomainEnd && !_grid.containsMultiDomainEntity(*_multiDomainIterator)) {
      ++_multiDomainIterator;
    }
    this->_entityWrapper.reset(_multiDomainIterator);
  }

  void increment() {
    ++_multiDomainIterator;
    incrementToNextValidPosition();
  }

  HierarchicIteratorWrapper& operator=(const HierarchicIteratorWrapper& rhs) {
    assert(_grid == rhs._grid);
    _multiDomainIterator = rhs._multiDomainIterator;
    _multiDomainEnd = rhs._multiDomainEnd;
    static_cast<EntityPointerWrapper<0,GridImp>&>(*this) = static_cast<const EntityPointerWrapper<0,GridImp>&>(rhs);
    return *this;
  }

  const GridImp& _grid;
  MultiDomainHierarchicIterator _multiDomainIterator;
  MultiDomainHierarchicIterator _multiDomainEnd;

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_HIERARCHICITERATOR_HH
