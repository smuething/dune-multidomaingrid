#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_HIERARCHICITERATOR_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_HIERARCHICITERATOR_HH

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<typename GridImp>
class HierarchicIteratorWrapper :
    public EntityPointerWrapper<0,GridImp>
{

  template<class, template<class> class>
  friend class HierarchicIterator;

  template<int, int, typename>
  friend class MakeableEntityWrapper;

  template<int, int, typename>
  friend class EntityWrapper;

  typedef typename GridImp::HostGridType::Traits::template Codim<0>::Entity::HierarchicIterator HostHierarchicIterator;
  typedef typename GridImp::Traits::LevelIndexSet IndexSet;

  explicit HierarchicIteratorWrapper(const GridImp& grid, const HostHierarchicIterator& hostIterator, const HostHierarchicIterator& hostEnd) :
    EntityPointerWrapper<0,GridImp>(hostIterator),
    _grid(grid),
    _hostIterator(hostIterator),
    _hostEnd(hostEnd)
  {
    incrementToNextValidPosition();
  }

  void incrementToNextValidPosition() {
    while(_hostIterator != _end && _grid.containsHostEntity(*_hostIterator)) {
      ++_hostIterator;
    }
  }

  void increment() {
    ++_hostIterator;
    incrementToNextValidPosition();
    this->_entityWrapper.reset(_hostIterator);
  }

  const GridImp& _grid;
  HostHierarchicIterator _hostIterator;
  const HostHierarchicIterator _hostEnd;

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_HIERARCHICITERATOR_HH
