#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_HIERARCHICITERATOR_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_HIERARCHICITERATOR_HH

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<int,int,typename>
class EntityWrapper;

template<int,typename>
class EntityPointerWrapper;

template<typename GridImp>
class HierarchicIteratorWrapper
{

  template< int cd, class Grid, class IteratorImp >
  friend class Dune::EntityIterator;

  template<int, int, typename>
  friend class EntityWrapper;

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename, typename>
  friend class Dune::EntityPointer;

  using MultiDomainIterator = typename GridImp::MultiDomainGrid::template Codim<0>::Entity::HierarchicIterator;

public:

  static const int codimension = 0;

  using EntityWrapper        = Dune::mdgrid::subdomain::EntityWrapper<0,GridImp::dimension,GridImp>;
  using Entity               = typename GridImp::template Codim<0>::Entity;
  using EntityPointerWrapper = Dune::mdgrid::subdomain::EntityPointerWrapper<0,GridImp>;
  using EntityPointer        = typename GridImp::template Codim<0>::EntityPointer;


  HierarchicIteratorWrapper()
    : _grid(nullptr)
  {}

  HierarchicIteratorWrapper(
    const GridImp* grid,
    const MultiDomainIterator& multiDomainIterator,
    const MultiDomainIterator& multiDomainEnd
    )
    : _grid(grid)
    , _multiDomainIterator(multiDomainIterator)
    , _multiDomainEnd(multiDomainEnd)
  {
    incrementToNextValidPosition();
  }

  bool equals(const HierarchicIteratorWrapper& r) const
  {
    return _grid == r._grid && _multiDomainIterator == r._multiDomainIterator;
  }

  Entity dereference() const
  {
    return {EntityWrapper(_grid,*_multiDomainIterator)};
  }

  void incrementToNextValidPosition() {
    while(_multiDomainIterator != _multiDomainEnd && !_grid->containsMultiDomainEntity(*_multiDomainIterator))
      {
        ++_multiDomainIterator;
      }
  }

  void increment() {
    ++_multiDomainIterator;
    incrementToNextValidPosition();
  }


  // TODO: Remove after 2.4
  operator EntityPointer() const
  {
    return {dereference()};
  }

  // TODO: Remove after 2.4
  bool equals(const EntityPointerWrapper& r) const
  {
    return _grid == r._grid && _multiDomainIterator == r._multiDomainEntityPointer;
  }


private:

  const GridImp* _grid;
  MultiDomainIterator _multiDomainIterator;
  MultiDomainIterator _multiDomainEnd;

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_HIERARCHICITERATOR_HH
