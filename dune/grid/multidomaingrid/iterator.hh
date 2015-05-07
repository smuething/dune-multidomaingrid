#ifndef DUNE_GRID_MULTIDOMAINGRID_ITERATOR_HH
#define DUNE_GRID_MULTIDOMAINGRID_ITERATOR_HH

#include <utility>

#include <dune/grid/common/gridenums.hh>

namespace Dune {

namespace mdgrid {

template<int codim, int dim, typename GridImp>
class EntityWrapper;

template<int codim, typename GridImp>
class EntityPointerWrapper;

template<typename HostGridView, int codim, PartitionIteratorType pitype, typename GridImp>
class IteratorWrapper
{

  template<typename, typename>
  friend class MultiDomainGrid;

  template< int cd, class Grid, class IteratorImp >
  friend class Dune::EntityIterator;

  template< class Grid, class IteratorImp >
  friend class Dune::EntityPointer;

  static const int codimension = codim;

  using HostIterator  = typename HostGridView::template Codim<codim>::template Partition<pitype>::Iterator;
  using Entity        = typename GridImp::template Codim<codim>::Entity;
  using EntityPointer = typename GridImp::template Codim<codim>::EntityPointer;
  using EntityWrapper = Dune::mdGrid::EntityWrapper<cd,GridImp::dimension,GridImp>;

  IteratorWrapper() = default;

  explicit IteratorWrapper(const HostIterator& hostIterator)
    : _hostIterator(hostIterator)
  {}

  explicit IteratorWrapper(HostIterator&& hostIterator)
    : _hostIterator(std::move(hostIterator))
  {}

  void increment() {
    ++_hostIterator;
  }

  bool equals(const IteratorWrapper& r) const
  {
    return _hostIterator == r._hostIterator;
  }

  Entity dereference() const
  {
    return {EntityWrapper(*_hostIterator)};
  }

  int level() const
  {
    return _hostIterator.level();
  }

  // TODO: Remove after 2.4
  operator EntityPointer() const
  {
    return EntityPointer(dereference());
  }

  // TODO: Remove after 2.4
  bool equals(const EntityPointerWrapper<codim,GridImp>& r) const
  {
    return _hostIterator == r._hostEntityPointer;
  }

  HostIterator _hostIterator;

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_GRID_MULTIDOMAINGRID_ITERATOR_HH
