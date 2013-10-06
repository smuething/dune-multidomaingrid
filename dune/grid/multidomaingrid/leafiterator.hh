#ifndef DUNE_MULTIDOMAINGRID_LEAFITERATOR_HH
#define DUNE_MULTIDOMAINGRID_LEAFITERATOR_HH

#include <dune/grid/common/gridenums.hh>

#include <dune/grid/multidomaingrid/entitypointer.hh>

namespace Dune {

namespace mdgrid {

template<typename HostGrid, typename MDGridTraits>
class MultiDomainGrid;

template<int codim, PartitionIteratorType pitype, typename GridImp>
class LeafIteratorWrapper :
    public EntityPointerWrapper<codim,GridImp>
{

  template<typename, typename>
  friend class MultiDomainGrid;

  template< int cd, class Grid, class IteratorImp >
  friend class Dune::EntityIterator;

  typedef typename GridImp::HostGridType::Traits::template Codim<codim>::template Partition<pitype>::LeafIterator HostLeafIterator;

  explicit LeafIteratorWrapper(const HostLeafIterator& hostIterator)
    : EntityPointerWrapper<codim,GridImp>(pitype_holder<pitype>(),hostIterator)
    , _hostIterator(hostIterator)
  {}

  void increment() {
    ++_hostIterator;
    this->_entityWrapper.reset(_hostIterator);
  }

  HostLeafIterator _hostIterator;

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_LEAFITERATOR_HH
