#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_LEAFITERATOR_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_LEAFITERATOR_HH

#include <dune/grid/common/gridenums.hh>

#include <dune/grid/multidomaingrid/subdomaingrid/entitypointer.hh>

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<typename MDGrid>
class SubDomainGrid;

template<int codim, PartitionIteratorType pitype, typename GridImp>
class LeafIteratorWrapper :
    public EntityPointerWrapper<codim,GridImp>
{

  template<typename MDGrid>
  friend class SubDomainGrid;

  template< int cd, class Grid, class IteratorImp >
  friend class Dune::EntityIterator;

  template<typename,typename>
  friend class Dune::EntityPointer;

  typedef typename GridImp::MultiDomainGrid::Traits::template Codim<codim>::template Partition<pitype>::LeafIterator MultiDomainLeafIterator;
  typedef typename GridImp::Traits::LeafIndexSet IndexSet;

  LeafIteratorWrapper(const IndexSet& indexSet, const MultiDomainLeafIterator multiDomainIterator, const MultiDomainLeafIterator endIterator) :
    EntityPointerWrapper<codim,GridImp>(indexSet._grid,multiDomainIterator),
    _indexSet(indexSet),
    _multiDomainIterator(multiDomainIterator),
    _end(endIterator)
  {
    incrementToNextValidPosition();
  }

  void incrementToNextValidPosition() {
    while(_multiDomainIterator != _end && !_indexSet.containsMultiDomainEntity(*_multiDomainIterator)) {
      ++_multiDomainIterator;
    }
    this->_entityWrapper.reset(_multiDomainIterator);
  }

  void increment() {
    ++_multiDomainIterator;
    incrementToNextValidPosition();
  }


  const LeafIteratorWrapper& operator=(const LeafIteratorWrapper& rhs) {
    assert(_indexSet == rhs._indexSet);
    _multiDomainIterator = rhs._multiDomainIterator;
    this->_entityWrapper.reset(_multiDomainIterator);
    _end = rhs._end;
    return *this;
  }

  const IndexSet& _indexSet;
  MultiDomainLeafIterator _multiDomainIterator;
  MultiDomainLeafIterator _end;

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_LEAFITERATOR_HH
