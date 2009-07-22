#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_LEAFITERATOR_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_LEAFITERATOR_HH

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

  template<int, PartitionIteratorType, class,
	   template<int,PartitionIteratorType,class> class>
  friend class LeafIterator;

  template<typename,typename>
  friend class EntityPointer;

  typedef typename GridImp::HostGridType::Traits::template Codim<codim>::template Partition<pitype>::LeafIterator HostLeafIterator;
  typedef typename GridImp::Traits::LeafIndexSet IndexSet;

  LeafIteratorWrapper(const IndexSet& indexSet, const HostLeafIterator& hostIterator, const HostLeafIterator& endIterator) :
    EntityPointerWrapper<codim,GridImp>(indexSet._grid,hostIterator),
    _indexSet(indexSet),
    _hostIterator(hostIterator),
    _end(endIterator)
  {
    incrementToNextValidPosition();
  }

  void incrementToNextValidPosition() {
    while(_hostIterator != _end && !_indexSet.containsHostEntity(*_hostIterator)) {
      ++_hostIterator;
    }
    this->_entityWrapper.reset(_hostIterator);
  }

  void increment() {
    ++_hostIterator;
    incrementToNextValidPosition();
  }


  const LeafIteratorWrapper& operator=(const LeafIteratorWrapper& rhs) {
    assert(_indexSet == rhs._indexSet);
    _hostIterator = rhs._hostIterator;
    _end = rhs._end;
  }

  const IndexSet& _indexSet;
  HostLeafIterator _hostIterator;
  HostLeafIterator _end;

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_LEAFITERATOR_HH
