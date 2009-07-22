#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_LEVELITERATOR_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_LEVELITERATOR_HH

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<int codim, PartitionIteratorType pitype, typename GridImp>
class LevelIteratorWrapper :
    public EntityPointerWrapper<codim,GridImp>
{

  template<typename>
  friend class SubDomainGrid;

  template<int, PartitionIteratorType, class,
	   template<int,PartitionIteratorType,class> class>
  friend class LevelIterator;

  template<typename,typename>
  friend class EntityPointer;

  typedef typename GridImp::HostGridType::Traits::template Codim<codim>::template Partition<pitype>::LevelIterator HostLevelIterator;
  typedef typename GridImp::Traits::LevelIndexSet IndexSet;

  LevelIteratorWrapper(const IndexSet& indexSet, const HostLevelIterator& hostIterator, const HostLevelIterator& endIterator) :
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

  const LevelIteratorWrapper& operator=(const LevelIteratorWrapper& rhs) {
    assert(_indexSet == rhs._indexSet);
    _hostIterator = rhs._hostIterator;
    _end = rhs._end;
  }

  const IndexSet& _indexSet;
  HostLevelIterator _hostIterator;
  HostLevelIterator _end;

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_LEVELITERATOR_HH
