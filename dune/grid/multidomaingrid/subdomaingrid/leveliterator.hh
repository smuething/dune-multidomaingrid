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
  friend class MultiDomainGrid;

  template<int, PartitionIteratorType, class,
	   template<int,PartitionIteratorType,class> class>
  friend class LevelIterator;

  typedef typename GridImp::HostGridType::Traits::template Codim<codim>::template Partition<pitype>::LevelIterator HostLevelIterator;
  typedef typename GridImp::Traits::LevelIndexSet IndexSet;

  LevelIteratorWrapper(const IndexSet& indexSet, const HostLevelIterator& hostIterator, const HostLevelIterator& endIterator) :
    EntityPointerWrapper<codim,GridImp>(hostIterator),
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
  }

  void increment() {
    ++_hostIterator;
    incrementToNextValidPosition();
    this->_entityWrapper.reset(_hostIterator);
  }

  const IndexSet& _indexSet;
  HostLeveIterator _hostIterator;
  const HostLevelIterator _end;

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_LEVELITERATOR_HH
