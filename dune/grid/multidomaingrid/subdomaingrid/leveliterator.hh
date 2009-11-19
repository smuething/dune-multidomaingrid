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

  typedef typename GridImp::MDGridType::Traits::template Codim<codim>::template Partition<pitype>::LevelIterator MultiDomainLevelIterator;
  typedef typename GridImp::Traits::LevelIndexSet IndexSet;

  LevelIteratorWrapper(const IndexSet& indexSet, const MultiDomainLevelIterator& multiDomainIterator, const MultiDomainLevelIterator& endIterator) :
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

  const LevelIteratorWrapper& operator=(const LevelIteratorWrapper& rhs) {
    assert(_indexSet == rhs._indexSet);
    _multiDomainIterator = rhs._multiDomainIterator;
    _end = rhs._end;
    return *this;
  }

  const IndexSet& _indexSet;
  MultiDomainLevelIterator _multiDomainIterator;
  MultiDomainLevelIterator _end;

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_LEVELITERATOR_HH
