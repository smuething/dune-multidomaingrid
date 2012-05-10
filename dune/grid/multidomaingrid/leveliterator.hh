#ifndef DUNE_MULTIDOMAINGRID_LEVELITERATOR_HH
#define DUNE_MULTIDOMAINGRID_LEVELITERATOR_HH

namespace Dune {

namespace mdgrid {

template<int codim, PartitionIteratorType pitype, typename GridImp>
class LevelIteratorWrapper :
    public EntityPointerWrapper<codim,GridImp>
{

  template<typename,typename>
  friend class MultiDomainGrid;

  template< int cd, class Grid, class IteratorImp >
  friend class EntityIterator;

  typedef typename GridImp::HostGridType::Traits::template Codim<codim>::template Partition<pitype>::LevelIterator HostLevelIterator;

  explicit LevelIteratorWrapper(const HostLevelIterator& hostIterator)
    : EntityPointerWrapper<codim,GridImp>(pitype_holder<pitype>(),hostIterator)
    , _hostIterator(hostIterator)
  {}

  void increment() {
    ++_hostIterator;
    this->_entityWrapper.reset(_hostIterator);
  }

  HostLevelIterator _hostIterator;

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_LEVELITERATOR_HH
