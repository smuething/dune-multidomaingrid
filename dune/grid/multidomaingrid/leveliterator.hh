#ifndef DUNE_MULTIDOMAINGRID_LEVELITERATOR_HH
#define DUNE_MULTIDOMAINGRID_LEVELITERATOR_HH

namespace Dune {

namespace mdgrid {

template<int codim, PartitionIteratorType pitype, typename GridImp>
class LevelIteratorWrapper :
    public EntityPointerWrapper<codim,GridImp>
{

  typedef typename GridImp::HostGridType::Traits::template Codim<codim>::template Partition<pitype>::LevelIterator HostLevelIterator;

  explicit LeafIteratorWrapper(const HostLevelIterator& hostIterator) :
    _hostIterator(hostIterator)
  {}

  void increment() {
    ++_hostIterator;
    _hostEntity.reset(_hostIterator);
  }

  HostLevelIterator _hostIterator;

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_LEVELITERATOR_HH
