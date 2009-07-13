#ifndef DUNE_MULTIDOMAINGRID_LEAFITERATOR_HH
#define DUNE_MULTIDOMAINGRID_LEAFITERATOR_HH

namespace Dune {

namespace mdgrid {

template<int codim, PartitionIteratorType pitype, typename GridImp>
class LeafIteratorWrapper :
    public EntityPointerWrapper<codim,GridImp>
{

  typedef typename GridImp::HostGridType::Traits::template Codim<codim>::template Partition<pitype>::LeafIterator HostLeafIterator;

  explicit LeafIteratorWrapper(const HostLeafIterator& hostIterator) :
    _hostIterator(hostIterator)
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
