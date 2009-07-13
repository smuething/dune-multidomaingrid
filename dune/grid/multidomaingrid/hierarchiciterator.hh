#ifndef DUNE_MULTIDOMAINGRID_HIERARCHICITERATOR_HH
#define DUNE_MULTIDOMAINGRID_HIERARCHICITERATOR_HH

namespace Dune {

namespace mdgrid {

template<typename GridImp>
class HierarchicIteratorWrapper :
    public EntityPointerWrapper<0,GridImp>
{

  typedef typename GridImp::HostGridType::Traits::template Codim<0>::Entity::HierarchicIterator HostHierarchicIterator;

  explicit HierarchicIteratorWrapper(const HostHierarchicIterator& hostIterator) :
    _hostIterator(hostIterator)
  {}

  void increment() {
    ++_hostIterator;
    this->_entityWrapper.reset(_hostIterator);
  }

  HostHierarchicIterator _hostIterator;

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_HIERARCHICITERATOR_HH
