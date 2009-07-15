#ifndef DUNE_MULTIDOMAINGRID_HIERARCHICITERATOR_HH
#define DUNE_MULTIDOMAINGRID_HIERARCHICITERATOR_HH

namespace Dune {

namespace mdgrid {

template<typename GridImp>
class HierarchicIteratorWrapper :
    public EntityPointerWrapper<0,GridImp>
{

  template<class, template<class> class>
  friend class HierarchicIterator;

  template<int, int, typename>
  friend class MakeableEntityWrapper;

  template<int, int, typename>
  friend class EntityWrapper;

  typedef typename GridImp::HostGridType::Traits::template Codim<0>::Entity::HierarchicIterator HostHierarchicIterator;

  explicit HierarchicIteratorWrapper(const HostHierarchicIterator& hostIterator) :
    EntityPointerWrapper<0,GridImp>(hostIterator),
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
