#ifndef DUNE_MULTIDOMAINGRID_HIERARCHICITERATOR_HH
#define DUNE_MULTIDOMAINGRID_HIERARCHICITERATOR_HH

namespace Dune {

namespace mdgrid {

template<int codim, int dim, typename GridImp>
class EntityWrapper;

template<int codim, typename GridImp>
class EntityPointerWrapper;

template<typename GridImp>
class HierarchicIteratorWrapper
{

public:

  static const int codimension = 0;

  using Entity               = typename GridImp::template Codim<0>::Entity;
  using EntityWrapper        = Dune::mdgrid::EntityWrapper<0,GridImp::dimension,GridImp>;
  using EntityPointerWrapper = Dune::mdgrid::EntityPointerWrapper<0,GridImp>;
  using EntityPointer        = typename GridImp::template Codim<0>::EntityPointer;
  using HostIterator         = typename GridImp::HostGridType::HierarchicIterator;

  HierarchicIteratorWrapper() = default;

  explicit HierarchicIteratorWrapper(const HostIterator& hostIterator)
    : _hostIterator(hostIterator)
  {}

  void increment() {
    ++_hostIterator;
  }

  bool equals(const HierarchicIteratorWrapper& r) const
  {
    return _hostIterator == r._hostIterator;
  }

  Entity dereference() const
  {
    return {EntityWrapper(*_hostIterator)};
  }

  int level() const
  {
    return _hostIterator.level();
  }

    // TODO: Remove after 2.4
  operator EntityPointer() const
  {
    return {dereference()};
  }

  // TODO: Remove after 2.4
  bool equals(const EntityPointerWrapper& r) const
  {
    return _hostIterator == r._hostEntityPointer;
  }

private:

  HostIterator _hostIterator;

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_HIERARCHICITERATOR_HH
