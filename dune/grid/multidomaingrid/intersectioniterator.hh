#ifndef DUNE_MULTIDOMAINGRID_INTERSECTIONITERATOR_HH
#define DUNE_MULTIDOMAINGRID_INTERSECTIONITERATOR_HH

#include <dune/grid/common/intersectioniterator.hh>

#include <dune/grid/multidomaingrid/hostgridaccessor.hh>

#include <dune/grid/multidomaingrid/subdomaingrid/intersectioniterator.hh>

namespace Dune {

namespace mdgrid {

template<typename GridImp, typename HostIntersectionIterator_>
class IntersectionIteratorWrapper {

  template<class, class, class>
  friend class Dune::IntersectionIterator;

  template<class, class>
  friend class Dune::Intersection;

  template<typename,typename,typename>
  friend class subdomain::IntersectionIteratorWrapper;

  template<typename,PartitionIteratorType>
  friend class LevelGridView;

  template<typename,PartitionIteratorType>
  friend class LeafGridView;

  using HostIntersectionIterator = HostIntersectionIterator_;
  using HostIntersection         = typename HostIntersectionIterator::Intersection;
  using IntersectionWrapper      = Dune::mdgrid::IntersectionWrapper<GridImp,HostIntersection>;
  using Intersection             = Dune::Intersection<GridImp,IntersectionWrapper>;

protected:

  IntersectionIteratorWrapper() = default;

  explicit IntersectionIteratorWrapper(const HostIntersectionIterator& hostIterator)
    : _hostIterator(hostIterator)
  {}

  explicit IntersectionIteratorWrapper(HostIntersectionIterator&& hostIterator)
    : _hostIterator(std::move(hostIterator))
  {}

private:

  const HostIntersectionIterator& hostIntersectionIterator() const {
    return _hostIterator;
  }

  bool equals(const IntersectionIteratorWrapper& rhs) const {
    return _hostIterator == rhs._hostIterator;
  }

  void increment() {
    ++_hostIterator;
  }

  Intersection dereference() const {
    return {IntersectionWrapper(*_hostIterator)};
  }

  HostIntersectionIterator _hostIterator;

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_INTERSECTIONITERATOR_HH
