#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_INTERSECTIONITERATOR_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_INTERSECTIONITERATOR_HH

#include <dune/grid/common/intersectioniterator.hh>

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<
  typename GridImp,
  typename IndexSet,
  typename MultiDomainIntersectionIterator
  >
class IntersectionIteratorWrapper {

  template<class, class, class>
  friend class Dune::IntersectionIterator;

  template<class, class>
  friend class Dune::Intersection;

  using MultiDomainIntersection = typename MultiDomainIntersectionIterator::Intersection;
  using IntersectionWrapper     = Dune::mdgrid::subdomain::IntersectionWrapper<
    GridImp,
    IndexSet,
    MultiDomainIntersection
    >;

public:

  using Intersection            = Dune::Intersection<GridImp,IntersectionWrapper>;

  IntersectionIteratorWrapper()
    : _indexSet(nullptr)
  {}

  IntersectionIteratorWrapper(const IndexSet* indexSet, const MultiDomainIntersectionIterator& multiDomainIterator)
    : _indexSet(indexSet)
    , _multiDomainIterator(multiDomainIterator)
  {}

  const typename GridImp::MultiDomainGrid::template ReturnImplementationType<MultiDomainIntersectionIterator>::ImplementationType::HostIntersectionIterator& hostIntersectionIterator() const {
    return GridImp::MultiDomainGrid::getRealImplementation(_multiDomainIterator).hostIntersectionIterator();
  }

  bool equals(const IntersectionIteratorWrapper& rhs) const {
    return _indexSet == rhs._indexSet && _multiDomainIterator == rhs._multiDomainIterator;
  }

  void increment() {
    ++_multiDomainIterator;
  }

  Intersection dereference() const {
    return {IntersectionWrapper(_indexSet,*_multiDomainIterator)};
  }

private:

  const IndexSet* _indexSet;
  MultiDomainIntersectionIterator _multiDomainIterator;

};


} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_INTERSECTIONITERATOR_HH
