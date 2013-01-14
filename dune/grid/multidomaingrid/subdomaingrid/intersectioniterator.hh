#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_INTERSECTIONITERATOR_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_INTERSECTIONITERATOR_HH

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<typename GridImp,
	 typename WrapperImp,
         typename IndexSet,
	 typename MultiDomainIntersectionIteratorType,
	 typename IntersectionType>
class IntersectionIteratorWrapper {

  template<class, class, class>
  friend class Dune::IntersectionIterator;

  template<class, class>
  friend class Dune::Intersection;

  typedef MultiDomainIntersectionIteratorType MultiDomainIntersectionIterator;
  typedef IntersectionType Intersection;

  typedef typename GridImp::Traits::template Codim<0>::EntityPointer EntityPointer;
  typedef typename GridImp::Traits::template Codim<0>::Entity Entity;
  typedef typename GridImp::Traits::template Codim<1>::Geometry Geometry;
  typedef typename GridImp::Traits::template Codim<1>::LocalGeometry LocalGeometry;

  typedef typename GridImp::ctype ctype;
  static const int dimension = GridImp::dimension;
  static const int dimensionworld = GridImp::dimensionworld;

  typedef FieldVector<ctype,dimensionworld> GlobalCoords;
  typedef FieldVector<ctype,dimension - 1> LocalCoords;

protected:

  IntersectionIteratorWrapper(const IndexSet& indexSet, const MultiDomainIntersectionIterator& multiDomainIterator) :
    _multiDomainIterator(multiDomainIterator),
    _intersection(typename GridImp::template ReturnImplementationType<IntersectionType>::ImplementationType(indexSet,NULL))
  {}

  IntersectionIteratorWrapper(const IntersectionIteratorWrapper& rhs)
    : _multiDomainIterator(rhs._multiDomainIterator)
    , _intersection(GridImp::getRealImplementation(rhs._intersection))
  {}

  const IntersectionIteratorWrapper& operator=(const IntersectionIteratorWrapper& rhs) {
    assert(GridImp::getRealImplementation(_intersection)._indexSet == GridImp::getRealImplementation(rhs._intersection)._indexSet);
    _multiDomainIterator = rhs._multiDomainIterator;
    GridImp::getRealImplementation(_intersection).clear();
    return *this;
  }

private:

  const typename GridImp::MultiDomainGrid::template ReturnImplementationType<MultiDomainIntersectionIterator>::ImplementationType::HostIntersectionIterator& hostIntersectionIterator() const {
    return GridImp::MultiDomainGrid::getRealImplementation(_multiDomainIterator).hostIntersectionIterator();
  }

  bool equals(const WrapperImp& rhs) const {
    return _multiDomainIterator == rhs._multiDomainIterator;
  }

  void increment() {
    GridImp::getRealImplementation(_intersection).clear();
    ++_multiDomainIterator;
  }

  const Intersection& dereference() const {
    if (!GridImp::getRealImplementation(_intersection).isSet()) {
      GridImp::getRealImplementation(_intersection).reset(*_multiDomainIterator);
    }
    return _intersection;
  }

private:

  MultiDomainIntersectionIterator _multiDomainIterator;
  mutable IntersectionType _intersection;

};

template<typename GridImp>
class LeafIntersectionIteratorWrapper :
    public IntersectionIteratorWrapper<GridImp,
				       LeafIntersectionIteratorWrapper<GridImp>,
                                       typename GridImp::Traits::LeafIndexSet,
				       typename GridImp::MultiDomainGrid::Traits::LeafIntersectionIterator,
				       typename GridImp::Traits::LeafIntersection>
{

  template<typename, typename, typename, typename, typename>
  friend class IntersectionIteratorWrapper;

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename>
  friend class SubDomainGrid;

  typedef typename GridImp::MultiDomainGrid::Traits::LeafIntersectionIterator MultiDomainIntersectionIterator;
  typedef typename GridImp::Traits::LeafIntersection Intersection;

  typedef IntersectionIteratorWrapper<GridImp,
				      LeafIntersectionIteratorWrapper<GridImp>,
                                      typename GridImp::Traits::LeafIndexSet,
				      MultiDomainIntersectionIterator,
				      Intersection> Base;

  LeafIntersectionIteratorWrapper(const GridImp& grid, const MultiDomainIntersectionIterator& multiDomainIterator) :
    Base(grid.leafIndexSet(),multiDomainIterator)
  {}

};


template<typename GridImp>
class LevelIntersectionIteratorWrapper :
    public IntersectionIteratorWrapper<GridImp,
				       LevelIntersectionIteratorWrapper<GridImp>,
                                       typename GridImp::Traits::LevelIndexSet,
				       typename GridImp::MultiDomainGrid::Traits::LevelIntersectionIterator,
				       typename GridImp::Traits::LevelIntersection>

{

  template<typename, typename, typename, typename, typename>
  friend class IntersectionIteratorWrapper;

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename>
  friend class SubDomainGrid;

  typedef typename GridImp::MultiDomainGrid::Traits::LevelIntersectionIterator MultiDomainIntersectionIterator;
  typedef typename GridImp::Traits::LevelIntersection Intersection;

  typedef IntersectionIteratorWrapper<GridImp,
				      LevelIntersectionIteratorWrapper<GridImp>,
                                      typename GridImp::Traits::LevelIndexSet,
				      MultiDomainIntersectionIterator,
				      Intersection> Base;

  LevelIntersectionIteratorWrapper(const GridImp& grid, int level, const MultiDomainIntersectionIterator& multiDomainIterator) :
    Base(grid.levelIndexSet(level),multiDomainIterator)
  {}

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_INTERSECTIONITERATOR_HH
