#ifndef DUNE_MULTIDOMAINGRID_INTERSECTIONITERATOR_HH
#define DUNE_MULTIDOMAINGRID_INTERSECTIONITERATOR_HH

namespace Dune {

namespace mdgrid {

template<typename, PartitionIteratorType>
class LevelGridView;

template<typename, PartitionIteratorType>
class LeafGridView;


template<typename GridImp,
	 typename WrapperImp,
	 typename HostIntersectionIteratorType,
	 typename IntersectionType>
class IntersectionIteratorWrapper {

  template<class, class, class>
  friend class Dune::IntersectionIterator;

  template<class, class>
  friend class Dune::Intersection;

  template<typename,typename,typename,typename,typename>
  friend class subdomain::IntersectionIteratorWrapper;

  typedef HostIntersectionIteratorType HostIntersectionIterator;
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

  explicit IntersectionIteratorWrapper(const HostIntersectionIterator& hostIterator) :
    _hostIterator(hostIterator),
    _intersection(typename GridImp::template ReturnImplementationType<IntersectionType>::ImplementationType(NULL))
  {}

  IntersectionIteratorWrapper(const IntersectionIteratorWrapper& rhs)
    : _hostIterator(rhs._hostIterator)
    , _intersection(typename GridImp::template ReturnImplementationType<IntersectionType>::ImplementationType(NULL))
  {}

  IntersectionIteratorWrapper& operator=(const IntersectionIteratorWrapper& rhs)
  {
    _hostIterator = rhs._hostIterator;
    GridImp::getRealImplementation(_intersection).reset(*_hostIterator);
    return *this;
  }

private:

  const HostIntersectionIterator& hostIntersectionIterator() const {
    return _hostIterator;
  }

  bool equals(const WrapperImp& rhs) const {
    return _hostIterator == rhs._hostIterator;
  }

  void increment() {
    GridImp::getRealImplementation(_intersection).clear();
    ++_hostIterator;
  }

  const Intersection& dereference() const {
    if (!GridImp::getRealImplementation(_intersection).isSet()) {
      GridImp::getRealImplementation(_intersection).reset(*_hostIterator);
    }
    return _intersection;
  }

  HostIntersectionIterator _hostIterator;
  mutable Intersection _intersection;

};

template<typename GridImp>
class LeafIntersectionIteratorWrapper :
    public IntersectionIteratorWrapper<GridImp,
				       LeafIntersectionIteratorWrapper<GridImp>,
				       typename detail::HostGridAccessor<GridImp>::Traits::LeafIntersectionIterator,
				       typename GridImp::Traits::LeafIntersection>
{

  template<typename, typename, typename, typename>
  friend class IntersectionIteratorWrapper;

  template<typename,typename,typename,typename,typename>
  friend class subdomain::IntersectionIteratorWrapper;

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename,typename>
  friend class MultiDomainGrid;

  template<typename, PartitionIteratorType>
  friend class LeafGridView;


  typedef typename GridImp::HostGridType::Traits::LeafIntersectionIterator HostIntersectionIterator;
  typedef typename GridImp::Traits::LeafIntersection Intersection;

  typedef IntersectionIteratorWrapper<GridImp,
				      LeafIntersectionIteratorWrapper<GridImp>,
				      HostIntersectionIterator,
				      Intersection> Base;

  explicit LeafIntersectionIteratorWrapper(const HostIntersectionIterator& hostIterator) :
    Base(hostIterator)
  {}

};


template<typename GridImp>
class LevelIntersectionIteratorWrapper :
    public IntersectionIteratorWrapper<GridImp,
				       LevelIntersectionIteratorWrapper<GridImp>,
				       typename detail::HostGridAccessor<GridImp>::Traits::LevelIntersectionIterator,
				       typename GridImp::Traits::LevelIntersection>

{

  template<typename, typename, typename, typename>
  friend class IntersectionIteratorWrapper;

  template<typename,typename,typename,typename,typename>
  friend class subdomain::IntersectionIteratorWrapper;

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename,typename>
  friend class MultiDomainGrid;

  template<typename, PartitionIteratorType>
  friend class LevelGridView;


  typedef typename GridImp::HostGridType::Traits::LevelIntersectionIterator HostIntersectionIterator;
  typedef typename GridImp::Traits::LevelIntersection Intersection;

  typedef IntersectionIteratorWrapper<GridImp,
				      LevelIntersectionIteratorWrapper<GridImp>,
				      HostIntersectionIterator,
				      Intersection> Base;

  explicit LevelIntersectionIteratorWrapper(const HostIntersectionIterator& hostIterator) :
    Base(hostIterator)
  {}

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_INTERSECTIONITERATOR_HH
