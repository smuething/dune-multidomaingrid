#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_INTERSECTION_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_INTERSECTION_HH

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<typename GridImp,
	 typename WrapperImp,
         typename IndexSet,
	 typename MultiDomainIntersectionType>
class IntersectionWrapper {

  template<class, class, class>
  friend class Dune::IntersectionIterator;

  template<class,class,class,class,class>
  friend class IntersectionIteratorWrapper;

  template<class, class>
  friend class Dune::Intersection;

  template<typename MDGrid>
  friend class SubDomainGrid;

  typedef MultiDomainIntersectionType MultiDomainIntersection;

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

  IntersectionWrapper(const IndexSet& indexSet, const MultiDomainIntersection* multiDomainIntersection) :
    _indexSet(indexSet),
    _multiDomainIntersection(multiDomainIntersection),
    _intersectionTypeTested(false)
  {}

  // copy constructor is required to make sure the wrapper does not pick up a pointer to the wrapped
  // intersection of a foreign IntersectionIterator
  IntersectionWrapper(const IntersectionWrapper& rhs) :
    _indexSet(rhs._indexSet),
    _multiDomainIntersection(NULL),
    _intersectionTypeTested(false)
  {}

  const IntersectionWrapper& operator=(const IntersectionWrapper& rhs);

private:

  const typename GridImp::MultiDomainGrid::template ReturnImplementationType<MultiDomainIntersection>::ImplementationType::HostIntersection& hostIntersection() const {
    return GridImp::MultiDomainGrid::getRealImplementation(*_multiDomainIntersection).hostIntersection();
  }


  bool equals(const WrapperImp& rhs) const {
    return *_multiDomainIntersection == *(rhs._multiDomainIntersection);
  }

  void checkIntersectionType() const {
    if (!_intersectionTypeTested) {
      if (_multiDomainIntersection->boundary()) {
        _intersectionType = GridImp::boundary;
        _intersectionTypeTested = true;
        return;
      }
      if (!_multiDomainIntersection->neighbor()) {
        _intersectionType = GridImp::processor;
        _intersectionTypeTested = true;
        return;
      }
      if (_indexSet.containsMultiDomainEntity(*(_multiDomainIntersection->outside()))) {
        _intersectionType = GridImp::neighbor;
        _intersectionTypeTested = true;
        return;
      } else {
        _intersectionType = GridImp::foreign;
        _intersectionTypeTested = true;
        return;
      }
      assert(false && "Should never get here - invalid intersection type!");
    }
  }

  bool boundary() const {
    checkIntersectionType();
    return
      _intersectionType == GridImp::boundary ||
      _intersectionType == GridImp::foreign;
  }

  int boundaryId() const {
    return _multiDomainIntersection->boundaryId();
  }

  std::size_t boundarySegmentIndex() const {
    return _multiDomainIntersection->boundarySegmentIndex();
  }

  bool neighbor() const {
    checkIntersectionType();
    return _intersectionType == GridImp::neighbor;
  }

  EntityPointer inside() const {
    return EntityPointerWrapper<0,GridImp>(_indexSet._grid,_multiDomainIntersection->inside());
  }

  EntityPointer outside() const {
    checkIntersectionType();
    assert(_intersectionType == GridImp::neighbor);
    return EntityPointerWrapper<0,GridImp>(_indexSet._grid,_multiDomainIntersection->outside());
  }

  bool conforming() const {
    return _multiDomainIntersection->conforming();
  }

  LocalGeometry geometryInInside() const {
    return LocalGeometry(hostIntersection().geometryInInside());
  }

  LocalGeometry geometryInOutside() const {
    checkIntersectionType();
    assert(_intersectionType == GridImp::neighbor);
    return LocalGeometry(hostIntersection().geometryInOutside());
  }

  Geometry geometry() const {
    return Geometry(hostIntersection().geometry());
  }

  GeometryType type() const {
    return _multiDomainIntersection->type();
  }

  int indexInInside() const {
    return _multiDomainIntersection->indexInInside();
  }

  int indexInOutside() const {
    checkIntersectionType();
    assert(_intersectionType == GridImp::neighbor);
    return _multiDomainIntersection->indexInOutside();
  }

  GlobalCoords outerNormal(const LocalCoords& local) const {
    return _multiDomainIntersection->outerNormal(local);
  }

  GlobalCoords integrationOuterNormal(const LocalCoords& local) const {
    return _multiDomainIntersection->integrationOuterNormal(local);
  }

  GlobalCoords unitOuterNormal(const LocalCoords& local) const {
    return _multiDomainIntersection->unitOuterNormal(local);
  }

  GlobalCoords centerUnitOuterNormal() const {
    return _multiDomainIntersection->centerUnitOuterNormal();
  }

  typename GridImp::IntersectionType intersectionType() const {
    checkIntersectionType();
    return _intersectionType;
  }

  const MultiDomainIntersection& multiDomainIntersection() const
  {
    return *_multiDomainIntersection;
  }

private:

  bool isSet() const {
    return _multiDomainIntersection != NULL;
  }

  void clear() {
    _multiDomainIntersection = NULL;
    _intersectionTypeTested = false;
  }

  void reset(const MultiDomainIntersection& multiDomainIntersection) {
    if (isSet()) {
      clear();
    }
    _multiDomainIntersection = &multiDomainIntersection;
  }

  const IndexSet& _indexSet;
  const MultiDomainIntersection* _multiDomainIntersection;
  mutable bool _intersectionTypeTested;
  mutable typename GridImp::IntersectionType _intersectionType;

};

template<typename GridImp>
class LeafIntersectionWrapper :
    public IntersectionWrapper<GridImp,
                               LeafIntersectionWrapper<GridImp>,
                               typename GridImp::Traits::LeafIndexSet,
                               typename GridImp::MultiDomainGrid::Traits::LeafIntersection>
{

  template<typename, typename, typename, typename, typename>
  friend class IntersectionIteratorWrapper;

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename>
  friend class SubDomainGrid;

  typedef typename GridImp::MultiDomainGrid::Traits::LeafIntersection MultiDomainIntersection;
  typedef typename GridImp::Traits::LeafIndexSet IndexSet;

  typedef IntersectionWrapper<GridImp,
                              LeafIntersectionWrapper<GridImp>,
                              IndexSet,
                              MultiDomainIntersection> Base;

  LeafIntersectionWrapper(const IndexSet& indexSet, const MultiDomainIntersection* multiDomainIntersection) :
    Base(indexSet,multiDomainIntersection)
  {}

};


template<typename GridImp>
class LevelIntersectionWrapper :
    public IntersectionWrapper<GridImp,
                               LevelIntersectionWrapper<GridImp>,
                               typename GridImp::Traits::LevelIndexSet,
                               typename GridImp::MultiDomainGrid::Traits::LevelIntersection>
{

  template<typename, typename, typename, typename, typename>
  friend class IntersectionIteratorWrapper;

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename>
  friend class SubDomainGrid;

  typedef typename GridImp::MultiDomainGrid::Traits::LevelIntersection MultiDomainIntersection;
  typedef typename GridImp::Traits::LevelIndexSet IndexSet;

  typedef IntersectionWrapper<GridImp,
                              LevelIntersectionWrapper<GridImp>,
                              IndexSet,
                              MultiDomainIntersection> Base;

  LevelIntersectionWrapper(const IndexSet& indexSet, const MultiDomainIntersection* multiDomainIntersection) :
    Base(indexSet,multiDomainIntersection)
  {}

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_INTERSECTION_HH
