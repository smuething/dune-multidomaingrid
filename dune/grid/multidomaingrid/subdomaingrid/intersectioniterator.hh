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

  template<class, template<class> class, template<class> class>
  friend class IntersectionIterator;

  template<class, template<class> class>
  friend class Intersection;

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
    _indexSet(indexSet),
    _multiDomainIterator(multiDomainIterator),
    _outsideTested(false)
  {}

  const IntersectionIteratorWrapper& operator=(const IntersectionIteratorWrapper& rhs) {
    assert(_indexSet == rhs._indexSet);
    _multiDomainIterator = rhs._multiDomainIterator;
    _outsideTested = false;
    _geometry.clear();
    _geometryInInside.clear();
    _geometryInOutside.clear();
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
    ++_multiDomainIterator;
    _geometry.clear();
    _geometryInInside.clear();
    _geometryInOutside.clear();
    _outsideTested = false;
  }

  const Intersection& dereference() const {
    return reinterpret_cast<const Intersection&>(*this);
  }

  void checkOutside() const {
    if (!_outsideTested) {
      if (_multiDomainIterator->boundary()) {
        _outsideType = otBoundary;
      } else {
        if (_indexSet.containsMultiDomainEntity(*(_multiDomainIterator->outside()))) {
          _outsideType = otNeighbor;
        } else {
          _outsideType = otForeignCell;
        }
      }
    }
  }

  bool boundary() const {
    checkOutside();
    return _outsideType != otNeighbor;
  }

  int boundaryId() const {
    return _multiDomainIterator->boundaryId();
  }

  bool neighbor() const {
    checkOutside();
    return _outsideType == otNeighbor;
  }

  EntityPointer inside() const {
    return EntityPointerWrapper<0,GridImp>(_indexSet._grid,_multiDomainIterator->inside());
  }

  EntityPointer outside() const {
    checkOutside();
    assert(_outsideType == otNeighbor);
    return EntityPointerWrapper<0,GridImp>(_indexSet._grid,_multiDomainIterator->outside());
  }

  bool conforming() const {
    return _multiDomainIterator->conforming();
  }

  const LocalGeometry& geometryInInside() const {
    if (!_geometryInInside.isSet()) {
      _geometryInInside.reset((*hostIntersectionIterator()).geometryInInside());
    }
    return _geometryInInside;
  }

  const LocalGeometry& intersectionSelfLocal() const {
    return geometryInInside();
  }

  const LocalGeometry& geometryInOutside() const {
    checkOutside();
    assert(_outsideType == otNeighbor);
    if (!_geometryInOutside.isSet()) {
      _geometryInOutside.reset((*hostIntersectionIterator()).geometryInOutside());
    }
    return _geometryInOutside;
  }

  const LocalGeometry& intersectionNeighborLocal() const {
    return geometryInOutside();
  }

  const Geometry& geometry() const {
    if (!_geometry.isSet()) {
      _geometry.reset((*hostIntersectionIterator()).geometry());
    }
    return _geometry;
  }

  const LocalGeometry& intersectionGlobal() const {
    return geometry();
  }

  GeometryType type() const {
    return _multiDomainIterator->type();
  }

  int indexInInside() const {
    return _multiDomainIterator->indexInInside();
  }

  int numberInSelf() const {
    return _multiDomainIterator->numberInSelf();
  }

  int indexInOutside() const {
    checkOutside();
    assert(_outsideType == otNeighbor);
    return _multiDomainIterator->indexInOutside();
  }

  int numberInNeighbor() const {
    checkOutside();
    assert(_outsideType == otNeighbor);
    return _multiDomainIterator->numberInNeighbor();
  }

  GlobalCoords outerNormal(const LocalCoords& local) const {
    return _multiDomainIterator->outerNormal(local);
  }

  GlobalCoords integrationOuterNormal(const LocalCoords& local) const {
    return _multiDomainIterator->integrationOuterNormal(local);
  }

  GlobalCoords unitOuterNormal(const LocalCoords& local) const {
    return _multiDomainIterator->unitOuterNormal(local);
  }

private:

  enum OutsideType { otNeighbor, otForeignCell, otBoundary };

  const IndexSet& _indexSet;
  MultiDomainIntersectionIterator _multiDomainIterator;
  MakeableGeometryWrapper<LocalGeometry::mydimension,LocalGeometry::coorddimension,GridImp> _geometryInInside, _geometryInOutside;
  MakeableGeometryWrapper<Geometry::mydimension,Geometry::coorddimension,GridImp> _geometry;
  mutable bool _outsideTested;
  mutable OutsideType _outsideType;


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
