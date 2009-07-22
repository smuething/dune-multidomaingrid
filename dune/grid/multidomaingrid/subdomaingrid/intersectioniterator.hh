#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_INTERSECTIONITERATOR_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_INTERSECTIONITERATOR_HH

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<typename GridImp,
	 typename WrapperImp,
         typename IndexSet,
	 typename HostIntersectionIteratorType,
	 typename IntersectionType>
class IntersectionIteratorWrapper {

  template<class, template<class> class, template<class> class>
  friend class IntersectionIterator;

  template<class, template<class> class>
  friend class Intersection;

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

  IntersectionIteratorWrapper(const IndexSet& indexSet, const HostIntersectionIterator& hostIterator) :
    _indexSet(indexSet),
    _hostIterator(hostIterator),
    _outsideTested(false)
  {}

  const IntersectionIteratorWrapper& operator=(const IntersectionIteratorWrapper& rhs) {
    assert(_indexSet == rhs._indexSet);
    _hostIterator = rhs._hostIterator;
    _outsideTested = false;
    _geometry.clear();
    _geometryInInside.clear();
    _geometryInOutside.clear();
  }

private:

  bool equals(const WrapperImp& rhs) const {
    return _hostIterator == rhs._hostIterator;
  }

  void increment() {
    ++_hostIterator;
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
      if (_hostIterator->boundary()) {
        _outsideType = otBoundary;
      } else {
        if (_indexSet.containsHostEntity(*(_hostIterator->outside()))) {
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
    return _hostIterator->boundaryId();
  }

  bool neighbor() const {
    checkOutside();
    return _outsideType == otNeighbor;
  }

  EntityPointer inside() const {
    return EntityPointerWrapper<0,GridImp>(_indexSet._grid,_hostIterator->inside());
  }

  EntityPointer outside() const {
    checkOutside();
    assert(_outsideType == otNeighbor);
    return EntityPointerWrapper<0,GridImp>(_indexSet._grid,_hostIterator->outside());
  }

  bool conforming() const {
    return _hostIterator->conforming();
  }

  const LocalGeometry& geometryInInside() const {
    if (!_geometryInInside.isSet()) {
      _geometryInInside.reset(_hostIterator->geometryInInside());
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
      _geometryInOutside.reset(_hostIterator->geometryInOutside());
    }
    return _geometryInOutside;
  }

  const LocalGeometry& intersectionNeighborLocal() const {
    return geometryInOutside();
  }

  const Geometry& geometry() const {
    if (!_geometry.isSet()) {
      _geometry.reset(_hostIterator->geometry());
    }
    return _geometry;
  }

  const LocalGeometry& intersectionGlobal() const {
    return geometry();
  }

  GeometryType type() const {
    return _hostIterator->type();
  }

  int indexInInside() const {
    return _hostIterator->indexInInside();
  }

  int numberInSelf() const {
    return _hostIterator->numberInSelf();
  }

  int indexInOutside() const {
    checkOutside();
    assert(_outsideType == otNeighbor);
    return _hostIterator->indexInOutside();
  }

  int numberInNeighbor() const {
    checkOutside();
    assert(_outsideType == otNeighbor);
    return _hostIterator->numberInNeighbor();
  }

  GlobalCoords outerNormal(const LocalCoords& local) const {
    return _hostIterator->outerNormal(local);
  }

  GlobalCoords integrationOuterNormal(const LocalCoords& local) const {
    return _hostIterator->integrationOuterNormal(local);
  }

  GlobalCoords unitOuterNormal(const LocalCoords& local) const {
    return _hostIterator->unitOuterNormal(local);
  }

private:

  enum OutsideType { otNeighbor, otForeignCell, otBoundary };

  const IndexSet& _indexSet;
  HostIntersectionIterator _hostIterator;
  MakeableGeometryWrapper<LocalGeometry::mydimension,LocalGeometry::coorddimension,GridImp> _geometryInInside, _geometryInOutside;
  MakeableGeometryWrapper<Geometry::mydimension,Geometry::coorddimension,GridImp> _geometry;
  mutable bool _outsideTested;
  mutable OutsideType _outsideType;


};

namespace detail {

  template<typename GridImp>
  struct HostGridTraits {
    typedef typename GridImp::HostGridType::Traits Traits;
  };

}

template<typename GridImp>
class LeafIntersectionIteratorWrapper :
    public IntersectionIteratorWrapper<GridImp,
				       LeafIntersectionIteratorWrapper<GridImp>,
                                       typename GridImp::Traits::LeafIndexSet,
				       typename detail::HostGridTraits<GridImp>::Traits::LeafIntersectionIterator,
				       typename GridImp::Traits::LeafIntersection>
{

  template<typename, typename, typename, typename, typename>
  friend class IntersectionIteratorWrapper;

  template<int, int, typename>
  friend class EntityWrapper;

  typedef typename GridImp::HostGridType::Traits::LeafIntersectionIterator HostIntersectionIterator;
  typedef typename GridImp::Traits::LeafIntersection Intersection;

  typedef IntersectionIteratorWrapper<GridImp,
				      LeafIntersectionIteratorWrapper<GridImp>,
                                      typename GridImp::Traits::LeafIndexSet,
				      HostIntersectionIterator,
				      Intersection> Base;

  LeafIntersectionIteratorWrapper(const GridImp& grid, const HostIntersectionIterator& hostIterator) :
    Base(grid.leafIndexSet(),hostIterator)
  {}

};


template<typename GridImp>
class LevelIntersectionIteratorWrapper :
    public IntersectionIteratorWrapper<GridImp,
				       LevelIntersectionIteratorWrapper<GridImp>,
                                       typename GridImp::Traits::LevelIndexSet,
				       typename detail::HostGridTraits<GridImp>::Traits::LevelIntersectionIterator,
				       typename GridImp::Traits::LevelIntersection>

{

  template<typename, typename, typename, typename, typename>
  friend class IntersectionIteratorWrapper;

  template<int, int, typename>
  friend class EntityWrapper;

  typedef typename GridImp::HostGridType::Traits::LevelIntersectionIterator HostIntersectionIterator;
  typedef typename GridImp::Traits::LevelIntersection Intersection;

  typedef IntersectionIteratorWrapper<GridImp,
				      LevelIntersectionIteratorWrapper<GridImp>,
                                      typename GridImp::Traits::LevelIndexSet,
				      HostIntersectionIterator,
				      Intersection> Base;

  LevelIntersectionIteratorWrapper(const GridImp& grid, int level, const HostIntersectionIterator& hostIterator) :
    Base(grid.levelIndexSet(level),hostIterator)
  {}

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_INTERSECTIONITERATOR_HH
