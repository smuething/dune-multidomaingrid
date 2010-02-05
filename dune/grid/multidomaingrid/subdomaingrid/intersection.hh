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

  template<class, template<class> class, template<class> class>
  friend class IntersectionIterator;

  template<class,class,class,class,class>
  friend class IntersectionIteratorWrapper;

  template<class, template<class> class>
  friend class Intersection;

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
    _outsideTested(false)
  {}

  const IntersectionWrapper& operator=(const IntersectionWrapper& rhs);

private:

  const typename GridImp::MultiDomainGrid::template ReturnImplementationType<MultiDomainIntersection>::ImplementationType::HostIntersection& hostIntersection() const {
    return GridImp::MultiDomainGrid::getRealImplementation(*_multiDomainIntersection).hostIntersection();
  }


  bool equals(const WrapperImp& rhs) const {
    return *_multiDomainIntersection == *(rhs._multiDomainIntersection);
  }

  void checkOutside() const {
    if (!_outsideTested) {
      if (_multiDomainIntersection->boundary()) {
        _outsideType = otBoundary;
      } else {
        if (_indexSet.containsMultiDomainEntity(*(_multiDomainIntersection->outside()))) {
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
    return _multiDomainIntersection->boundaryId();
  }

  std::size_t boundarySegmentIndex() const {
    return _multiDomainIntersection->boundarySegmentIndex();
  }

  bool neighbor() const {
    checkOutside();
    return _outsideType == otNeighbor;
  }

  EntityPointer inside() const {
    return EntityPointerWrapper<0,GridImp>(_indexSet._grid,_multiDomainIntersection->inside());
  }

  EntityPointer outside() const {
    checkOutside();
    assert(_outsideType == otNeighbor);
    return EntityPointerWrapper<0,GridImp>(_indexSet._grid,_multiDomainIntersection->outside());
  }

  bool conforming() const {
    return _multiDomainIntersection->conforming();
  }

  const LocalGeometry& geometryInInside() const {
    if (!_geometryInInside.isSet()) {
      _geometryInInside.reset(hostIntersection().geometryInInside());
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
      _geometryInOutside.reset(hostIntersection().geometryInOutside());
    }
    return _geometryInOutside;
  }

  const LocalGeometry& intersectionNeighborLocal() const {
    return geometryInOutside();
  }

  const Geometry& geometry() const {
    if (!_geometry.isSet()) {
      _geometry.reset(hostIntersection().geometry());
    }
    return _geometry;
  }

  const LocalGeometry& intersectionGlobal() const {
    return geometry();
  }

  GeometryType type() const {
    return _multiDomainIntersection->type();
  }

  int indexInInside() const {
    return _multiDomainIntersection->indexInInside();
  }

  int numberInSelf() const {
    return _multiDomainIntersection->numberInSelf();
  }

  int indexInOutside() const {
    checkOutside();
    assert(_outsideType == otNeighbor);
    return _multiDomainIntersection->indexInOutside();
  }

  int numberInNeighbor() const {
    checkOutside();
    assert(_outsideType == otNeighbor);
    return _multiDomainIntersection->numberInNeighbor();
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

private:

  bool isSet() const {
    return _multiDomainIntersection != NULL;
  }

  void clear() {
    _multiDomainIntersection = NULL;
    _geometryInInside.clear();
    _geometryInOutside.clear();
    _geometry.clear();
    _outsideTested = false;
  }

  void reset(const MultiDomainIntersection& multiDomainIntersection) {
    if (isSet()) {
      clear();
    }
    _multiDomainIntersection = &multiDomainIntersection;
  }

  enum OutsideType { otNeighbor, otForeignCell, otBoundary };

  const IndexSet& _indexSet;
  const MultiDomainIntersection* _multiDomainIntersection;
  MakeableGeometryWrapper<LocalGeometry::mydimension,LocalGeometry::coorddimension,GridImp> _geometryInInside, _geometryInOutside;
  MakeableGeometryWrapper<Geometry::mydimension,Geometry::coorddimension,GridImp> _geometry;
  mutable bool _outsideTested;
  mutable OutsideType _outsideType;

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
