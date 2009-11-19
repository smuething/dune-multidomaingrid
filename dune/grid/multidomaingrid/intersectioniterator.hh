#ifndef DUNE_MULTIDOMAINGRID_INTERSECTIONITERATOR_HH
#define DUNE_MULTIDOMAINGRID_INTERSECTIONITERATOR_HH

namespace Dune {

namespace mdgrid {

template<typename GridImp,
	 typename WrapperImp,
	 typename HostIntersectionIteratorType,
	 typename IntersectionType>
class IntersectionIteratorWrapper {

  template<class, template<class> class, template<class> class>
  friend class IntersectionIterator;

  template<class, template<class> class>
  friend class Intersection;

  template<typename,typename,typename,typename,typename>
  friend class subdomain::IntersectionIteratorWrapper;

  typedef HostIntersectionIteratorType HostIntersectionIterator;
  typedef IntersectionType Intersection;

  typedef typename GridImp::Traits::template Codim<0>::EntityPointer EntityPointer;
  typedef typename GridImp::Traits::template Codim<0>::Entity Entity;
  typedef typename GridImp::Traits::template Codim<1>::Geometry Geometry;
  typedef typename GridImp::Traits::template Codim<1>::LocalGeometry LocalGeometry;

  typedef typename
GridImp::ctype ctype;
  static const int dimension = GridImp::dimension;
  static const int dimensionworld = GridImp::dimensionworld;

  typedef FieldVector<ctype,dimensionworld> GlobalCoords;
  typedef FieldVector<ctype,dimension - 1> LocalCoords;

protected:

  explicit IntersectionIteratorWrapper(const HostIntersectionIterator& hostIterator) :
    _hostIterator(hostIterator)
  {}

private:

  const HostIntersectionIterator& hostIntersectionIterator() const {
    return _hostIterator;
  }

  bool equals(const WrapperImp& rhs) const {
    return _hostIterator == rhs._hostIterator;
  }

  void increment() {
    ++_hostIterator;
    _geometry.clear();
    _geometryInInside.clear();
    _geometryInOutside.clear();
  }

  const Intersection& dereference() const {
    return reinterpret_cast<const Intersection&>(*this);
  }

  bool boundary() const {
    return _hostIterator->boundary();
  }

  int boundaryId() const {
    return _hostIterator->boundaryId();
  }

  bool neighbor() const {
    return _hostIterator->neighbor();
  }

  EntityPointer inside() const {
    return EntityPointerWrapper<0,GridImp>(_hostIterator->inside());
  }

  EntityPointer outside() const {
    return EntityPointerWrapper<0,GridImp>(_hostIterator->outside());
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
    return _hostIterator->indexInOutside();
  }

  int numberInNeighbor() const {
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
  HostIntersectionIterator _hostIterator;
  MakeableGeometryWrapper<LocalGeometry::mydimension,LocalGeometry::coorddimension,GridImp> _geometryInInside, _geometryInOutside;
MakeableGeometryWrapper<Geometry::mydimension,Geometry::coorddimension,GridImp> _geometry;

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
