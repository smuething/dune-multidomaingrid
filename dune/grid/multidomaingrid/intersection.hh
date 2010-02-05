#ifndef DUNE_MULTIDOMAINGRID_INTERSECTION_HH
#define DUNE_MULTIDOMAINGRID_INTERSECTION_HH

namespace Dune {

namespace mdgrid {

template<typename GridImp,
	 typename WrapperImp,
	 typename HostIntersectionType>
class IntersectionWrapper {

  template<class, template<class> class>
  friend class Intersection;

  template<typename,typename,typename,typename,typename>
  friend class subdomain::IntersectionIteratorWrapper;

  typedef HostIntersectionType HostIntersection;

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

  explicit IntersectionWrapper(const HostIntersection* hostIntersection) :
    _hostIntersection(hostIntersection)
  {}

private:

  const HostIntersection& hostIntersection() const {
    return *_hostIntersection;
  }

  bool equals(const WrapperImp& rhs) const {
    return *_hostIntersection == *(rhs._hostIntersection);
  }

  bool boundary() const {
    return _hostIntersection->boundary();
  }

  int boundaryId() const {
    return _hostIntersection->boundaryId();
  }

  std::size_t boundarySegmentIndex() const {
    return _hostIntersection->boundarySegmentIndex();
  }

  bool neighbor() const {
    return _hostIntersection->neighbor();
  }

  EntityPointer inside() const {
    return EntityPointerWrapper<0,GridImp>(_hostIntersection->inside());
  }

  EntityPointer outside() const {
    return EntityPointerWrapper<0,GridImp>(_hostIntersection->outside());
  }

  bool conforming() const {
    return _hostIntersection->conforming();
  }

  const LocalGeometry& geometryInInside() const {
    if (!_geometryInInside.isSet()) {
      _geometryInInside.reset(_hostIntersection->geometryInInside());
    }
    return _geometryInInside;
  }

  const LocalGeometry& intersectionSelfLocal() const {
    return geometryInInside();
  }

  const LocalGeometry& geometryInOutside() const {
    if (!_geometryInOutside.isSet()) {
      _geometryInOutside.reset(_hostIntersection->geometryInOutside());
    }
    return _geometryInOutside;
  }

  const LocalGeometry& intersectionNeighborLocal() const {
    return geometryInOutside();
  }

  const Geometry& geometry() const {
    if (!_geometry.isSet()) {
      _geometry.reset(_hostIntersection->geometry());
    }
    return _geometry;
  }

  const LocalGeometry& intersectionGlobal() const {
    return geometry();
  }

  GeometryType type() const {
    return _hostIntersection->type();
  }

  int indexInInside() const {
    return _hostIntersection->indexInInside();
  }

  int numberInSelf() const {
    return _hostIntersection->numberInSelf();
  }

  int indexInOutside() const {
    return _hostIntersection->indexInOutside();
  }

  int numberInNeighbor() const {
    return _hostIntersection->numberInNeighbor();
  }

  GlobalCoords outerNormal(const LocalCoords& local) const {
    return _hostIntersection->outerNormal(local);
  }

  GlobalCoords integrationOuterNormal(const LocalCoords& local) const {
    return _hostIntersection->integrationOuterNormal(local);
  }

  GlobalCoords unitOuterNormal(const LocalCoords& local) const {
    return _hostIntersection->unitOuterNormal(local);
  }

  GlobalCoords centerUnitOuterNormal() const {
    return _hostIntersection->centerUnitOuterNormal();
  }

private:

  bool isSet() const {
    return _hostIntersection != NULL;
  }

  void clear() {
    _hostIntersection = NULL;
    _geometry.clear();
    _geometryInInside.clear();
    _geometryInOutside.clear();
  }

  void reset(const HostIntersection& hostIntersection) {
    if (isSet()) {
      clear();
    }
    _hostIntersection = &hostIntersection;
  }

  const HostIntersection* _hostIntersection;
  MakeableGeometryWrapper<LocalGeometry::mydimension,LocalGeometry::coorddimension,GridImp> _geometryInInside, _geometryInOutside;
  MakeableGeometryWrapper<Geometry::mydimension,Geometry::coorddimension,GridImp> _geometry;

};

template<typename GridImp>
class LeafIntersectionWrapper :
    public IntersectionWrapper<GridImp,
                               LeafIntersectionWrapper<GridImp>,
                               typename detail::HostGridAccessor<GridImp>::Traits::LeafIntersection>
{

  template<typename, typename, typename>
  friend class IntersectionWrapper;

  template<typename, typename, typename, typename>
  friend class IntersectionIteratorWrapper;

  template<typename,typename,typename,typename>
  friend class subdomain::IntersectionWrapper;

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename,typename>
  friend class MultiDomainGrid;

  typedef typename GridImp::HostGridType::Traits::LeafIntersection HostIntersection;
  typedef typename GridImp::Traits::LeafIntersection Intersection;

  typedef IntersectionWrapper<GridImp,
                              LeafIntersectionWrapper<GridImp>,
                              HostIntersection> Base;

  explicit LeafIntersectionWrapper(const HostIntersection* hostIntersection) :
    Base(hostIntersection)
  {}

};


template<typename GridImp>
class LevelIntersectionWrapper :
    public IntersectionWrapper<GridImp,
                               LevelIntersectionWrapper<GridImp>,
                               typename detail::HostGridAccessor<GridImp>::Traits::LevelIntersection>
{

  template<typename, typename, typename>
  friend class IntersectionWrapper;

  template<typename, typename, typename, typename>
  friend class IntersectionIteratorWrapper;

  template<typename,typename,typename,typename>
  friend class subdomain::IntersectionWrapper;

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename,typename>
  friend class MultiDomainGrid;

  typedef typename GridImp::HostGridType::Traits::LevelIntersection HostIntersection;
  typedef typename GridImp::Traits::LevelIntersection Intersection;

  typedef IntersectionWrapper<GridImp,
                              LevelIntersectionWrapper<GridImp>,
                              HostIntersection> Base;

  explicit LevelIntersectionWrapper(const HostIntersection* hostIntersection) :
    Base(hostIntersection)
  {}

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_INTERSECTION_HH
