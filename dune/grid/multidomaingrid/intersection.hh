#ifndef DUNE_MULTIDOMAINGRID_INTERSECTION_HH
#define DUNE_MULTIDOMAINGRID_INTERSECTION_HH

#include <dune/grid/common/intersection.hh>

#include <dune/grid/multidomaingrid/entitypointer.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/intersection.hh>

namespace Dune {

namespace mdgrid {

template<typename GridImp,
	 typename WrapperImp,
	 typename HostIntersectionType>
class IntersectionWrapper {

  template<class, class>
  friend class Dune::Intersection;

  template<typename,typename,typename,typename>
  friend class subdomain::IntersectionWrapper;

  template<typename,typename,typename,typename>
  friend class IntersectionIteratorWrapper;

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

  // copy constructor is required to make sure the wrapper does not pick up a pointer to the wrapped
  // intersection of a foreign IntersectionIterator
  IntersectionWrapper(const IntersectionWrapper& rhs) :
    _hostIntersection(NULL)
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

  LocalGeometry geometryInInside() const {
    return LocalGeometry(_hostIntersection->geometryInInside());
  }

  LocalGeometry geometryInOutside() const {
    return LocalGeometry(_hostIntersection->geometryInOutside());
  }

  Geometry geometry() const {
    return Geometry(_hostIntersection->geometry());
  }

  GeometryType type() const {
    return _hostIntersection->type();
  }

  int indexInInside() const {
    return _hostIntersection->indexInInside();
  }

  int indexInOutside() const {
    return _hostIntersection->indexInOutside();
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
  }

  void reset(const HostIntersection& hostIntersection) {
    if (isSet()) {
      clear();
    }
    _hostIntersection = &hostIntersection;
  }

  const HostIntersection* _hostIntersection;

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

  template<typename, typename>
  friend class Dune::Intersection;

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

  LeafIntersectionWrapper(const LeafIntersectionWrapper& rhs):
    Base(rhs)
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

  template<typename, typename>
  friend class Dune::Intersection;

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

  LevelIntersectionWrapper(const LevelIntersectionWrapper& rhs):
    Base(rhs)
  {}

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_INTERSECTION_HH
