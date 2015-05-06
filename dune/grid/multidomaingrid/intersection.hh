#ifndef DUNE_MULTIDOMAINGRID_INTERSECTION_HH
#define DUNE_MULTIDOMAINGRID_INTERSECTION_HH

#include <dune/grid/common/intersection.hh>

#include <dune/grid/multidomaingrid/entitypointer.hh>
#include <dune/grid/multidomaingrid/hostgridaccessor.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/intersection.hh>

namespace Dune {

namespace mdgrid {

template<int codim, int dim, typename GridImp>
class EntityWrapper;


template<typename GridImp, typename HostIntersection_>
class IntersectionWrapper {

  template<class, class>
  friend class Dune::Intersection;

  template<typename,typename,typename>
  friend class subdomain::IntersectionWrapper;

  template<typename,typename>
  friend class IntersectionIteratorWrapper;

  using HostIntersection = HostIntersection_;
  using EntityPointer    = typename GridImp::Traits::template Codim<0>::EntityPointer;
  using Entity           = typename GridImp::Traits::template Codim<0>::Entity;
  using Geometry         = typename GridImp::Traits::template Codim<1>::Geometry;
  using LocalGeometry    = typename GridImp::Traits::template Codim<1>::LocalGeometry;
  using EntityWrapper    = Dune::mdgrid::EntityWrapper<0,GridImp::dimension,GridImp>;

  static const int dimension = GridImp::dimension;
  static const int dimensionworld = GridImp::dimensionworld;

  using ctype            = typename GridImp::ctype;
  using GlobalCoords     = FieldVector<ctype,dimensionworld>;
  using LocalCoords      = FieldVector<ctype,dimension - 1>;



protected:

  IntersectionWrapper() = default;

  explicit IntersectionWrapper(const HostIntersection& hostIntersection)
    : _hostIntersection(hostIntersection)
  {}

  explicit IntersectionWrapper(HostIntersection&& hostIntersection)
    : _hostIntersection(std::move(hostIntersection))
  {}

private:

  const HostIntersection& hostIntersection() const {
    return _hostIntersection;
  }

  bool equals(const IntersectionWrapper& rhs) const {
    return _hostIntersection == rhs._hostIntersection;
  }

  bool boundary() const {
    return _hostIntersection.boundary();
  }

  int boundaryId() const {
    return _hostIntersection.boundaryId();
  }

  std::size_t boundarySegmentIndex() const {
    return _hostIntersection.boundarySegmentIndex();
  }

  bool neighbor() const {
    return _hostIntersection.neighbor();
  }

  Entity inside() const {
    return {EntityWrapper(_hostIntersection.inside())};
  }

  Entity outside() const {
    return {EntityWrapper(_hostIntersection.outside())};
  }

  bool conforming() const {
    return _hostIntersection.conforming();
  }

  LocalGeometry geometryInInside() const {
    return LocalGeometry(_hostIntersection.geometryInInside());
  }

  LocalGeometry geometryInOutside() const {
    return LocalGeometry(_hostIntersection.geometryInOutside());
  }

  Geometry geometry() const {
    return Geometry(_hostIntersection.geometry());
  }

  GeometryType type() const {
    return _hostIntersection.type();
  }

  int indexInInside() const {
    return _hostIntersection.indexInInside();
  }

  int indexInOutside() const {
    return _hostIntersection.indexInOutside();
  }

  GlobalCoords outerNormal(const LocalCoords& local) const {
    return _hostIntersection.outerNormal(local);
  }

  GlobalCoords integrationOuterNormal(const LocalCoords& local) const {
    return _hostIntersection.integrationOuterNormal(local);
  }

  GlobalCoords unitOuterNormal(const LocalCoords& local) const {
    return _hostIntersection.unitOuterNormal(local);
  }

  GlobalCoords centerUnitOuterNormal() const {
    return _hostIntersection.centerUnitOuterNormal();
  }

private:

  HostIntersection _hostIntersection;

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_INTERSECTION_HH
