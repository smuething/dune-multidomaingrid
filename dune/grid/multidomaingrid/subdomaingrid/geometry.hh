#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_GEOMETRY_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_GEOMETRY_HH

#include <dune/grid/common/geometry.hh>
#include <dune/common/version.hh>

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<int mydim, int coorddim, typename GridImp>
class GeometryWrapper
{

  template<int,int,typename>
  friend class EntityWrapper;

  template<typename,typename,typename,typename>
  friend class IntersectionWrapper;

public:

  typedef typename GridImp::ctype ctype;
  static const int dimension = GridImp::dimension;
  static const int dimensionworld = GridImp::dimensionworld;
  static const int mydimension = mydim;
  static const int coorddimension = coorddim;

private:

  typedef FieldVector<ctype,coorddimension> GlobalCoords;
  typedef FieldVector<ctype,mydimension> LocalCoords;
  typedef typename GridImp::HostGridType::Traits::template Codim<dimension-mydim>::Geometry HostGeometry; //TODO: fix this

public:

  typedef typename HostGeometry::JacobianInverseTransposed JacobianInverseTransposed;
  typedef typename HostGeometry::JacobianTransposed JacobianTransposed;

  GeometryType type() const {
    return _hostGeometry.type();
  }

  int corners() const {
    return _hostGeometry.corners();
  }

  bool affine() const {
    return _hostGeometry.affine();
  }

  GlobalCoords corner(int i) const {
    return _hostGeometry.corner(i);
  }

  GlobalCoords global(const LocalCoords& local) const {
    return _hostGeometry.global(local);
  }

  LocalCoords local(const GlobalCoords& global) const {
    return _hostGeometry.local(global);
  }

  bool checkInside(const LocalCoords& local) const {
    return _hostGeometry.checkInside(local);
  }

  ctype integrationElement(const LocalCoords& local) const {
    return _hostGeometry.integrationElement(local);
  }

  ctype volume() const {
    return _hostGeometry.volume();
  }

  GlobalCoords center() const {
    return _hostGeometry.center();
  }

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
  const JacobianTransposed
#else
  const JacobianTransposed&
#endif
  jacobianTransposed(const LocalCoords& local) const {
    return _hostGeometry.jacobianTransposed(local);
  }

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
  const JacobianInverseTransposed
#else
  const JacobianInverseTransposed&
#endif
  jacobianInverseTransposed(const LocalCoords& local) const {
    return _hostGeometry.jacobianInverseTransposed(local);
  }

private:

  const HostGeometry _hostGeometry;

  GeometryWrapper(const HostGeometry& hostGeometry)
    : _hostGeometry(hostGeometry)
  {}

};

} // namespace subdomain

} // namespace mdgrid


#if !DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
namespace FacadeOptions {

template< int mydim, int coorddim, class GridImp >
struct StoreGeometryReference< mydim, coorddim, GridImp, mdgrid::subdomain::GeometryWrapper >
{
        static const bool v = false;
};

} // namespace FacadeOptions
#endif

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_GEOMETRY_HH
