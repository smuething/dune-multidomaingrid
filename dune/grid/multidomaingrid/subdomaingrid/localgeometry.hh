#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_LOCALGEOMETRY_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_LOCALGEOMETRY_HH

#include <dune/grid/common/geometry.hh>

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<int mydim, int coorddim, typename GridImp>
class LocalGeometryWrapper
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
  typedef typename GridImp::HostGridType::Traits::template Codim<dimension-mydim>::LocalGeometry HostLocalGeometry; //TODO: fix this

public:

  typedef typename HostLocalGeometry::JacobianInverseTransposed JacobianInverseTransposed;
  typedef typename HostLocalGeometry::JacobianTransposed JacobianTransposed;

  GeometryType type() const {
    return _hostLocalGeometry.type();
  }

  int corners() const {
    return _hostLocalGeometry.corners();
  }

  bool affine() const {
    return _hostLocalGeometry.affine();
  }

  GlobalCoords corner(int i) const {
    return _hostLocalGeometry.corner(i);
  }

  GlobalCoords global(const LocalCoords& local) const {
    return _hostLocalGeometry.global(local);
  }

  LocalCoords local(const GlobalCoords& global) const {
    return _hostLocalGeometry.local(global);
  }

  bool checkInside(const LocalCoords& local) const {
    return _hostLocalGeometry.checkInside(local);
  }

  ctype integrationElement(const LocalCoords& local) const {
    return _hostLocalGeometry.integrationElement(local);
  }

  ctype volume() const {
    return _hostLocalGeometry.volume();
  }

  GlobalCoords center() const {
    return _hostLocalGeometry.center();
  }

  const JacobianTransposed
  jacobianTransposed(const LocalCoords& local) const {
    return _hostLocalGeometry.jacobianTransposed(local);
  }

  const JacobianInverseTransposed
  jacobianInverseTransposed(const LocalCoords& local) const {
    return _hostLocalGeometry.jacobianInverseTransposed(local);
  }

private:

  const HostLocalGeometry _hostLocalGeometry;

  LocalGeometryWrapper(const HostLocalGeometry& hostLocalGeometry)
    : _hostLocalGeometry(hostLocalGeometry)
  {}

};

} // namespace subdomain

} // namespace mdgrid


namespace FacadeOptions {

template< int mydim, int coorddim, class GridImp >
struct StoreGeometryReference< mydim, coorddim, GridImp, mdgrid::subdomain::LocalGeometryWrapper >
{
        static const bool v = false;
};

} // namespace FacadeOptions

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_LOCALGEOMETRY_HH
