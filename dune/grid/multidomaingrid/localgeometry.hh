#ifndef DUNE_MULTIDOMAINGRID_LOCALGEOMETRY_HH
#define DUNE_MULTIDOMAINGRID_LOCALGEOMETRY_HH

#include <dune/grid/common/geometry.hh>

namespace Dune {

namespace mdgrid {


template<int mydim, int coorddim, typename GridImp>
class LocalGeometryWrapper
{

  template<int,int,typename>
  friend class EntityWrapper;

  template<typename,typename>
  friend class IntersectionWrapper;

  template<typename,typename,typename,typename>
  friend class SubDomainInterface;

public:

  typedef typename GridImp::ctype ctype;
  static const int dimension = GridImp::dimension;
  static const int dimensionworld = GridImp::dimensionworld;
  static const int mydimension = mydim;
  static const int coorddimension = coorddim;

private:

  typedef FieldVector<ctype,coorddimension> GlobalCoords;
  typedef FieldVector<ctype,mydimension> LocalCoords;
  typedef typename GridImp::HostGrid::Traits::template Codim<dimension-mydim>::LocalGeometry HostLocalGeometry; //TODO: fix this

public:

  typedef typename HostLocalGeometry::JacobianInverseTransposed JacobianInverseTransposed;
  typedef typename HostLocalGeometry::JacobianTransposed JacobianTransposed;

  GeometryType type() const {
    return _wrappedLocalGeometry.type();
  }

  int corners() const {
    return _wrappedLocalGeometry.corners();
  }

  bool affine() const {
    return _wrappedLocalGeometry.affine();
  }

  GlobalCoords corner(int i) const {
    return _wrappedLocalGeometry.corner(i);
  }

  GlobalCoords global(const LocalCoords& local) const {
    return _wrappedLocalGeometry.global(local);
  }

  LocalCoords local(const GlobalCoords& global) const {
    return _wrappedLocalGeometry.local(global);
  }

  bool checkInside(const LocalCoords& local) const {
    return _wrappedLocalGeometry.checkInside(local);
  }

  ctype integrationElement(const LocalCoords& local) const {
    return _wrappedLocalGeometry.integrationElement(local);
  }

  ctype volume() const {
    return _wrappedLocalGeometry.volume();
  }

  GlobalCoords center() const {
    return _wrappedLocalGeometry.center();
  }

  const JacobianTransposed
  jacobianTransposed(const LocalCoords& local) const {
    return _wrappedLocalGeometry.jacobianTransposed(local);
  }

  const JacobianInverseTransposed
  jacobianInverseTransposed(const LocalCoords& local) const {
    return _wrappedLocalGeometry.jacobianInverseTransposed(local);
  }

private:

  const HostLocalGeometry _wrappedLocalGeometry;

  LocalGeometryWrapper(const HostLocalGeometry& wrappedLocalGeometry)
    : _wrappedLocalGeometry(wrappedLocalGeometry)
  {}


};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_LOCALGEOMETRY_HH
