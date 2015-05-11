#ifndef DUNE_MULTIDOMAINGRID_GEOMETRY_HH
#define DUNE_MULTIDOMAINGRID_GEOMETRY_HH

#include <dune/grid/common/geometry.hh>

namespace Dune {

namespace mdgrid {

template<int mydim, int coorddim, typename GridImp>
class GeometryWrapper
{

  template<int,int,typename>
  friend class EntityWrapperBase;

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
  typedef typename GridImp::HostGrid::Traits::template Codim<dimension-mydim>::Geometry HostGeometry; //TODO: fix this

public:

  typedef typename HostGeometry::JacobianInverseTransposed JacobianInverseTransposed;
  typedef typename HostGeometry::JacobianTransposed JacobianTransposed;

  GeometryType type() const {
    return _wrappedGeometry.type();
  }

  int corners() const {
    return _wrappedGeometry.corners();
  }

  bool affine() const {
    return _wrappedGeometry.affine();
  }

  GlobalCoords corner(int i) const {
    return _wrappedGeometry.corner(i);
  }

  GlobalCoords global(const LocalCoords& local) const {
    return _wrappedGeometry.global(local);
  }

  LocalCoords local(const GlobalCoords& global) const {
    return _wrappedGeometry.local(global);
  }

  bool checkInside(const LocalCoords& local) const {
    return _wrappedGeometry.checkInside(local);
  }

  ctype integrationElement(const LocalCoords& local) const {
    return _wrappedGeometry.integrationElement(local);
  }

  ctype volume() const {
    return _wrappedGeometry.volume();
  }

  GlobalCoords center() const {
    return _wrappedGeometry.center();
  }

  const JacobianTransposed
  jacobianTransposed(const LocalCoords& local) const {
    return _wrappedGeometry.jacobianTransposed(local);
  }

  const JacobianInverseTransposed
  jacobianInverseTransposed(const LocalCoords& local) const {
    return _wrappedGeometry.jacobianInverseTransposed(local);
  }

private:

  GeometryWrapper(const HostGeometry& wrappedGeometry)
    : _wrappedGeometry(wrappedGeometry)
  {}

  const HostGeometry _wrappedGeometry;

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_GEOMETRY_HH
