#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_LOCALGEOMETRY_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_LOCALGEOMETRY_HH

#include <dune/grid/common/geometry.hh>

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<int mydim, int coorddim, typename GridImp>
class LocalGeometryWrapper;


template<int mydim, int coorddim, typename GridImp>
class MakeableLocalGeometryWrapper :
    public Dune::Geometry<mydim, coorddim, GridImp, LocalGeometryWrapper>
{

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename, typename, typename, typename>
  friend class IntersectionWrapper;

  typedef typename GridImp::HostGridType::Traits::template Codim<GridImp::dimension-mydim>::LocalGeometry HostLocalGeometry;

  MakeableLocalGeometryWrapper() :
    GridImp::template Codim<GridImp::dimension-mydim>::LocalGeometry(LocalGeometryWrapper<mydim,coorddim,GridImp>())
  {}

  void reset(const HostLocalGeometry& geometry) const {
    this->getRealImp().reset(geometry);
  }

  void clear() const {
    this->getRealImp().clear();
  }

  bool isSet() const {
    return this->getRealImp().isSet();
  }

};


template<int mydim, int coorddim, typename GridImp>
class LocalGeometryWrapper
{

  template<int, int, typename>
  friend class MakeableLocalGeometryWrapper;

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

  GeometryType type() const {
    return _hostLocalGeometry->type();
  }

  int corners() const {
    return _hostLocalGeometry->corners();
  }

  const GlobalCoords& operator[](int i) const {
    return _hostLocalGeometry->operator[](i);
  }

  bool affine() const {
    return _hostLocalGeometry->affine();
  }

  GlobalCoords corner(int i) const {
    return _hostLocalGeometry->corner(i);
  }

  GlobalCoords global(const LocalCoords& local) const {
    return _hostLocalGeometry->global(local);
  }

  LocalCoords local(const GlobalCoords& global) const {
    return _hostLocalGeometry->local(global);
  }

  bool checkInside(const LocalCoords& local) const {
    return _hostLocalGeometry->checkInside(local);
  }

  ctype integrationElement(const LocalCoords& local) const {
    return _hostLocalGeometry->integrationElement(local);
  }

  ctype volume() const {
    return _hostLocalGeometry->volume();
  }

  GlobalCoords center() const {
    return _hostLocalGeometry->center();
  }

  const FieldMatrix<ctype,mydimension,coorddimension>&
  jacobianTransposed(const LocalCoords& local) const {
    return _hostLocalGeometry->jacobianTransposed(local);
  }

  const FieldMatrix<ctype,coorddimension,mydimension>&
  jacobianInverseTransposed(const LocalCoords& local) const {
    return _hostLocalGeometry->jacobianInverseTransposed(local);
  }

private:

  typedef const HostLocalGeometry* ConstHostLocalGeometryPointer;

  mutable ConstHostLocalGeometryPointer _hostLocalGeometry;

  void reset(const HostLocalGeometry& localGeometry) const {
    _hostLocalGeometry = &localGeometry;
  }

  bool isSet() const {
    return _hostLocalGeometry != NULL;
  }

  void clear() const {
    _hostLocalGeometry = NULL;
  }

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_LOCALGEOMETRY_HH
