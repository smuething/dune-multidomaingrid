#ifndef DUNE_MULTIDOMAINGRID_LOCALGEOMETRY_HH
#define DUNE_MULTIDOMAINGRID_LOCALGEOMETRY_HH

#include <dune/grid/common/geometry.hh>

namespace Dune {

namespace mdgrid {

template<int mydim, int coorddim, typename GridImp>
class LocalGeometryWrapper;


template<int mydim, int coorddim, typename GridImp>
class MakeableLocalGeometryWrapper :
    public Dune::Geometry<mydim, coorddim, GridImp, LocalGeometryWrapper>
{

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename, typename, typename>
  friend class IntersectionWrapper;

  template<typename, typename, typename, typename>
  friend class SubDomainInterface;

  typedef typename GridImp::HostGridType::Traits::template Codim<GridImp::dimension-mydim>::LocalGeometry HostLocalGeometry;

  MakeableLocalGeometryWrapper() :
    GridImp::template Codim<GridImp::dimension-mydim>::LocalGeometry(LocalGeometryWrapper<mydim,coorddim,GridImp>())
  {}

  void reset(const HostLocalGeometry& geometry) const {
    GridImp::getRealImplementation(*this).reset(geometry);
  }

  void clear() const {
    GridImp::getRealImplementation(*this).clear();
  }

  bool isSet() const {
    return GridImp::getRealImplementation(*this).isSet();
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
    return _wrappedLocalGeometry->type();
  }

  int corners() const {
    return _wrappedLocalGeometry->corners();
  }

  const GlobalCoords& operator[](int i) const {
    return _wrappedLocalGeometry->operator[](i);
  }

  bool affine() const {
    return _wrappedLocalGeometry->affine();
  }

  GlobalCoords corner(int i) const {
    return _wrappedLocalGeometry->corner(i);
  }

  GlobalCoords global(const LocalCoords& local) const {
    return _wrappedLocalGeometry->global(local);
  }

  LocalCoords local(const GlobalCoords& global) const {
    return _wrappedLocalGeometry->local(global);
  }

  bool checkInside(const LocalCoords& local) const {
    return _wrappedLocalGeometry->checkInside(local);
  }

  ctype integrationElement(const LocalCoords& local) const {
    return _wrappedLocalGeometry->integrationElement(local);
  }

  ctype volume() const {
    return _wrappedLocalGeometry->volume();
  }

  GlobalCoords center() const {
    return _wrappedLocalGeometry->center();
  }

  const FieldMatrix<ctype,mydimension,coorddimension>&
  jacobianTransposed(const LocalCoords& local) const {
    return _wrappedLocalGeometry->jacobianTransposed(local);
  }

  const FieldMatrix<ctype,coorddimension,mydimension>&
  jacobianInverseTransposed(const LocalCoords& local) const {
    return _wrappedLocalGeometry->jacobianInverseTransposed(local);
  }

private:

  typedef const HostLocalGeometry* ConstHostLocalGeometryPointer;

  mutable ConstHostLocalGeometryPointer _wrappedLocalGeometry;

  void reset(const HostLocalGeometry& localGeometry) const {
    _wrappedLocalGeometry = &localGeometry;
  }

  bool isSet() const {
    return _wrappedLocalGeometry != NULL;
  }

  void clear() const {
    _wrappedLocalGeometry = NULL;
  }

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_LOCALGEOMETRY_HH
