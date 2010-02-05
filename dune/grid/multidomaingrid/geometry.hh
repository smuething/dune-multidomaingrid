#ifndef DUNE_MULTIDOMAINGRID_GEOMETRY_HH
#define DUNE_MULTIDOMAINGRID_GEOMETRY_HH

#include <dune/grid/common/geometry.hh>

namespace Dune {

namespace mdgrid {

template<int mydim, int coorddim, typename GridImp>
class GeometryWrapper;


template<int mydim, int coorddim, typename GridImp>
class MakeableGeometryWrapper :
    public Dune::Geometry<mydim, coorddim, GridImp, GeometryWrapper>
{

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename, typename, typename>
  friend class IntersectionWrapper;

  template<typename, typename, typename, typename, typename>
  friend class SubDomainInterfaceIterator;

  typedef typename GridImp::HostGridType::Traits::template Codim<GridImp::dimension-mydim>::Geometry HostGeometry;

  // MakeableGeometryWrapper(const HostGeometry& geometry) :
  //   GridImp::template Codim<GridImp::dimension-mydim>::Geometry(GeometryWrapper<mydim,coorddim,GridImp>(geometry))
  // {}

  MakeableGeometryWrapper() :
    GridImp::template Codim<GridImp::dimension-mydim>::Geometry(GeometryWrapper<mydim,coorddim,GridImp>())
  {}

  void reset(const HostGeometry& geometry) const {
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
class GeometryWrapper
{

  template<int, int, typename>
  friend class MakeableGeometryWrapper;

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

  GeometryType type() const {
    return _wrappedGeometry->type();
  }

  int corners() const {
    return _wrappedGeometry->corners();
  }

  const GlobalCoords& operator[](int i) const {
    return _wrappedGeometry->operator[](i);
  }

  bool affine() const {
    return _wrappedGeometry->affine();
  }

  GlobalCoords corner(int i) const {
    return _wrappedGeometry->corner(i);
  }

  GlobalCoords global(const LocalCoords& local) const {
    return _wrappedGeometry->global(local);
  }

  LocalCoords local(const GlobalCoords& global) const {
    return _wrappedGeometry->local(global);
  }

  bool checkInside(const LocalCoords& local) const {
    return _wrappedGeometry->checkInside(local);
  }

  ctype integrationElement(const LocalCoords& local) const {
    return _wrappedGeometry->integrationElement(local);
  }

  ctype volume() const {
    return _wrappedGeometry->volume();
  }

  GlobalCoords center() const {
    return _wrappedGeometry->center();
  }

  const FieldMatrix<ctype,mydimension,coorddimension>&
  jacobianTransposed(const LocalCoords& local) const {
    return _wrappedGeometry->jacobianTransposed(local);
  }

  const FieldMatrix<ctype,coorddimension,mydimension>&
  jacobianInverseTransposed(const LocalCoords& local) const {
    return _wrappedGeometry->jacobianInverseTransposed(local);
  }

private:

  typedef const HostGeometry* ConstHostGeometryPointer;

  mutable ConstHostGeometryPointer _wrappedGeometry;

  void reset(const HostGeometry& geometry) const {
    _wrappedGeometry = &geometry;
  }

  bool isSet() const {
    return _wrappedGeometry != NULL;
  }

  void clear() const {
    _wrappedGeometry = NULL;
  }

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_GEOMETRY_HH
