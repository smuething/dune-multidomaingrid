#ifndef DUNE_MULTIDOMAINGRID_GEOMETRY_HH
#define DUNE_MULTIDOMAINGRID_GEOMETRY_HH

namespace Dune {

namespace mdgrid {

template<int mydim, int coorddim, typename GridImp>
class GeometryWrapper :
    public Dune::Geometry<mydim, coorddim, GridImp, GeometryWrapper>
{

public:
  
  typedef typename GridImp::ctype ctype;
  static const int dimension = GridImp::dimension;
  static const int dimensionworld = GridImp::dimensionworld;
  static const int mydimension = mydim;
  static const int coorddimension = coorddim;

private:

  typedef FieldVector<ctype,coorddimension> GlobalCoords;
  typedef FieldVector<ctype,mydimension> LocalCoords;
  typedef typename GridImp::HostGridType::Traits::template Codim<dimension-mydim>::Geometry HostGridGeometry;

public:

  GeometryType type() const {
    return _wrappedGeometry.type();
  }

  int corners() const {
    return _wrappedGeometry.corners();
  }

  const GlobalCoords& operator[](int i) const {
    return _wrappedGeometry[i];
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

  const FieldMatrix<ctype,mydimension,coorddimension>&
  jacobianTransposed(const LocalCoords& local) const {
    return _wrappedGeometry.jacobianTransposed(local);
  }

  const FieldMatrix<ctype,coorddimension,mydimension>&
  jacobianInverseTransposed(const LocalCoords& local) const {
    return _wrappedGeometry.jacobianInverseTransposed(local);
  }

private:

  const HostGridGeometry& _wrappedGeometry;

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_GEOMETRY_HH
