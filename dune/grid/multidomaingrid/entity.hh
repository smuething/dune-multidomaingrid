#ifndef DUNE_MULTIDOMAINGRID_ENTITY_HH
#define DUNE_MULTIDOMAINGRID_ENTITY_HH

namespace Dune {

namespace mdgrid {

template<int codim, int dim, typename GridImp>
class MakeableEntity :
    public GridImp::template Codim<codim>::Entity
{
  
  typedef typename GridImp::HostGridType::Traits::template Codim<codim>::EntityPointer HostEntityPointer;
  
  MakeableEntity(const GridImp& grid, const HostEntityPointer& hostEntityPointer) :
    GridImp::template Codim<codim>::Entity(EntityWrapper<codim,const GridImp>(grid,hostEntityPointer))
  {}

  void reset(const HostEntityPointer& hostEntityPointer) {
    realEntity().reset(hostEntityPointer);
  }

  void compactify() {
    realEntity().compactify();
  }

};


template<int codim, int dim, typename GridImp>
class EntityWrapper :
    public EntityDefaultImplementation<codim,dim,GridImp,EntityWrapper>
{

public:

  typedef typename GridImp::template Codim<codim>::Geometry Geometry;
  
  EntityWrapper(const GridImp& grid, const HostEntityPointer& e) :
    _hostEntityPointer(&hostEntityPointer)
  {}

  int level() const {
    return _hostEntityPointer->level();
  }

  PartitionType partitionType() const {
    return _hostEntityPointer->partitionType();
  }

  template<int cc>
  int count() const {
    return _hostEntityPointer->template count<cc>();
  }

  const Geometry& geometry const() {
    if (_geometry == NULL) {
      _geometry.reset(new MakeableInterfaceObject<Geometry>(_hostEntityPointer->geometry()));
    }
    return *_geometry;
  }

  private:

  HostEntityPointer _hostEntityPointer;
  mutable boost::scoped_ptr<MakeableInterfaceObject<Geometry> > _geometry;

  void reset(const HostEntityPointer& hostEntityPointer) {
    if (_hostEntityPointer != hostEntityPointer) {
      _geometry.reset(NULL);
      _hostEntityPointer = hostEntityPointer;
    }
  }

  void compactify() {
    _geometry.reset(NULL);
    _hostEntityPointer.compactify();
  }

};


template<int dim, typename GridImp>
class EntityWrapper<0,dim,GridImp> :
    public EntityDefaultImplementation<0,dim,GridImp,EntityWrapper>
{

public:

  typedef typename GridImp::template Codim<0>::Geometry Geometry;
  typedef typename GridImp::template Codim<0>::LocalGeometry LocalGeometry;

  int level() const {
    return _hostEntityPointer->level();
  }

  PartitionType partitionType() const {
    return _hostEntityPointer->partitionType();
  }

  const Geometry& geometry() const {
    if (_geometry == NULL) {
      _geometry.reset(new MakeableInterfaceObject<Geometry>(_hostEntityPointer->geometry()));
    }
    return *_geometry;
  }

  template<int cc>
  int count() const {
    return _hostEntityPointer->template count<cc>();
  }

  template<int cc>
  typename GridImp::template Codim<cc>::EntityPointer entity(int i) const {
  }

  template<int cc>
  typename GridImp::template Codim<cc>::EntityPointer subEntity(int i) const {
  }

  LeafIntersectionIterator ileafbegin() const {
  }

  LeafIntersectionIterator ileafend() const {
  }

  LevelIntersectionIterator ilevelbegin() const {
  }

  LevelIntersectionIterator ilevelend() const {
  }

  EntityPointer father() const {
  }

  bool isLeaf() const {
    return _hostEntityPointer->isLeaf();
  }

  bool isRegular() const {
    return _hostEntityPointer->isRegular();
  }

  const LocalGeometry& geometryInFather() const {
    if (_fatherGeometry == NULL) {
      _fatherGeometry.reset(new MakeableInterfaceObject<LocalGeometry>(_hostEntityPointer->geometryInFather()));
    }
    return *_fatherGeometry;
  }

  HierarchicIterator hbegin(int maxlevel) const {
  }

  HierarchicIterator hend(int maxlevel) const {
  }

  bool isNew() const {
    return _hostEntityPointer->isNew();
  }

  bool mightVanish() const {
    return _hostEntityPointer->mightVanish();
  }

private:
  HostEntityPointer _hostEntityPointer;
  mutable boost::scoped_ptr<MakeableInterfaceObject<Geometry> > _geometry;
  mutable boost::scoped_ptr<MakeableInterfaceObject<LocalGeometry> > _fatherGeometry;

  void reset(const HostEntityPointer& hostEntityPointer) {
    if (_hostEntityPointer != hostEntityPointer) {
      _geometry.reset(NULL);
      _fatherGeometry.reset(NULL);
      _hostEntityPointer = hostEntityPointer;
    }
  }

  void compactify() {
    _geometry.reset(NULL);
    _fatherGeometry.reset(NULL);
    _hostEntityPointer.compactify();
  }

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_ENTITY_HH
