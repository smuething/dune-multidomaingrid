#ifndef DUNE_MULTIDOMAINGRID_ENTITY_HH
#define DUNE_MULTIDOMAINGRID_ENTITY_HH

namespace Dune {

namespace mdgrid {

template<int codim, int dim, typename GridImp>
class EntityWrapper;

template<int codim, typename GridImp>
class EntityPointerWrapper;

template<typename GridImp>
class LeafIntersectionIteratorWrapper;

template<typename GridImp>
class LevelIntersectionIteratorWrapper;

template<typename GridImp>
class HierarchicIteratorWrapper;

template<int codim, int dim, typename GridImp>
class MakeableEntityWrapper :
    public GridImp::template Codim<codim>::Entity
{
  
  typedef typename GridImp::HostGridType::Traits::template Codim<codim>::EntityPointer HostEntityPointer;
  
  MakeableEntityWrapper(const GridImp& grid, const HostEntityPointer& hostEntityPointer) :
    GridImp::template Codim<codim>::Entity(EntityWrapper<codim,dim,const GridImp>(grid,hostEntityPointer))
  {}

  void reset(const HostEntityPointer& hostEntityPointer) {
    this->getRealImp().reset(hostEntityPointer);
  }

  void compactify() {
    this->getRealImp().compactify();
  }

};


template<int codim, int dim, typename GridImp>
class EntityWrapper :
    public EntityDefaultImplementation<codim,dim,GridImp,EntityWrapper>
{

  typedef typename GridImp::HostGridType::Traits::template Codim<codim>::EntityPointer HostEntityPointer;

public:

  typedef typename GridImp::template Codim<codim>::Geometry Geometry;
  
  EntityWrapper(const GridImp& grid, const HostEntityPointer& e) :
    _hostEntityPointer(&e)
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

  const Geometry& geometry() const {
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
  typedef typename GridImp::Traits::LeafIntersectionIterator LeafIntersectionIterator;
  typedef typename GridImp::Traits::LevelIntersectionIterator LevelIntersectionIterator;
  typedef typename GridImp::Traits::HierarchicIterator HierarchicIterator;
  typedef typename GridImp::Traits::template Codim<0>::EntityPointer EntityPointer;



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
    return LeafIntersectionIteratorWrapper<GridImp>(_hostEntityPointer->ileafbegin());
  }

  LeafIntersectionIterator ileafend() const {
    return LeafIntersectionIteratorWrapper<GridImp>(_hostEntityPointer->ileafend());
  }

  LevelIntersectionIterator ilevelbegin() const {
    return LevelIntersectionIteratorWrapper<GridImp>(_hostEntityPointer->ilevelbegin());
  }

  LevelIntersectionIterator ilevelend() const {
    return LevelIntersectionIteratorWrapper<GridImp>(_hostEntityPointer->ilevelend());
  }

  EntityPointer father() const {
    return EntityPointerWrapper<0,GridImp>(_hostEntityPointer->father());
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

  HierarchicIterator hbegin(int maxLevel) const {
    return HierarchicIteratorWrapper<GridImp>(_hostEntityPointer->hbegin(maxLevel));
  }

  HierarchicIterator hend(int maxLevel) const {
    return HierarchicIteratorWrapper<GridImp>(_hostEntityPointer->hend(maxLevel));
  }

  bool isNew() const {
    return _hostEntityPointer->isNew();
  }

  bool mightVanish() const {
    return _hostEntityPointer->mightVanish();
  }

private:

  typedef typename GridImp::HostGridType::Traits::template Codim<0>::EntityPointer HostEntityPointer;

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
