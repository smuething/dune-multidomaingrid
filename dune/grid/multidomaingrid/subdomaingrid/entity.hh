#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_ENTITY_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_ENTITY_HH

namespace Dune {

namespace mdgrid {

template<typename, typename>
class MultiDomainGrid;

namespace subdomain {

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

template<typename MDGrid>
class SubDomainGrid;

template<int codim, int dim, typename GridImp>
class MakeableEntityWrapper :
    public GridImp::template Codim<codim>::Entity
{

  template<int, typename>
  friend class EntityPointerWrapper;

  template<int, PartitionIteratorType, typename>
  friend class LeafIteratorWrapper;

  template<int, PartitionIteratorType, typename>
  friend class LevelIteratorWrapper;

  template<typename>
  friend class SubDomainGrid;

  template<typename>
  friend class HierarchicIteratorWrapper;

  typedef typename GridImp::MDGridType::Traits::template Codim<codim>::EntityPointer MultiDomainEntityPointer;


  MakeableEntityWrapper(const GridImp& grid, const MultiDomainEntityPointer& multiDomainEntityPointer) :
    GridImp::template Codim<codim>::Entity(EntityWrapper<codim,dim,const GridImp>(grid,multiDomainEntityPointer))
  {}

  void reset(const MultiDomainEntityPointer& multiDomainEntityPointer) {
    this->getRealImp().reset(multiDomainEntityPointer);
  }

  void compactify() {
    this->getRealImp().compactify();
  }

  const MultiDomainEntityPointer& multiDomainEntityPointer() const {
    return this->getRealImp().multiDomainEntityPointer();
  }

};


template<int codim, int dim, typename GridImp>
class EntityWrapper :
    public EntityDefaultImplementation<codim,dim,GridImp,EntityWrapper>
{

  template<int, int, typename>
  friend class MakeableEntityWrapper;

  template<int, typename>
  friend class EntityPointerWrapper;

  template<typename>
  friend class SubDomainGrid;

  template<typename,typename>
  friend class Dune::mdgrid::MultiDomainGrid;

  template<int, int, typename, template<int,int,typename> class>
  friend class Entity;

  typedef typename GridImp::MultiDomainGrid::Traits::template Codim<codim>::EntityPointer MultiDomainEntityPointer;
  typedef typename GridImp::HostGridType::Traits::template Codim<codim>::EntityPointer HostEntityPointer;
  typedef typename GridImp::HostGridType::Traits::template Codim<codim>::Entity HostEntity;

public:

  typedef typename GridImp::template Codim<codim>::Geometry Geometry;

  EntityWrapper(const GridImp& grid, const MultiDomainEntityPointer& e) :
    _grid(grid),
    _multiDomainEntityPointer(e)
  {}

  // copy constructor. The default constructor does not work correctly here,
  // as it will initialise the geometry with a pointer to the host geometry
  // of rhs, causing geometry access to fail if rhs gets destructed.
  EntityWrapper(const EntityWrapper& rhs) :
    _grid(rhs._grid),
    _multiDomainEntityPointer(rhs._multiDomainEntityPointer)
  {}

  int level() const {
    return _multiDomainEntityPointer->level();
  }

  PartitionType partitionType() const {
    return _multiDomainEntityPointer->partitionType();
  }

  template<int cc>
  int count() const {
    return _multiDomainEntityPointer->template count<cc>();
  }

  const Geometry& geometry() const {
    if (!_geometry.isSet()) {
      _geometry.reset(hostEntity().geometry());
    }
    return _geometry;
  }

private:

  const GridImp& _grid;
  MultiDomainEntityPointer _multiDomainEntityPointer;
  MakeableGeometryWrapper<Geometry::mydimension,Geometry::coorddimension,GridImp> _geometry;

  const EntityWrapper& operator=(const EntityWrapper& rhs) {
    assert(_grid == rhs._grid);
    reset(rhs._multiDomainEntityPointer);
    return *this;
  }

  void reset(const MultiDomainEntityPointer& multiDomainEntityPointer) {
    if (_multiDomainEntityPointer != multiDomainEntityPointer) {
      _geometry.clear();
      _multiDomainEntityPointer = multiDomainEntityPointer;
    }
  }

  void compactify() {
    _geometry.clear();
    _multiDomainEntityPointer.compactify();
  }

  const MultiDomainEntityPointer& multiDomainEntityPointer() const {
    return _multiDomainEntityPointer;
  }

  const HostEntity& hostEntity() const {
    return _grid._grid.hostEntity(*_multiDomainEntityPointer);
  }

};


template<int dim, typename GridImp>
class EntityWrapper<0,dim,GridImp> :
    public EntityDefaultImplementation<0,dim,GridImp,EntityWrapper>
{

  template<int, int, typename>
  friend class MakeableEntityWrapper;

  template<int, typename>
  friend class EntityPointerWrapper;

  template<typename>
  friend class SubDomainGrid;

  template<typename,typename>
  friend class Dune::mdgrid::MultiDomainGrid;

  template<int, int, typename, template<int,int,typename> class>
  friend class Entity;

  typedef typename GridImp::MDGridType::Traits::template Codim<0>::EntityPointer MultiDomainEntityPointer;
  typedef typename GridImp::HostGridType::Traits::template Codim<0>::EntityPointer HostEntityPointer;
  typedef typename GridImp::HostGridType::Traits::template Codim<0>::Entity HostEntity;

public:

  typedef typename GridImp::template Codim<0>::Geometry Geometry;
  typedef typename GridImp::template Codim<0>::LocalGeometry LocalGeometry;
  typedef typename GridImp::Traits::LeafIntersectionIterator LeafIntersectionIterator;
  typedef typename GridImp::Traits::LevelIntersectionIterator LevelIntersectionIterator;
  typedef typename GridImp::Traits::HierarchicIterator HierarchicIterator;
  typedef typename GridImp::Traits::template Codim<0>::EntityPointer EntityPointer;


  EntityWrapper(const GridImp& grid, const MultiDomainEntityPointer& e) :
    _grid(grid),
    _multiDomainEntityPointer(e)
  {}

  // copy constructor. The default constructor does not work correctly here,
  // as it will initialise the geometry with a pointer to the host geometry
  // of rhs, causing geometry access to fail if rhs gets destructed.
  EntityWrapper(const EntityWrapper& rhs) :
    _grid(rhs._grid),
    _multiDomainEntityPointer(rhs._multiDomainEntityPointer)
  {}

  int level() const {
    return _multiDomainEntityPointer->level();
  }

  PartitionType partitionType() const {
    return _multiDomainEntityPointer->partitionType();
  }

  const Geometry& geometry() const {
    if (!_geometry.isSet()) {
      _geometry.reset(hostEntity().geometry());
    }
    return _geometry;
  }

  template<int cc>
  int count() const {
    return hostEntity().template count<cc>();
  }

  template<int cc>
  typename GridImp::template Codim<cc>::EntityPointer subEntity(int i) const {
    return EntityPointerWrapper<cc,GridImp>(_grid,_multiDomainEntityPointer->subEntity<cc>(i));
  }

  LeafIntersectionIterator ileafbegin() const {
    return LeafIntersectionIteratorWrapper<GridImp>(_grid,_multiDomainEntityPointer->ileafbegin());
  }

  LeafIntersectionIterator ileafend() const {
    return LeafIntersectionIteratorWrapper<GridImp>(_grid,_multiDomainEntityPointer->ileafend());
  }

  LevelIntersectionIterator ilevelbegin() const {
    return LevelIntersectionIteratorWrapper<GridImp>(_grid,level(),_multiDomainEntityPointer->ilevelbegin());
  }

  LevelIntersectionIterator ilevelend() const {
    return LevelIntersectionIteratorWrapper<GridImp>(_grid,level(),_multiDomainEntityPointer->ilevelend());
  }

  EntityPointer father() const {
    return EntityPointerWrapper<0,GridImp>(_grid,_multiDomainEntityPointer->father());
  }

  bool hasFather() const {
    return _multiDomainEntityPointer->hasFather();
  }

  bool isLeaf() const {
    return _multiDomainEntityPointer->isLeaf();
  }

  bool isRegular() const {
    return _multiDomainEntityPointer->isRegular();
  }

  const LocalGeometry& geometryInFather() const {
    if (!_fatherGeometry.isSet()) {
      _fatherGeometry.reset(hostEntity().geometryInFather());
    }
    return _fatherGeometry;
  }

  HierarchicIterator hbegin(int maxLevel) const {
    return HierarchicIteratorWrapper<GridImp>(_grid,
                                              _multiDomainEntityPointer->hbegin(maxLevel),
                                              _multiDomainEntityPointer->hend(maxLevel));
  }

  HierarchicIterator hend(int maxLevel) const {
    return HierarchicIteratorWrapper<GridImp>(_grid,
                                              _multiDomainEntityPointer->hend(maxLevel),
                                              _multiDomainEntityPointer->hend(maxLevel));
  }

  bool isNew() const {
    return _multiDomainEntityPointer->isNew();
  }

  bool mightVanish() const {
    return _multiDomainEntityPointer->mightVanish();
  }

private:

  const GridImp& _grid;
  MultiDomainEntityPointer _multiDomainEntityPointer;
  MakeableGeometryWrapper<Geometry::mydimension,Geometry::coorddimension,GridImp> _geometry;
  MakeableLocalGeometryWrapper<LocalGeometry::mydimension,LocalGeometry::coorddimension,GridImp> _fatherGeometry;

  const EntityWrapper& operator=(const EntityWrapper& rhs) {
    assert(_grid == rhs._grid);
    reset(rhs._multiDomainEntityPointer);
    return *this;
  }

  void reset(const MultiDomainEntityPointer& multiDomainEntityPointer) {
    if (_multiDomainEntityPointer != multiDomainEntityPointer) {
      _geometry.clear();
      _fatherGeometry.clear();
      _multiDomainEntityPointer = multiDomainEntityPointer;
    }
  }

  void compactify() {
    _geometry.clear();
    _fatherGeometry.clear();
    _multiDomainEntityPointer.compactify();
  }

  const MultiDomainEntityPointer& multiDomainEntityPointer() const {
    return _multiDomainEntityPointer;
  }

  const HostEntity& hostEntity() const {
    return _grid._grid.hostEntity(*_multiDomainEntityPointer);
  }

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_ENTITY_HH
