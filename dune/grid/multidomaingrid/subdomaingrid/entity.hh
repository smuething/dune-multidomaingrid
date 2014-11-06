#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_ENTITY_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_ENTITY_HH

#include <dune/grid/common/entity.hh>
#include <dune/grid/common/gridenums.hh>

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

template<typename, PartitionIteratorType>
class LevelGridView;

template<typename, PartitionIteratorType>
class LeafGridView;

template<typename HostES>
class EntitySeedWrapper
{

  typedef HostES HostEntitySeed;

  template<typename>
  friend class SubDomainGrid;

  template<typename,typename>
  friend class ::Dune::mdgrid::MultiDomainGrid;

  template<int, int, typename>
  friend class EntityWrapper;

public:

  static const std::size_t codimension = HostEntitySeed::codimension;

  EntitySeedWrapper()
  {}

  EntitySeedWrapper(const HostEntitySeed& hostEntitySeed)
    : _hostEntitySeed(hostEntitySeed)
  {}

  const HostEntitySeed& hostEntitySeed() const
  {
    return _hostEntitySeed;
  }

  bool isValid() const
  {
    return _hostEntitySeed.isValid();
  }

  HostEntitySeed _hostEntitySeed;

};

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

  template<typename, PartitionIteratorType>
  friend class LevelGridView;

  template<typename, PartitionIteratorType>
  friend class LeafGridView;


  typedef typename GridImp::MDGridType::Traits::template Codim<codim>::EntityPointer MultiDomainEntityPointer;

  template<typename MultiDomainIteratorOrEntityPointer>
  MakeableEntityWrapper(const GridImp& grid, const MultiDomainIteratorOrEntityPointer& multiDomainEntityPointer) :
    GridImp::template Codim<codim>::Entity(EntityWrapper<codim,dim,const GridImp>(grid,multiDomainEntityPointer))
  {}

  MakeableEntityWrapper& operator=(const MakeableEntityWrapper& rhs)
  {
    reset(rhs.multiDomainEntityPointer());
    return *this;
  }

  template<typename MultiDomainIteratorOrEntityPointer>
  void reset(const MultiDomainIteratorOrEntityPointer& multiDomainEntityPointer) {
    GridImp::getRealImplementation(*this).reset(multiDomainEntityPointer);
  }

  void compactify() {
    GridImp::getRealImplementation(*this).compactify();
  }

  const MultiDomainEntityPointer& multiDomainEntityPointer() const {
    return GridImp::getRealImplementation(*this).multiDomainEntityPointer();
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

  typedef EntitySeedWrapper<typename HostEntity::EntitySeed> EntitySeed;

  typedef typename GridImp::template Codim<codim>::Geometry Geometry;

  template<typename MultiDomainIteratorOrEntityPointer>
  EntityWrapper(const GridImp& grid, const MultiDomainIteratorOrEntityPointer& e) :
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
  int DUNE_DEPRECATED_MSG("Use subEntities instead") count() const {
    return _multiDomainEntityPointer->template count<cc>();
  }

  unsigned int subEntities(unsigned int codimSubEntitiy) const {
    return _multiDomainEntityPointer->subEntities(codimSubEntitiy);
  }

  Geometry geometry() const {
    return Geometry(hostEntity().geometry());
  }

  EntitySeed seed() const {
    return EntitySeed(hostEntity().seed());
  }

private:

  const GridImp& _grid;
  MultiDomainEntityPointer _multiDomainEntityPointer;

  const EntityWrapper& operator=(const EntityWrapper& rhs) {
    assert(_grid == rhs._grid);
    reset(rhs._multiDomainEntityPointer);
    return *this;
  }

  template<typename MultiDomainIteratorOrEntityPointer>
  void reset(const MultiDomainIteratorOrEntityPointer& multiDomainEntityPointer) {
    _multiDomainEntityPointer = multiDomainEntityPointer;
  }

  void compactify() {
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

  template<typename, PartitionIteratorType>
  friend class LevelGridView;

  template<typename, PartitionIteratorType>
  friend class LeafGridView;


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

  typedef EntitySeedWrapper<typename HostEntity::EntitySeed> EntitySeed;

  template<typename MultiDomainIteratorOrEntityPointer>
  EntityWrapper(const GridImp& grid, const MultiDomainIteratorOrEntityPointer& e) :
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

  Geometry geometry() const {
    return Geometry(hostEntity().geometry());
  }

  template<int cc>
  int DUNE_DEPRECATED_MSG("Use subEntities instead") count() const {
    return hostEntity().template count<cc>();
  }

  unsigned int subEntities(unsigned int codim) const {
    return hostEntity().subEntities(codim);
  }

  template<int cc>
  typename GridImp::template Codim<cc>::EntityPointer subEntity(int i) const {
    return EntityPointerWrapper<cc,GridImp>(_grid,_multiDomainEntityPointer->template subEntity<cc>(i));
  }

  bool equals(const EntityWrapper& other) const
  {
    return _multiDomainEntityPointer == other.multiDomainEntityPointer();
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

  LocalGeometry geometryInFather() const {
    return LocalGeometry(hostEntity().geometryInFather());
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

  LevelIntersectionIterator ilevelbegin() const {
    return LevelIntersectionIteratorWrapper<GridImp>(_grid,
                                                     _multiDomainEntityPointer->level(),
                                                     _multiDomainEntityPointer->ilevelbegin());
  }

  LevelIntersectionIterator ilevelend() const {
    return LevelIntersectionIteratorWrapper<GridImp>(_grid,
                                                     _multiDomainEntityPointer->level(),
                                                     _multiDomainEntityPointer->ilevelend());
  }

  LeafIntersectionIterator ileafbegin() const {
    return LeafIntersectionIteratorWrapper<GridImp>(_grid,
                                                     _multiDomainEntityPointer->ileafbegin());
  }

  LeafIntersectionIterator ileafend() const {
    return LeafIntersectionIteratorWrapper<GridImp>(_grid,
                                                    _multiDomainEntityPointer->ileafend());
  }


  bool isNew() const {
    return _multiDomainEntityPointer->isNew();
  }

  bool mightVanish() const {
    return _multiDomainEntityPointer->mightVanish();
  }

  EntitySeed seed() const {
    return EntitySeed(hostEntity().seed());
  }

  bool hasBoundaryIntersections () const
  {
    typedef typename GridImp::template Codim<0>::Entity Entity;
    Entity entity(*this);

    {
      typedef typename GridImp::LevelIntersectionIterator IntersectionIterator;
      IntersectionIterator end = _grid.levelGridView(level()).iend(entity);
      for(IntersectionIterator it = _grid.levelGridView(level()).ibegin(entity); it != end; ++it)
      {
        if( it->boundary() ) return true;
      }
    }

    {
      typedef typename GridImp::LeafIntersectionIterator IntersectionIterator;
      IntersectionIterator end = _grid.leafGridView().iend(entity);
      for(IntersectionIterator it = _grid.leafGridView().ibegin(entity); it != end; ++it)
      {
        if( it->boundary() ) return true;
      }
    }

    return false;
  }

  const MultiDomainEntityPointer& multiDomainEntityPointer() const {
    return _multiDomainEntityPointer;
  }

private:

  const GridImp& _grid;
  MultiDomainEntityPointer _multiDomainEntityPointer;

  const EntityWrapper& operator=(const EntityWrapper& rhs) {
    assert(_grid == rhs._grid);
    reset(rhs._multiDomainEntityPointer);
    return *this;
  }

  template<typename MultiDomainIteratorOrEntityPointer>
  void reset(const MultiDomainIteratorOrEntityPointer& multiDomainEntityPointer) {
    _multiDomainEntityPointer = multiDomainEntityPointer;
  }

  void compactify() {
    _multiDomainEntityPointer.compactify();
  }

  const HostEntity& hostEntity() const {
    return _grid._grid.hostEntity(*_multiDomainEntityPointer);
  }

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_ENTITY_HH
