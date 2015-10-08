#ifndef DUNE_MULTIDOMAINGRID_ENTITY_HH
#define DUNE_MULTIDOMAINGRID_ENTITY_HH

#include <dune/grid/common/entity.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/common/version.hh>

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

template<typename HostGrid, typename MDGridTraits>
class MultiDomainGrid;

template<typename, PartitionIteratorType>
class LevelGridView;

template<typename, PartitionIteratorType>
class LeafGridView;


template<typename HostES>
class EntitySeedWrapper
{

  typedef HostES HostEntitySeed;

  template<typename, typename>
  friend class MultiDomainGrid;

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

  template<typename,typename>
  friend class MultiDomainGrid;

  template<typename>
  friend class HierarchicIteratorWrapper;

  template<typename, PartitionIteratorType>
  friend class LevelGridView;

  template<typename, PartitionIteratorType>
  friend class LeafGridView;


  typedef typename GridImp::HostGridType::Traits::template Codim<codim>::EntityPointer HostEntityPointer;

  template<typename HostIteratorOrEntityPointer>
  explicit MakeableEntityWrapper(const HostIteratorOrEntityPointer& hostEntityPointer) :
    GridImp::template Codim<codim>::Entity(EntityWrapper<codim,dim,const GridImp>(hostEntityPointer))
  {}

  MakeableEntityWrapper& operator=(const MakeableEntityWrapper& rhs)
  {
    reset(rhs.hostEntityPointer());
    return *this;
  }

  template<typename HostIteratorOrEntityPointer>
  void reset(const HostIteratorOrEntityPointer& hostEntityPointer) {
    GridImp::getRealImplementation(*this).reset(hostEntityPointer);
  }

  void compactify() {
    GridImp::getRealImplementation(*this).compactify();
  }

  const HostEntityPointer& hostEntityPointer() const {
    return GridImp::getRealImplementation(*this).hostEntityPointer();
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

  template<typename,typename>
  friend class MultiDomainGrid;

  typedef typename GridImp::HostGridType::Traits::template Codim<codim>::EntityPointer HostEntityPointer;
  typedef typename GridImp::HostGridType::Traits::template Codim<codim>::Entity HostEntity;

public:

  typedef typename GridImp::template Codim<codim>::Geometry Geometry;

  typedef EntitySeedWrapper<typename HostEntity::EntitySeed> EntitySeed;

  template<typename HostIteratorOrEntityPointer>
  explicit EntityWrapper(const HostIteratorOrEntityPointer& e) :
    _hostEntityPointer(e)
  {}

  explicit EntityWrapper(const HostEntity& e) :
    _hostEntityPointer(e)
  {}

  // copy constructor. The default constructor does not work correctly here,
  // as it will initialise the geometry with a pointer to the host geometry
  // of rhs, causing geometry access to fail if rhs gets destructed.
  EntityWrapper(const EntityWrapper& rhs) :
    _hostEntityPointer(rhs._hostEntityPointer)
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

#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
  int subEntities(int cc) const {
    return _hostEntityPointer->subEntities(cc);
  }
#endif

  Geometry geometry() const {
    return Geometry(_hostEntityPointer->geometry());
  }

  EntitySeed seed() const {
    return EntitySeed(_hostEntityPointer->seed());
  }

  private:

  HostEntityPointer _hostEntityPointer;

  template<typename HostIteratorOrEntityPointer>
  void reset(const HostIteratorOrEntityPointer& hostEntityPointer) {
    _hostEntityPointer = hostEntityPointer;
  }

  void compactify() {
    _hostEntityPointer.compactify();
  }

  const HostEntityPointer& hostEntityPointer() const {
    return _hostEntityPointer;
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

  template<typename,typename>
  friend class MultiDomainGrid;

  template<typename, PartitionIteratorType>
  friend class LevelGridView;

  template<typename, PartitionIteratorType>
  friend class LeafGridView;


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


  template<typename HostIteratorOrEntityPointer>
  explicit EntityWrapper(const HostIteratorOrEntityPointer& e) :
    _hostEntityPointer(e)
  {}

  // copy constructor. The default constructor does not work correctly here,
  // as it will initialise the geometry with a pointer to the host geometry
  // of rhs, causing geometry access to fail if rhs gets destructed.
  EntityWrapper(const EntityWrapper& rhs) :
    _hostEntityPointer(rhs._hostEntityPointer)
  {}

  int level() const {
    return _hostEntityPointer->level();
  }

  PartitionType partitionType() const {
    return _hostEntityPointer->partitionType();
  }

  Geometry geometry() const {
    return Geometry(_hostEntityPointer->geometry());
  }

  template<int cc>
  int count() const {
    return _hostEntityPointer->template count<cc>();
  }

#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
  int subEntities(int cc) const {
    return _hostEntityPointer->subEntities(cc);
  }
#endif

  template<int cc>
  typename GridImp::template Codim<cc>::EntityPointer subEntity(int i) const {
    return EntityPointerWrapper<cc,GridImp>(_hostEntityPointer->template subEntity<cc>(i));
  }

  EntityPointer father() const {
    return EntityPointerWrapper<0,GridImp>(_hostEntityPointer->father());
  }

  bool hasFather() const {
    return _hostEntityPointer->hasFather();
  }

  bool isLeaf() const {
    return _hostEntityPointer->isLeaf();
  }

  bool isRegular() const {
    return _hostEntityPointer->isRegular();
  }

  LocalGeometry geometryInFather() const {
    return LocalGeometry(_hostEntityPointer->geometryInFather());
  }

  HierarchicIterator hbegin(int maxLevel) const {
    return HierarchicIteratorWrapper<GridImp>(_hostEntityPointer->hbegin(maxLevel));
  }

  HierarchicIterator hend(int maxLevel) const {
    return HierarchicIteratorWrapper<GridImp>(_hostEntityPointer->hend(maxLevel));
  }

  LevelIntersectionIterator ilevelbegin() const {
    return LevelIntersectionIteratorWrapper<GridImp>(_hostEntityPointer->ilevelbegin());
  }

  LevelIntersectionIterator ilevelend() const {
    return LevelIntersectionIteratorWrapper<GridImp>(_hostEntityPointer->ilevelend());
  }

  LeafIntersectionIterator ileafbegin() const {
    return LeafIntersectionIteratorWrapper<GridImp>(_hostEntityPointer->ileafbegin());
  }

  LeafIntersectionIterator ileafend() const {
    return LeafIntersectionIteratorWrapper<GridImp>(_hostEntityPointer->ileafend());
  }

  bool isNew() const {
    return _hostEntityPointer->isNew();
  }

  bool mightVanish() const {
    return _hostEntityPointer->mightVanish();
  }

  EntitySeed seed() const {
    return EntitySeed(_hostEntityPointer->seed());
  }

  bool hasBoundaryIntersections () const
  {
    return _hostEntityPointer->hasBoundaryIntersections();
  }

private:

  HostEntityPointer _hostEntityPointer;

  template<typename HostIteratorOrEntityPointer>
  void reset(const HostIteratorOrEntityPointer& hostEntityPointer) {
    _hostEntityPointer = hostEntityPointer;
  }

  void compactify() {
    _hostEntityPointer.compactify();
  }

  const HostEntityPointer& hostEntityPointer() const {
    return _hostEntityPointer;
  }

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_ENTITY_HH
