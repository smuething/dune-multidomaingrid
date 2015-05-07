#ifndef DUNE_MULTIDOMAINGRID_ENTITY_HH
#define DUNE_MULTIDOMAINGRID_ENTITY_HH

#include <dune/common/deprecated.hh>

#include <dune/grid/common/entity.hh>
#include <dune/grid/common/gridenums.hh>

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

template<typename, typename>
class IntersectionIteratorWrapper;

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
class EntityWrapperBase :
    public EntityDefaultImplementation<codim,dim,GridImp,EntityWrapper>
{

  template<int, typename>
  friend class EntityPointerWrapper;

  template<typename,typename>
  friend class MultiDomainGrid;

  typedef typename GridImp::HostGrid::Traits::template Codim<codim>::EntityPointer HostEntityPointer;
  typedef typename GridImp::HostGrid::Traits::template Codim<codim>::Entity HostEntity;

public:

  typedef typename GridImp::template Codim<codim>::Geometry Geometry;

  typedef EntitySeedWrapper<typename HostEntity::EntitySeed> EntitySeed;

  EntityWrapperBase()
  {}

  explicit EntityWrapperBase(const HostEntity& e)
    : _hostEntity(e)
  {}

  explicit EntityWrapperBase(HostEntity&& e)
    : _hostEntity(std::move(e))
  {}

  int level() const {
    return _hostEntity.level();
  }

  PartitionType partitionType() const {
    return _hostEntity.partitionType();
  }

  template<int cc>
  int DUNE_DEPRECATED_MSG("Use subEntities instead") count() const {
    return _hostEntity.template count<cc>();
  }

  unsigned int subEntities(unsigned int codimSubEntitiy) const {
    return _hostEntity.subEntities(codimSubEntitiy);
  }

  Geometry geometry() const {
    return Geometry(_hostEntity.geometry());
  }

  EntitySeed seed() const {
    return EntitySeed(_hostEntity.seed());
  }

  bool equals(const EntityWrapperBase& other) const
  {
    return hostEntity() == other.hostEntity();
  }

private:

  HostEntity _hostEntity;

protected:

  const HostEntity& hostEntity() const {
    return _hostEntity;
  }

};


template<int codim, int dim, typename GridImp>
class EntityWrapper :
    public EntityWrapperBase<codim,dim,GridImp>
{

  using Base = EntityWrapperBase<codim,dim,GridImp>;

  template<int, typename>
  friend class EntityPointerWrapper;

  template<typename,typename>
  friend class MultiDomainGrid;

public:

  // inherit constructors
  using Base::Base;

};

template<int dim, typename GridImp>
class EntityWrapper<0,dim,GridImp> :
    public EntityWrapperBase<0,dim,GridImp>
{

  using Base = EntityWrapperBase<0,dim,GridImp>;

  template<int, typename>
  friend class EntityPointerWrapper;

  template<typename,typename>
  friend class MultiDomainGrid;

  template<typename, PartitionIteratorType>
  friend class LevelGridView;

  template<typename, PartitionIteratorType>
  friend class LeafGridView;

  using Base::hostEntity;

public:

  using LocalGeometry = typename GridImp::template Codim<0>::LocalGeometry;
  using LeafIntersectionIterator = typename GridImp::Traits::LeafIntersectionIterator;
  using LevelIntersectionIterator = typename GridImp::Traits::LevelIntersectionIterator;
  using HierarchicIterator = typename GridImp::Traits::HierarchicIterator;

  // inherit constructors
  using Base::Base;

  template<int cc>
  int DUNE_DEPRECATED_MSG("Use subEntities instead") count() const {
    return hostEntity().template count<cc>();
  }

  unsigned int subEntities(unsigned int codim) const {
    return hostEntity().subEntities(codim);
  }

  template<int cc>
  typename GridImp::template Codim<cc>::Entity subEntity(int i) const {
    return {EntityWrapper<cc,dim,GridImp>(hostEntity().template subEntity<cc>(i))};
  }

  typename GridImp::template Codim<0>::Entity father() const {
    return {EntityWrapper(hostEntity().father())};
  }

  bool hasFather() const {
    return hostEntity().hasFather();
  }

  bool isLeaf() const {
    return hostEntity().isLeaf();
  }

  bool isRegular() const {
    return hostEntity().isRegular();
  }

  LocalGeometry geometryInFather() const {
    return LocalGeometry(hostEntity().geometryInFather());
  }

  HierarchicIterator hbegin(int maxLevel) const {
    return HierarchicIteratorWrapper<GridImp>(hostEntity().hbegin(maxLevel));
  }

  HierarchicIterator hend(int maxLevel) const {
    return HierarchicIteratorWrapper<GridImp>(hostEntity().hend(maxLevel));
  }

  LevelIntersectionIterator ilevelbegin() const {
    return IntersectionIteratorWrapper<
      GridImp,
      typename GridImp::HostGridType::LevelGridView::IntersectionIterator
      >(
        hostEntity().ilevelbegin()
        );
  }

  LevelIntersectionIterator ilevelend() const {
    return IntersectionIteratorWrapper<
      GridImp,
      typename GridImp::HostGridType::LevelGridView::IntersectionIterator
      >(
        hostEntity().ilevelend()
        );
  }

  LeafIntersectionIterator ileafbegin() const {
    return IntersectionIteratorWrapper<
      GridImp,
      typename GridImp::HostGridType::LeafGridView::IntersectionIterator
      >(
        hostEntity().ileafbegin()
        );
  }

  LeafIntersectionIterator ileafend() const {
    return IntersectionIteratorWrapper<
      GridImp,
      typename GridImp::HostGridType::LeafGridView::IntersectionIterator
      >(
        hostEntity().ileafend()
        );
  }

  bool isNew() const {
    return hostEntity().isNew();
  }

  bool mightVanish() const {
    return hostEntity().mightVanish();
  }

  bool hasBoundaryIntersections () const
  {
    return hostEntity().hasBoundaryIntersections();
  }

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_ENTITY_HH
