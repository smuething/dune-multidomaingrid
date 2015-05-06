#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_ENTITY_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_ENTITY_HH

#include <utility>

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
class EntityWrapperBase
  : public EntityDefaultImplementation<codim,dim,GridImp,EntityWrapper>
{

protected:

  using MultiDomainEntity = typename GridImp::MultiDomainGrid::Traits::template Codim<codim>::Entity;
  using HostEntity        = typename GridImp::HostGridType::Traits::template Codim<codim>::Entity;

public:

  using EntitySeed    = EntitySeedWrapper<typename HostEntity::EntitySeed>;
  using Geometry      = typename GridImp::template Codim<codim>::Geometry;
  using EntityPointer = typename GridImp::Traits::template Codim<0>::EntityPointer;

  EntityWrapperBase()
    : _grid(nullptr)
  {}

  EntityWrapperBase(const GridImp* grid, const MultiDomainEntity& e)
    : _grid(grid)
    , _multiDomainEntity(e)
  {}

  EntityWrapperBase(const GridImp* grid, MultiDomainEntity&& e)
    : _grid(grid)
    , _multiDomainEntity(std::move(e))
  {}

  int level() const {
    return multiDomainEntity().level();
  }

  PartitionType partitionType() const {
    return multiDomainEntity().partitionType();
  }

  template<int cc>
  int DUNE_DEPRECATED_MSG("Use subEntities instead") count() const {
    return multiDomainEntity().template count<cc>();
  }

  unsigned int subEntities(unsigned int codimSubEntitiy) const {
    return multiDomainEntity().subEntities(codimSubEntitiy);
  }

  Geometry geometry() const {
    return Geometry(hostEntity().geometry());
  }

  EntitySeed seed() const {
    return EntitySeed(hostEntity().seed());
  }

  bool equals(const EntityWrapperBase& other) const
  {
    return grid() == other.grid() && multiDomainEntity() == other.multiDomainEntity();
  }

private:

  const GridImp* _grid;
  MultiDomainEntity _multiDomainEntity;

public:

  const MultiDomainEntity& multiDomainEntity() const {
    return _multiDomainEntity;
  }

  const HostEntity& hostEntity() const {
    return grid()._grid.hostEntity(_multiDomainEntity);
  }

  const GridImp& grid() const
  {
    assert(_grid);
    return *_grid;
  }

};


template<int codim, int dim, typename GridImp>
class EntityWrapper
  : public EntityWrapperBase<codim,dim,GridImp>
{

  using Base = EntityWrapperBase<codim,dim,GridImp>;

  template<int, typename>
  friend class EntityPointerWrapper;

  template<typename>
  friend class SubDomainGrid;

  template<typename,typename>
  friend class Dune::mdgrid::MultiDomainGrid;

  template<int, int, typename, template<int,int,typename> class>
  friend class Entity;

public:

  // inherit constructors
  using Base::Base;

};


template<int dim, typename GridImp>
class EntityWrapper<0,dim,GridImp>
  : public EntityWrapperBase<0,dim,GridImp>
{

  using Base = EntityWrapperBase<0,dim,GridImp>;

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

  using Base::multiDomainEntity;
  using Base::hostEntity;
  using Base::grid;

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
    return {EntityWrapper<cc,dim,GridImp>(&grid(),multiDomainEntity().template subEntity<cc>(i))};
  }

  typename GridImp::template Codim<0>::Entity father() const {
    return {EntityWrapper<0,dim,GridImp>(&grid(),multiDomainEntity().father())};
  }

  bool hasFather() const {
    return multiDomainEntity().hasFather();
  }

  bool isLeaf() const {
    return multiDomainEntity().isLeaf();
  }

  bool isRegular() const {
    return multiDomainEntity().isRegular();
  }

  LocalGeometry geometryInFather() const {
    return LocalGeometry(hostEntity().geometryInFather());
  }

  HierarchicIterator hbegin(int maxLevel) const {
    return HierarchicIteratorWrapper<GridImp>(
      grid(),
      multiDomainEntity().hbegin(maxLevel),
      multiDomainEntity().hend(maxLevel)
      );
  }

  HierarchicIterator hend(int maxLevel) const {
    return HierarchicIteratorWrapper<GridImp>(
      &grid(),
      multiDomainEntity().hend(maxLevel),
      multiDomainEntity().hend(maxLevel)
      );
  }

  LevelIntersectionIterator ilevelbegin() const {
    return LevelIntersectionIteratorWrapper<GridImp>(
      grid(),
      multiDomainEntity.level(),
      multiDomainEntity.ilevelbegin()
      );
  }

  LevelIntersectionIterator ilevelend() const {
    return LevelIntersectionIteratorWrapper<GridImp>(
      grid(),
      multiDomainEntity.level(),
      multiDomainEntity.ilevelend()
      );
  }

  LeafIntersectionIterator ileafbegin() const {
    return LeafIntersectionIteratorWrapper<GridImp>(
      grid(),
      multiDomainEntity.ileafbegin()
      );
  }

  LeafIntersectionIterator ileafend() const {
    return LeafIntersectionIteratorWrapper<GridImp>(
      grid(),
      multiDomainEntity.ileafend()
      );
  }

  bool isNew() const {
    return multiDomainEntity().isNew();
  }

  bool mightVanish() const {
    return multiDomainEntity().mightVanish();
  }

  bool hasBoundaryIntersections () const
  {
    for(const auto& is : intersections(grid().levelGridView(this->level()),*this))
      {
        if(is.boundary())
          return true;
      }

    if (isLeaf())
      for(const auto& is : intersections(grid().leafGridView(),*this))
        {
          if(is.boundary())
            return true;
        }

    return false;
  }

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_ENTITY_HH
