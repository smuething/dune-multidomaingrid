#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_HH

#include <string>
#include <boost/shared_ptr.hpp>

#include <dune/common/exceptions.hh>
#include <dune/grid/multidomaingrid/subdomainset.hh>

#include <dune/grid/multidomaingrid/subdomaingrid/geometry.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/entity.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/entitypointer.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/leafiterator.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/leveliterator.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/hierarchiciterator.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/intersectioniterator.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/idsets.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/indexsets.hh>


namespace Dune {

namespace mdgrid {


// forward declaration in correct namespace
template<typename,typename>
class MultiDomainGrid;

namespace subdomain {

template<typename MDGrid>
class SubDomainGrid;

template<typename MDGrid>
struct SubDomainGridFamily {

  template <int dim, int dimw, class GridImp,
            template<int,int,class> class GeometryImp,
            template<int,int,class> class EntityImp,
            template<int,class> class EntityPointerImp,
            template<int,PartitionIteratorType,class> class LevelIteratorImp,
            template<class> class LeafIntersectionImp,
            template<class> class LevelIntersectionImp,
            template<class> class LeafIntersectionIteratorImp,
            template<class> class LevelIntersectionIteratorImp,
            template<class> class HierarchicIteratorImp,
            template<int,PartitionIteratorType,class> class LeafIteratorImp,
            class LevelIndexSetImp, class LeafIndexSetImp,
            class GlobalIdSetImp, class GIDType, class LocalIdSetImp, class LIDType, class CCType,
            template<class,PartitionIteratorType> class LevelGridViewTraits = DefaultLevelGridViewTraits,
            template<class,PartitionIteratorType> class LeafGridViewTraits = DefaultLeafGridViewTraits
            >
  struct SubDomainGridTraits
  {
    /** \brief The type that implementing the grid. */
    typedef GridImp Grid;

    /** \brief The type of the intersection at the leafs of the grid. */
    typedef Dune::Intersection<const GridImp, LeafIntersectionImp>  LeafIntersection;
    /** \brief The type of the intersection at the levels of the grid. */
    typedef Dune::Intersection<const GridImp, LevelIntersectionImp> LevelIntersection;
    /** \brief The type of the intersection iterator at the leafs of the grid. */
    typedef Dune::IntersectionIterator<const GridImp, LeafIntersectionIteratorImp, LeafIntersectionImp>   LeafIntersectionIterator;
    /** \brief The type of the intersection iterator at the levels of the grid. */
    typedef Dune::IntersectionIterator<const GridImp, LevelIntersectionIteratorImp, LevelIntersectionImp> LevelIntersectionIterator;

    /** \brief The type of the  hierarchic iterator. */
    typedef Dune::HierarchicIterator<const GridImp, HierarchicIteratorImp> HierarchicIterator;

    /**
     * \brief Traits associated with a specific codim.
     * \tparam cd The codimension.
     */
    template <int cd>
    struct Codim
    {
      //! IMPORTANT: Codim<codim>::Geometry == Geometry<dim-codim,dimw>
      /** \brief The type of the geometry associated with the entity.*/
      typedef Dune::Geometry<dim-cd, dimw, const GridImp, GeometryImp> Geometry;
      /** \brief The type of the local geometry associated with the entity.*/
      typedef Dune::Geometry<dim-cd, dim, const GridImp, GeometryImp> LocalGeometry;
      /** \brief The type of the entity. */
      // we could - if needed - introduce another struct for dimglobal of Geometry
      typedef Dune::Entity<cd, dim, const GridImp, EntityImp> Entity;

      /** \brief The type of the iterator over all level entities of this codim. */
      typedef Dune::LevelIterator<cd,All_Partition,const GridImp,LevelIteratorImp> LevelIterator;

      /** \brief The type of the iterator over all leaf entities of this codim. */
      typedef Dune::LeafIterator<cd,All_Partition,const GridImp,LeafIteratorImp> LeafIterator;

      /** \brief The type of the entity pointer for entities of this codim.*/
      typedef Dune::EntityPointer<const GridImp,EntityPointerImp<cd,const GridImp> > EntityPointer;

      /**
       * \brief Traits associated with a specific grid partition type.
       * \tparam pitype The type of the grid partition.
       */
      template <PartitionIteratorType pitype>
      struct Partition
      {
        /** \brief The type of the iterator over the level entities of this codim on this partition. */
        typedef Dune::LevelIterator<cd,pitype,const GridImp,LevelIteratorImp> LevelIterator;
        /** \brief The type of the iterator over the leaf entities of this codim on this partition. */
        typedef Dune::LeafIterator<cd,pitype,const GridImp,LeafIteratorImp> LeafIterator;
      };
    private:
      friend class Dune::Entity<cd, dim, const GridImp, EntityImp>;




      typedef EntityPointerImp<cd,const GridImp> EntityPointerImpl;
    };

    /**
     * \brief Traits associated with a specific grid partition type.
     * \tparam pitype The type of the grid partition.
     */
    template <PartitionIteratorType pitype>
    struct Partition
    {
      /** \brief The type of the level grid view associated with this partition type. */
      typedef Dune::GridView<LevelGridViewTraits<const GridImp,pitype> >
      LevelGridView;

      /** \brief The type of the leaf grid view associated with this partition type. */
      typedef Dune::GridView<LeafGridViewTraits<const GridImp,pitype> >
      LeafGridView;
    };

    /** \brief The type of the level index set. */
    typedef LevelIndexSetImp LevelIndexSet;
    /** \brief The type of the leaf index set. */
    typedef LeafIndexSetImp LeafIndexSet;
    /** \brief The type of the global id set. */
    typedef IdSet<const GridImp,GlobalIdSetImp,GIDType> GlobalIdSet;
    /** \brief The type of the local id set. */
    typedef IdSet<const GridImp,LocalIdSetImp,LIDType> LocalIdSet;

    /** \brief The type of the collective communication. */
    typedef CCType CollectiveCommunication;
  };

  typedef SubDomainGridTraits<
    MDGrid::dimension,
    MDGrid::dimensionworld,
    SubDomainGrid<MDGrid>,
    GeometryWrapper,
    EntityWrapper,
    EntityPointerWrapper,
    LevelIteratorWrapper,
    LeafIntersectionIteratorWrapper, // leaf intersection
    LevelIntersectionIteratorWrapper, // level intersection
    LeafIntersectionIteratorWrapper, // leaf intersection iterator
    LevelIntersectionIteratorWrapper, // level intersection iterator
    HierarchicIteratorWrapper,
    LeafIteratorWrapper,
    IndexSetWrapper<const SubDomainGrid<MDGrid>, typename MDGrid::LevelIndexSetImp>,
    IndexSetWrapper<const SubDomainGrid<MDGrid>, typename MDGrid::LeafIndexSetImp>,
    IdSetWrapper<const SubDomainGrid<MDGrid>, typename MDGrid::HostGridType::Traits::GlobalIdSet>,
    typename MDGrid::Traits::GlobalIdSet::IdType,
    IdSetWrapper<const SubDomainGrid<MDGrid>, typename MDGrid::HostGridType::Traits::LocalIdSet>,
    typename MDGrid::Traits::LocalIdSet::IdType,
    typename MDGrid::HostGridType::CollectiveCommunication
    > Traits;

};


template<typename MDGrid>
class SubDomainGrid :
    public GridDefaultImplementation<MDGrid::dimension,
				     MDGrid::dimensionworld,
				     typename MDGrid::ctype,
				     SubDomainGridFamily<MDGrid> > {

  template<typename,typename>
  friend class ::Dune::mdgrid::MultiDomainGrid;

  template<int codim, int dim, typename GridImp>
  friend class EntityWrapper;

  template<int codim, int dim, typename GridImp>
  friend class MakeableEntityWrapper;

  template<int codim, typename GridImp>
  friend class EntityPointerWrapper;

  template<int codim, PartitionIteratorType pitype, typename GridImp>
  friend class LeafIteratorWrapper;

  template<int codim, PartitionIteratorType pitype, typename GridImp>
  friend class LevelIteratorWrapper;

  template<typename GridImp>
  friend class HierarchicIteratorWrapper;

  template<int mydim, int coorddim, typename GridImp>
  friend class GeometryWrapper;

  template<int mydim, int coorddim, typename GridImp>
  friend class MakeableGeometryWrapper;

  template<typename GridImp, typename WrappedIndexSet>
  friend class IndexSetWrapper;

  template<typename GridImp, typename WrappedIdSet>
  friend class IdSetWrapper;

  template<typename GridImp>
  friend struct ::Dune::mdgrid::detail::HostGridAccessor;

  template<typename GridImp>
  friend class LeafIntersectionIteratorWrapper;

  template<typename GridImp>
  friend class LevelIntersectionIteratorWrapper;

  typedef SubDomainGrid<MDGrid> GridImp;
  typedef MDGrid MDGridType;

  typedef typename MDGrid::HostGridType HostGridType;

  typedef IndexSetWrapper<const SubDomainGrid<MDGrid>, typename MDGrid::LevelIndexSetImp> LevelIndexSetImp;

  typedef IndexSetWrapper<const SubDomainGrid<MDGrid>, typename MDGrid::LeafIndexSetImp> LeafIndexSetImp;

  typedef IdSetWrapper<const SubDomainGrid<MDGrid>, typename HostGridType::Traits::GlobalIdSet> GlobalIdSetImp;

  typedef IdSetWrapper<const SubDomainGrid<MDGrid>, typename HostGridType::Traits::LocalIdSet> LocalIdSetImp;

public:

  typedef SubDomainGridFamily<MDGrid> GridFamily;
  typedef typename GridFamily::Traits Traits;
  typedef typename MDGrid::ctype ctype;
  typedef typename MDGrid::SubDomainType SubDomainType;
  typedef MDGridType MultiDomainGrid;

public:

  std::string name() const {
    return "SubDomainGrid";
  }

  std::size_t maxLevel() const {
    return _grid.maxLevel();
  }

  template<int codim>
  typename Traits::template Codim<codim>::LevelIterator lbegin(int level) const {
    return LevelIteratorWrapper<codim,All_Partition,const GridImp>(*_levelIndexSets[level],
                                                                   _grid.template lbegin<codim>(level),
                                                                   _grid.template lend<codim>(level));
  }

  template<int codim>
  typename Traits::template Codim<codim>::LevelIterator lend(int level) const {
    return LevelIteratorWrapper<codim,All_Partition,const GridImp>(*_levelIndexSets[level],
                                                                   _grid.template lend<codim>(level),
                                                                   _grid.template lend<codim>(level));
  }

  template<int codim, PartitionIteratorType PiType>
  typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin(int level) const {
    return LevelIteratorWrapper<codim,PiType,const GridImp>(*_levelIndexSets[level],
                                                            _grid.template lbegin<codim,PiType>(level),
                                                            _grid.template lend<codim,PiType>(level));
  }

  template<int codim, PartitionIteratorType PiType>
  typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend(int level) const {
    return LevelIteratorWrapper<codim,PiType,const GridImp>(*_levelIndexSets[level],
                                                            _grid.template lend<codim,PiType>(level),
                                                            _grid.template lend<codim,PiType>(level));
  }

  template<int codim>
  typename Traits::template Codim<codim>::LeafIterator leafbegin() const {
    return LeafIteratorWrapper<codim,All_Partition,const GridImp>(_leafIndexSet,
                                                                  _grid.template leafbegin<codim>(),
                                                                  _grid.template leafend<codim>());
  }

  template<int codim>
  typename Traits::template Codim<codim>::LeafIterator leafend() const {
    return LeafIteratorWrapper<codim,All_Partition,const GridImp>(_leafIndexSet,
                                                                  _grid.template leafend<codim>(),
                                                                  _grid.template leafend<codim>());
  }

  template<int codim, PartitionIteratorType PiType>
  typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const {
    return LeafIteratorWrapper<codim,PiType,const GridImp>(_leafIndexSet,
                                                           _grid.template leafbegin<codim,PiType>(),
                                                           _grid.template leafend<codim,PiType>());
  }

  template<int codim, PartitionIteratorType PiType>
  typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const {
    return LeafIteratorWrapper<codim,PiType,const GridImp>(_leafIndexSet,
                                                           _grid.template leafend<codim,PiType>(),
                                                           _grid.template leafend<codim,PiType>());
  }

  int size(int level, int codim) const {
    assert(level <= maxLevel());
    return _levelIndexSets[level]->size(codim);
  }

  int size(int codim) const {
    return _leafIndexSet.size(codim);
  }

  int size(int level, GeometryType type) const {
    assert(level <= maxLevel());
    return _levelIndexSets[level]->size(type);
  }

  int size(GeometryType type) const {
    // TODO: check this (perhaps return sum over levelindexsets?)
    return _leafIndexSet.size(type);
  }

  const typename Traits::GlobalIdSet& globalIdSet() const {
    return _globalIdSet;
  }

  const typename Traits::LocalIdSet& localIdSet() const {
    return _localIdSet;
  }

  const typename Traits::LevelIndexSet& levelIndexSet(int level) const {
    if (!_grid.supportLevelIndexSets()) {
      DUNE_THROW(GridError,"level index set support not enabled for this grid");
    }
    assert(level <= maxLevel());
    return *_levelIndexSets[level];
  }

  const typename Traits::LeafIndexSet& leafIndexSet() const {
    return _leafIndexSet;
  }

  void globalRefine(int refCount) {
    DUNE_THROW(NotImplemented,"grid modification only allowed on the MultiDomainGrid");
  }

  bool mark(int refCount, const typename Traits::template Codim<0>::Entity& e) {
    DUNE_THROW(NotImplemented,"grid modification only allowed on the MultiDomainGrid");
  }

  int getMark(const typename Traits::template Codim<0>::Entity& e) {
    DUNE_THROW(NotImplemented,"grid modification only allowed on the MultiDomainGrid");
  }

  bool preAdapt() {
    DUNE_THROW(NotImplemented,"grid modification only allowed on the MultiDomainGrid");
  }

  bool adapt() {
    DUNE_THROW(NotImplemented,"grid modification only allowed on the MultiDomainGrid");
  }

  void postAdapt() {
    DUNE_THROW(NotImplemented,"grid modification only allowed on the MultiDomainGrid");
  }

  int overlapSize(int level, int codim) const {
    return _grid.overlapSize(level,codim);
  }

  int overlapSize(int codim) const {
    return _grid.overlapSize(codim);
  }

  int ghostSize(int level, int codim) const {
    return _grid.ghostSize(level,codim);
  }

  int ghostSize(int codim) const {
    return _grid.ghostSize(codim);
  }

  const typename Traits::CollectiveCommunication& comm() const {
    return _grid.comm();
  }

  const MDGrid& multiDomainGrid() const {
    return _grid;
  }

  SubDomainType domain() const {
    return _subDomain;
  }

  void update() const {
    if (_grid.supportLevelIndexSets()) {
      while (_levelIndexSets.size() <= maxLevel()) {
        _levelIndexSets.push_back(make_shared_ptr(new LevelIndexSetImp(*this,_grid.levelIndexSet(_levelIndexSets.size()))));
      }
    }
  }

  bool operator==(const SubDomainGrid& rhs) const {
    return (&_grid == &rhs._grid && _subDomain == rhs._subDomain);
  }

  template<int cc>
  typename Traits::template Codim<cc>::EntityPointer subDomainEntityPointer(const typename MDGrid::Traits::template Codim<cc>::Entity& mdEntity) const {
    return EntityPointerWrapper<cc,const SubDomainGrid<MDGrid> >(*this,typename MDGrid::Traits::template Codim<cc>::EntityPointer(mdEntity));
  }

  template<typename EntityType>
  typename Traits::template Codim<EntityType::codimension>::EntityPointer subDomainEntityPointer(const EntityType& mdEntity) const {
    return subDomainEntityPointer<EntityType::codimension>(mdEntity);
  }

  template<typename EntityType>
  const typename MDGrid::template MultiDomainEntity<EntityType>::type& multiDomainEntity(const EntityType& e) const {
    return *(getRealImplementation(e).multiDomainEntityPointer());
  }

  template<typename EntityType>
  typename MDGrid::template MultiDomainEntityPointer<EntityType>::type multiDomainEntityPointer(const EntityType& e) const {
    return getRealImplementation(e).multiDomainEntityPointer();
  }

  template<typename EntityType>
  const typename MDGrid::template HostEntity<EntityType>::type& hostEntity(const EntityType& e) const {
    return *(getRealImplementation(e).hostEntityPointer());
  }

  template<typename EntityType>
  typename MDGrid::template HostEntityPointer<EntityType>::type hostEntityPointer(const EntityType& e) const {
    return getRealImplementation(e).hostEntityPointer();
  }

  typename Traits::LeafIntersectionIterator getSubDomainInterfaceIterator(const typename MDGrid::LeafSubDomainInterfaceIterator it) const {
    assert(_subDomain == it.domain1() || _subDomain == it.domain2());
    if (_subDomain == it.domain1())
      return LeafIntersectionIteratorWrapper<const GridImp>(*this,it.firstHostIntersection());
    else
      return LeafIntersectionIteratorWrapper<const GridImp>(*this,it.secondHostIntersection());
  }

  typename Traits::LevelIntersectionIterator getSubDomainInterfaceIterator(const typename MDGrid::LevelSubDomainInterfaceIterator it) const {
    assert(_subDomain == it.domain1() || _subDomain == it.domain2());
    if (_subDomain == it.domain1())
      return LevelIntersectionIteratorWrapper<const GridImp>(*this,it.firstHostIntersection());
    else
      return LevelIntersectionIteratorWrapper<const GridImp>(*this,it.secondHostIntersection());
  }

private:

  const MDGrid& _grid;
  SubDomainType _subDomain;
  GlobalIdSetImp _globalIdSet;
  LocalIdSetImp _localIdSet;
  LeafIndexSetImp _leafIndexSet;
  mutable std::vector<boost::shared_ptr<LevelIndexSetImp> > _levelIndexSets;

  SubDomainGrid(const MDGrid& grid, SubDomainType subDomain) :
    _grid(grid),
    _subDomain(subDomain),
    _globalIdSet(*this,grid._hostGrid.globalIdSet()),
    _localIdSet(*this,grid._hostGrid.localIdSet()),
    _leafIndexSet(*this,grid.leafIndexSet())
  {
    update();
  }
  /*
  void reset(const SubDomainGrid& rhs) {
    assert(_grid == rhs._grid);
    _subDomain = rhs._subDomain;
  }

  void reset(const SubDomainType subDomain) {
    _subDomain = subDomain;
  }
  */

  template<typename EntityType>
  bool containsMultiDomainEntity(const EntityType& e) const {
    return levelIndexSet(e.level()).containsMultiDomainEntity(e);
  }

  template<typename EntityType>
  bool containsHostEntity(const EntityType& e) const {
    return levelIndexSet(e.level()).containsHostEntity(e);
  }

  SubDomainGrid(const SubDomainGrid& rv);
  SubDomainGrid& operator=(const SubDomainGrid& rv);

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_MULTIDOMAINGRID_HH
