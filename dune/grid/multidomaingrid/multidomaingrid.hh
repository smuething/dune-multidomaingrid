#ifndef DUNE_MULTIDOMAINGRID_MULTIDOMAINGRID_HH
#define DUNE_MULTIDOMAINGRID_MULTIDOMAINGRID_HH

#include <string>
#include <boost/shared_ptr.hpp>

namespace Dune {
namespace mdgrid {
namespace detail {

  template<typename GridImp>
  struct HostGridAccessor {
    typedef typename GridImp::HostGridType::Traits Traits;
    typedef typename GridImp::HostGridType Type;
  };

}
}
}

#include <dune/grid/multidomaingrid/subdomainset.hh>

#include <dune/grid/multidomaingrid/subdomaingrid/subdomaingrid.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/subdomaingridpointer.hh>

#include <dune/grid/multidomaingrid/geometry.hh>
#include <dune/grid/multidomaingrid/entity.hh>
#include <dune/grid/multidomaingrid/entitypointer.hh>
#include <dune/grid/multidomaingrid/leafiterator.hh>
#include <dune/grid/multidomaingrid/leveliterator.hh>
#include <dune/grid/multidomaingrid/hierarchiciterator.hh>
#include <dune/grid/multidomaingrid/intersectioniterator.hh>
#include <dune/grid/multidomaingrid/idsets.hh>
#include <dune/grid/multidomaingrid/indexsets.hh>

#include <dune/grid/multidomaingrid/subdomaininterfaceiterator.hh>


namespace Dune {

namespace mdgrid {

template<typename HostGrid>
class MultiDomainGrid;

template<typename HostGrid>
struct MultiDomainGridFamily {

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
  struct MultiDomainGridTraits
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

  typedef MultiDomainGridTraits<
    HostGrid::dimension,
    HostGrid::dimensionworld,
    MultiDomainGrid<HostGrid>,
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
    IndexSetWrapper<const MultiDomainGrid<HostGrid>, typename HostGrid::LevelGridView>,
    IndexSetWrapper<const MultiDomainGrid<HostGrid>, typename HostGrid::LeafGridView>,
    IdSetWrapper<const MultiDomainGrid<HostGrid>, typename HostGrid::Traits::GlobalIdSet>,
    typename HostGrid::Traits::GlobalIdSet::IdType,
    IdSetWrapper<const MultiDomainGrid<HostGrid>, typename HostGrid::Traits::LocalIdSet>,
    typename HostGrid::Traits::LocalIdSet::IdType,
    CollectiveCommunication<HostGrid>
    > Traits;

};

template<typename HostGrid>
class MultiDomainGrid :
    public GridDefaultImplementation<HostGrid::dimension,
				     HostGrid::dimensionworld,
				     typename HostGrid::ctype,
				     MultiDomainGridFamily<HostGrid> > {


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
  friend struct detail::HostGridAccessor;

  template<typename GridImp>
  friend class LeafIntersectionIteratorWrapper;

  template<typename GridImp>
  friend class LevelIntersectionIteratorWrapper;

  template<typename>
  friend class subdomain::SubDomainGrid;

  template<typename>
  friend struct subdomain::SubDomainGridFamily;

  template <typename>
  friend class LeafSubDomainInterfaceIterator;

  template <typename>
  friend class LevelSubDomainInterfaceIterator;

  typedef MultiDomainGrid<HostGrid> GridImp;
  typedef HostGrid HostGridType;

  typedef IndexSetWrapper<const MultiDomainGrid<HostGrid>, typename HostGrid::LevelGridView> LevelIndexSetImp;

  typedef IndexSetWrapper<const MultiDomainGrid<HostGrid>, typename HostGrid::LeafGridView> LeafIndexSetImp;

  typedef IdSetWrapper<const MultiDomainGrid<HostGrid>, typename HostGrid::Traits::GlobalIdSet> GlobalIdSetImp;

  typedef IdSetWrapper<const MultiDomainGrid<HostGrid>, typename HostGrid::Traits::LocalIdSet> LocalIdSetImp;

  enum State { fixed, marking, preUpdate, postUpdate };

  typedef MultiDomainGrid<HostGrid> ThisType;

public:

  typedef MultiDomainGridFamily<HostGrid> GridFamily;
  typedef typename GridFamily::Traits Traits;
  typedef typename HostGrid::ctype ctype;

  static const std::size_t maxNumberOfSubDomains = 3;

  typedef IntegralTypeSubDomainSet<maxNumberOfSubDomains> SubDomainSet;
  typedef typename SubDomainSet::DomainType SubDomainType;

  typedef subdomain::SubDomainGrid<ThisType> SubDomainGrid;
  typedef subdomain::SubDomainGridPointer<SubDomainGrid> SubDomainGridPointer;

  typedef LeafSubDomainInterfaceIterator<const ThisType> LeafSubDomainInterfaceIteratorType;
  typedef LevelSubDomainInterfaceIterator<const ThisType> LevelSubDomainInterfaceIteratorType;

  explicit MultiDomainGrid(HostGrid& hostGrid) :
    _hostGrid(hostGrid),
    _leafIndexSet(*this,hostGrid.leafView()),
    _globalIdSet(*this),
    _localIdSet(*this),
    _state(fixed),
    _subDomainGrid(*this,0)
  {
    updateIndexSets();
    _subDomainGrid.update();
  }

  std::string name() const {
    return "MultiDomainGrid";
  }

  std::size_t maxLevel() const {
    return _hostGrid.maxLevel();
  }

  template<int codim>
  typename Traits::template Codim<codim>::LevelIterator lbegin(int level) const {
    return LevelIteratorWrapper<codim,All_Partition,const GridImp>(_hostGrid.template lbegin<codim>(level));
  }

  template<int codim>
  typename Traits::template Codim<codim>::LevelIterator lend(int level) const {
    return LevelIteratorWrapper<codim,All_Partition,const GridImp>(_hostGrid.template lend<codim>(level));
  }

  template<int codim, PartitionIteratorType PiType>
  typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin(int level) const {
    return LevelIteratorWrapper<codim,PiType,const GridImp>(_hostGrid.template lbegin<codim,PiType>(level));
  }

  template<int codim, PartitionIteratorType PiType>
  typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend(int level) const {
    return LevelIteratorWrapper<codim,PiType,const GridImp>(_hostGrid.template lend<codim,PiType>(level));
  }

  template<int codim>
  typename Traits::template Codim<codim>::LeafIterator leafbegin() const {
    return LeafIteratorWrapper<codim,All_Partition,const GridImp>(_hostGrid.template leafbegin<codim>());
  }

  template<int codim>
  typename Traits::template Codim<codim>::LeafIterator leafend() const {
    return LeafIteratorWrapper<codim,All_Partition,const GridImp>(_hostGrid.template leafend<codim>());
  }

  template<int codim, PartitionIteratorType PiType>
  typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const {
    return LeafIteratorWrapper<codim,PiType,const GridImp>(_hostGrid.template leafbegin<codim,PiType>());
  }

  template<int codim, PartitionIteratorType PiType>
  typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const {
    return LeafIteratorWrapper<codim,PiType,const GridImp>(_hostGrid.template leafend<codim,PiType>());
  }

  LeafSubDomainInterfaceIteratorType leafSubDomainInterfaceBegin(SubDomainType subDomain1, SubDomainType subDomain2) const {
    return LeafSubDomainInterfaceIteratorType(*this,subDomain1,subDomain2);
  }

  LeafSubDomainInterfaceIteratorType leafSubDomainInterfaceEnd(SubDomainType subDomain1, SubDomainType subDomain2) const {
    return LeafSubDomainInterfaceIteratorType(*this,subDomain1,subDomain2,true);
  }

  LevelSubDomainInterfaceIteratorType levelSubDomainInterfaceBegin(SubDomainType subDomain1, SubDomainType subDomain2, int level) const {
    return LevelSubDomainInterfaceIteratorType(*this,subDomain1,subDomain2,level);
  }

  LevelSubDomainInterfaceIteratorType levelSubDomainInterfaceEnd(SubDomainType subDomain1, SubDomainType subDomain2, int level) const {
    return LevelSubDomainInterfaceIteratorType(*this,subDomain1,subDomain2,level,true);
  }

  int size(int level, int codim) const {
    return _hostGrid.size(level,codim);
  }

  int size(int codim) const {
    return _hostGrid.size(codim);
  }

  int size(int level, GeometryType type) const {
    return _hostGrid.size(level,type);
  }

  int size(GeometryType type) const {
    return _hostGrid.size(type);
  }

  const typename Traits::GlobalIdSet& globalIdSet() const {
    return _globalIdSet;
  }

  const typename Traits::LocalIdSet& localIdSet() const {
    return _localIdSet;
  }

  const typename Traits::LevelIndexSet& levelIndexSet(int level) const {
    assert(level <= maxLevel());
    return *_levelIndexSets[level];
  }

  const typename Traits::LeafIndexSet& leafIndexSet() const {
    return _leafIndexSet;
  }

  void globalRefine(int refCount) {
    _hostGrid.globalRefine(refCount);
    updateIndexSets();
  }

  bool mark(int refCount, const typename Traits::template Codim<0>::Entity& e) {
    return _hostGrid.mark(refCount, hostEntity(e));
  }

  int getMark(const typename Traits::template Codim<0>::Entity& e) {
    return _hostGrid.getMark(hostEntity(e));
  }

  bool preAdapt() {
    assert(_state == fixed);
    return _hostGrid.preAdapt();
  }

  bool adapt() {
    assert(_state == fixed);
    bool r = _hostGrid.adapt();
    updateIndexSets();
    return r;
  }

  void postAdapt() {
    assert(_state == fixed);
    _hostGrid.postAdapt();
  }

  int overlapSize(int level, int codim) const {
    return _hostGrid.overlapSize(level,codim);
  }

  int overlapSize(int codim) const {
    return _hostGrid.overlapSize(codim);
  }

  int ghostSize(int level, int codim) const {
    return _hostGrid.ghostSize(level,codim);
  }

  int ghostSize(int codim) const {
    return _hostGrid.ghostSize(codim);
  }

  const typename Traits::CollectiveCommunication& comm() const {
    return _hostGrid.comm();
  }

  void startSubDomainMarking() {
    assert(_state == fixed);
    _tmpLeafIndexSet.reset(new LeafIndexSetImp(_leafIndexSet));
    _state = marking;
  }

  void preUpdateSubDomains() {
    assert(_state == marking);
    for (int l = 0; l <= maxLevel(); ++l) {
      _tmpLevelIndexSets.push_back(new LevelIndexSetImp(*this,_hostGrid.levelView(l)));
    }
    _tmpLeafIndexSet->update(_tmpLevelIndexSets,true);
    _state = preUpdate;
  }

  void updateSubDomains() {
    assert(_state == preUpdate);
    _leafIndexSet.swap(*_tmpLeafIndexSet);
    for (int l = 0; l <= maxLevel(); ++l) {
      _levelIndexSets[l]->swap(*_tmpLevelIndexSets[l]);
    }
    _state = postUpdate;
  }

  void postUpdateSubDomains() {
    assert(_state == postUpdate);
    _tmpLevelIndexSets.clear();
    _tmpLeafIndexSet.reset(NULL);
    _state = fixed;
  }

  void addToSubDomain(SubDomainType subDomain, const typename Traits::template Codim<0>::Entity& e) {
    assert(_state == marking);
    assert(e.isLeaf());
    _tmpLeafIndexSet->addToSubDomain(subDomain,e);
  }


  void removeFromSubDomain(SubDomainType subDomain, const typename Traits::template Codim<0>::Entity& e) {
    assert(_state == marking);
    assert(e.isLeaf());
    _tmpLeafIndexSet->removeFromSubDomain(subDomain,e);
  }


  void assignToSubDomain(SubDomainType subDomain, const typename Traits::template Codim<0>::Entity& e) {
    assert(_state == marking);
    assert(e.isLeaf());
    _tmpLeafIndexSet->assignToSubDomain(subDomain,e);
  }

  const SubDomainGrid& subDomain(SubDomainType subDomain) const {
    _subDomainGrid.reset(subDomain);
    _subDomainGrid.update();
    return _subDomainGrid;
  }

  SubDomainGridPointer subDomainPointer(SubDomainType subDomain) const {
    return SubDomainGridPointer(*this,subDomain);
  }

private:

  HostGrid& _hostGrid;

  std::vector<boost::shared_ptr<LevelIndexSetImp> > _levelIndexSets;
  LeafIndexSetImp _leafIndexSet;

  std::vector<boost::shared_ptr<LevelIndexSetImp> > _tmpLevelIndexSets;
  boost::scoped_ptr<LeafIndexSetImp> _tmpLeafIndexSet;

  GlobalIdSetImp _globalIdSet;
  LocalIdSetImp _localIdSet;

  State _state;

  mutable SubDomainGrid _subDomainGrid;

  template<typename Entity>
  struct HostEntity {
    typedef typename HostGrid::Traits::template Codim<Entity::codimension>::Entity type;
  };

  template<typename Entity>
  struct HostEntityPointer {
    typedef typename HostGrid::Traits::template Codim<Entity::codimension>::EntityPointer type;
  };

  /*
  template<typename EntityType>
  typename HostEntity<EntityType>::type& hostEntity(EntityType& e) const {
    return getRealImplementation(e).wrappedEntity();
    }*/

  template<typename EntityType>
  const typename HostEntity<EntityType>::type& hostEntity(const EntityType& e) const {
    return *(getRealImplementation(e).hostEntityPointer());
  }

  template<typename EntityType>
  const typename HostEntityPointer<EntityType>::type& hostEntityPointer(const EntityType& e) const {
    return getRealImplementation(e).hostEntityPointer();
  }

  void updateIndexSets() {
    // make sure we have enough LevelIndexSets
    while (_levelIndexSets.size() <= maxLevel()) {
      _levelIndexSets.push_back(new LevelIndexSetImp(*this,_hostGrid.levelView(_levelIndexSets.size())));
    }

    _leafIndexSet.reset(true);
    _leafIndexSet.update(_levelIndexSets,true);

    _globalIdSet.update(_hostGrid.globalIdSet());
    _localIdSet.update(_hostGrid.localIdSet());
  }

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_MULTIDOMAINGRID_HH
