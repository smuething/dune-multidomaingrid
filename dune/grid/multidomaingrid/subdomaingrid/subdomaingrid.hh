#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_HH

#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/grid/multidomaingrid/hostgridaccessor.hh>
#include <dune/grid/multidomaingrid/subdomainset.hh>

#include <dune/grid/multidomaingrid/subdomaingrid/geometry.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/localgeometry.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/entity.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/entitypointer.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/leafiterator.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/leveliterator.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/hierarchiciterator.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/intersection.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/intersectioniterator.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/idsets.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/indexsets.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/gridview.hh>


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
            template<int,int,class> class LocalGeometryImp,
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
            template<class,PartitionIteratorType> class LevelGridViewTraits,
            template<class,PartitionIteratorType> class LeafGridViewTraits
            >
  struct SubDomainGridTraits
  {
    /** \brief The type that is implementing the grid. */
    typedef GridImp Grid;

    /** \brief The type of the intersection at the leafs of the grid. */
    typedef Dune::Intersection<const GridImp, LeafIntersectionImp<const GridImp> >  LeafIntersection;
    /** \brief The type of the intersection at the levels of the grid. */
    typedef Dune::Intersection<const GridImp, LevelIntersectionImp<const GridImp> > LevelIntersection;
    /** \brief The type of the intersection iterator at the leafs of the grid. */
    typedef Dune::IntersectionIterator<const GridImp, LeafIntersectionIteratorImp<const GridImp>, LeafIntersectionImp<const GridImp> >   LeafIntersectionIterator;
    /** \brief The type of the intersection iterator at the levels of the grid. */
    typedef Dune::IntersectionIterator<const GridImp, LevelIntersectionIteratorImp<const GridImp>, LevelIntersectionImp<const GridImp> > LevelIntersectionIterator;

    /** \brief The type of the  hierarchic iterator. */
    typedef Dune::EntityIterator< 0, const GridImp, HierarchicIteratorImp< const GridImp > > HierarchicIterator;

    /**
     * \brief Traits associated with a specific codimension.
     * \tparam cd The codimension.
     */
    template <int cd>
    struct Codim
    {
      //! IMPORTANT: Codim<codim>::Geometry == Geometry<dim-codim,dimw>
      /** \brief The type of the geometry associated with the entity.*/
      typedef Dune::Geometry<dim-cd, dimw, const GridImp, GeometryImp> Geometry;
      /** \brief The type of the local geometry associated with the entity.*/
      typedef Dune::Geometry<dim-cd, dim, const GridImp, LocalGeometryImp> LocalGeometry;
      /** \brief The type of the entity. */
      // we could - if needed - introduce another struct for dimglobal of Geometry
      typedef Dune::Entity<cd, dim, const GridImp, EntityImp> Entity;

      /** \brief The type of the entity pointer for entities of this codim.*/
      typedef Dune::EntityPointer<const GridImp,EntityPointerImp<cd,const GridImp> > EntityPointer;

      typedef EntitySeedWrapper<typename MDGrid::HostGridType::template Codim<cd>::EntitySeed> EntitySeed;

      /**
       * \brief Traits associated with a specific grid partition type.
       * \tparam pitype The type of the grid partition.
       */
      template <PartitionIteratorType pitype>
      struct Partition
      {
        /** \brief The type of the iterator over the level entities of this codim on this partition. */
        typedef Dune::EntityIterator< cd, const GridImp, LevelIteratorImp< cd, pitype, const GridImp > > LevelIterator;
        /** \brief The type of the iterator over the leaf entities of this codim on this partition. */
        typedef Dune::EntityIterator< cd, const GridImp, LeafIteratorImp< cd, pitype, const GridImp > > LeafIterator;
      };

      /** \brief The type of the iterator over all leaf entities of this codim. */
      typedef typename Partition< All_Partition >::LeafIterator LeafIterator;

      /** \brief The type of the entity pointer for entities of this codim.*/
      typedef typename Partition< All_Partition >::LevelIterator LevelIterator;

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
    LocalGeometryWrapper,
    EntityWrapper,
    EntityPointerWrapper,
    LevelIteratorWrapper,
    LeafIntersectionWrapper, // leaf intersection
    LevelIntersectionWrapper, // level intersection
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
    typename MDGrid::HostGridType::CollectiveCommunication,
    LevelGridViewTraits,
    LeafGridViewTraits
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

  template<int mydim, int coorddim, typename GridImp>
  friend class LocalGeometryWrapper;

  template<int mydim, int coorddim, typename GridImp>
  friend class MakeableLocalGeometryWrapper;

  template<typename GridImp, typename WrappedIndexSet>
  friend class IndexSetWrapper;

  template<typename GridImp, typename WrappedIdSet>
  friend class IdSetWrapper;

  template<typename GridImp>
  friend struct ::Dune::mdgrid::detail::HostGridAccessor;

  template<typename,typename,typename,typename,typename>
  friend class IntersectionIteratorWrapper;

  template<typename GridImp>
  friend class LeafIntersectionIteratorWrapper;

  template<typename GridImp>
  friend class LeafIntersectionWrapper;

  template<typename GridImp>
  friend class LevelIntersectionIteratorWrapper;

  template<typename GridImp>
  friend class LevelIntersectionWrapper;

  template<typename, PartitionIteratorType>
  friend class LevelGridView;

  template<typename, PartitionIteratorType>
  friend class LeafGridView;

  typedef GridDefaultImplementation<MDGrid::dimension,
                                    MDGrid::dimensionworld,
                                    typename MDGrid::ctype,
                                    SubDomainGridFamily<MDGrid>
                                    > BaseT;

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

  /** \brief The type used for coordinates */
  typedef typename MDGrid::ctype ctype;

  /** \brief The type used for subdomain numbers */
  typedef typename MDGrid::SubDomainIndex SubDomainIndex;
  typedef SubDomainIndex SubDomainIndexType DUNE_DEPRECATED_MSG("Use SubDomainIndex instead.");
  typedef SubDomainIndex SubDomainType DUNE_DEPRECATED_MSG("Use SubDomainIndex instead.");

  typedef MDGridType MultiDomainGrid;

  enum IntersectionType { neighbor, foreign, boundary, processor };

public:

  /** @name Dune grid interface methods */
  /*@{*/

  //! Reconstruct EntityPointer from EntitySeed
  template<typename EntitySeed>
  typename Traits::template Codim<EntitySeed::codimension>::EntityPointer
  entityPointer(const EntitySeed& entitySeed) const
  {
    return
      EntityPointerWrapper<EntitySeed::codimension,const GridImp>(
        *this,
        typename MDGrid::template Codim<EntitySeed::codimension>::EntityPointer(
          _grid.entityPointer(entitySeed)
        )
      );
  }

  int maxLevel() const {
    return _grid.maxLevel();
  }

  template<int codim>
  typename Traits::template Codim<codim>::LevelIterator lbegin(int level) const {
    return LevelIteratorWrapper<codim,All_Partition,const GridImp>(levelIndexSet(level),
                                                                   _grid.template lbegin<codim>(level),
                                                                   _grid.template lend<codim>(level));
  }

  template<int codim>
  typename Traits::template Codim<codim>::LevelIterator lend(int level) const {
    return LevelIteratorWrapper<codim,All_Partition,const GridImp>(levelIndexSet(level),
                                                                   _grid.template lend<codim>(level),
                                                                   _grid.template lend<codim>(level));
  }

  template<int codim, PartitionIteratorType PiType>
  typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin(int level) const {
    return LevelIteratorWrapper<codim,PiType,const GridImp>(levelIndexSet(level),
                                                            _grid.template lbegin<codim,PiType>(level),
                                                            _grid.template lend<codim,PiType>(level));
  }

  template<int codim, PartitionIteratorType PiType>
  typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend(int level) const {
    return LevelIteratorWrapper<codim,PiType,const GridImp>(levelIndexSet(level),
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

  //! Use MultiDomainGrid::globalRefine() instead of this method.
  /**
   * Like all grid modification methods, globalRefine() may ONLY be called on the underlying
   * MultiDomainGrid.
   *
   * \throws NotImplemented calling globalRefine() will always throw this exception.
   */
  void globalRefine(int refCount) {
    DUNE_THROW(NotImplemented,"grid modification only allowed on the MultiDomainGrid");
  }

  //! Mark the entity e for grid refinement across ALL subdomains.
  /**
   * This method marks the passed entity for refinement on the underlying MultiDomainGrid.
   * When using this method, keep in mind that
   * a) this will have an effect across all subdomains and also on the MultiDomainGrid.
   * b) the exact semantics of refinement depend on the host grid.
   *
   * \param refCount the refinement mark to set, for exact semantics see the documentation
   *                 of the host grid.
   * \param e        the entity to refine / coarsen.
   */
  bool mark(int refCount, const typename Traits::template Codim<0>::Entity& e) {
    return _grid.mark(refCount,multiDomainEntity(e));
  }

  //! Retrieve the refinement mark of entity e.
  /**
   * This method returns the refinement mark set on the entity by a call to mark(). As the
   * mark might have been set on any subdomain containing e, the caller should not assume that
   * entities always carry the mark assigned within the current subdomain.
.
   * \param e        the entity for which to return the refinement mark.
   */
  int getMark(const typename Traits::template Codim<0>::Entity& e) {
    return _grid.getMark(multiDomainEntity(e));
  }

  //! Use MultiDomainGrid::preAdapt() instead of this method.
  /**
   * Like all grid modification methods, preAdapt() may ONLY be called on the underlying
   * MultiDomainGrid.
   *
   * \throws NotImplemented calling preAdapt() will always throw this exception.
   */
  bool preAdapt() {
    DUNE_THROW(NotImplemented,"grid modification only allowed on the MultiDomainGrid");
  }

  //! Use MultiDomainGrid::adapt() instead of this method.
  /**
   * Like all grid modification methods, adapt() may ONLY be called on the underlying
   * MultiDomainGrid.
   *
   * \throws NotImplemented calling adapt() will always throw this exception.
   */
  bool adapt() {
    DUNE_THROW(NotImplemented,"grid modification only allowed on the MultiDomainGrid");
  }

  //! Use MultiDomainGrid::postAdapt() instead of this method.
  /**
   * Like all grid modification methods, postAdapt() may ONLY be called on the underlying
   * MultiDomainGrid.
   *
   * \throws NotImplemented calling postAdapt() will always throw this exception.
   */
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

  template<typename DataHandleImp, typename DataTypeImp>
  void communicate (CommDataHandleIF<DataHandleImp,DataTypeImp> &data,
                    InterfaceType iftype,
                    CommunicationDirection dir,
                    int level) const
  {
    DataHandleWrapper<CommDataHandleIF<DataHandleImp,DataTypeImp> > datahandle(data,*this);
    _grid._hostGrid.communicate(datahandle,iftype,dir,level);
  }

  template<typename DataHandleImp, typename DataTypeImp>
  void communicate (CommDataHandleIF<DataHandleImp,DataTypeImp> &data,
                    InterfaceType iftype,
                    CommunicationDirection dir) const
  {
    DataHandleWrapper<CommDataHandleIF<DataHandleImp,DataTypeImp> > datahandle(data,*this);
    _grid._hostGrid.communicate(datahandle,iftype,dir);
  }

  size_t numBoundarySegments() const
  {
    return _grid.numBoundarySegments();
  }

  /*@}*/

  /** \brief Get the MultiDomainGrid that we are part of */
  const MDGrid& multiDomainGrid() const {
    return _grid;
  }

  /** \brief Return our subdomain number */
  SubDomainIndex domain() const {
    return _subDomain;
  }

  void update() const {
    if (_grid.supportLevelIndexSets()) {
      while (_levelIndexSets.size() <= static_cast<std::size_t>(maxLevel())) {
        _levelIndexSets.push_back(make_shared_ptr(new LevelIndexSetImp(*this,_grid.levelIndexSet(_levelIndexSets.size()))));
      }
    }
  }

  bool operator==(const SubDomainGrid& rhs) const {
    return (&_grid == &rhs._grid && _subDomain == rhs._subDomain);
  }

  /** @name Entity conversion methods */
  /*@{*/
  template<int cc>
  typename Traits::template Codim<cc>::EntityPointer subDomainEntityPointer(const typename MDGrid::Traits::template Codim<cc>::Entity& mdEntity) const {
    return EntityPointerWrapper<cc,const SubDomainGrid<MDGrid> >(*this,typename MDGrid::Traits::template Codim<cc>::EntityPointer(mdEntity));
  }

  template<typename EntityType>
  typename Traits::template Codim<EntityType::codimension>::EntityPointer subDomainEntityPointer(const EntityType& mdEntity) const {
    return subDomainEntityPointer<EntityType::codimension>(mdEntity);
  }

  template<typename EntityType>
  static const typename MDGrid::template MultiDomainEntity<EntityType>::type& multiDomainEntity(const EntityType& e) {
    return *(SubDomainGrid::getRealImplementation(e).multiDomainEntityPointer());
  }

  template<typename EntityType>
  static typename MDGrid::template MultiDomainEntityPointer<EntityType>::type multiDomainEntityPointer(const EntityType& e) {
    return SubDomainGrid::getRealImplementation(e).multiDomainEntityPointer();
  }

  template<typename EntityType>
  static const typename MDGrid::template HostEntity<EntityType>::type& hostEntity(const EntityType& e) {
    return SubDomainGrid::getRealImplementation(e).hostEntity();
  }

  template<typename EntityType>
  static typename MDGrid::template HostEntityPointer<EntityType>::type hostEntityPointer(const EntityType& e) {
    return typename MDGrid::template HostEntityPointer<EntityType>::type(SubDomainGrid::getRealImplementation(e).hostEntity());
  }

  template<typename IntersectionType>
  struct MultiDomainIntersection
  {
    typedef typename BaseT::template ReturnImplementationType<IntersectionType>::ImplementationType::MultiDomainIntersection Type;
  };

  template<typename IntersectionType>
  static const typename MultiDomainIntersection<IntersectionType>::Type& multiDomainIntersection(const IntersectionType& is) {
    return SubDomainGrid::getRealImplementation(is).multiDomainIntersection();
  }

  /*@}*/

  typename Traits::LeafIntersectionIterator subDomainIntersectionIterator(const typename MDGrid::LeafSubDomainInterfaceIterator it) const {
    assert(_subDomain == it->subDomain1() || _subDomain == it->subDomain2());
    if (_subDomain == it->subDomain1())
      return LeafIntersectionIteratorWrapper<const GridImp>(*this,it->firstMultiDomainIntersectionIterator());
    else
      return LeafIntersectionIteratorWrapper<const GridImp>(*this,it->secondMultiDomainIntersectionIterator());
  }

  typename Traits::LevelIntersectionIterator subDomainIntersectionIterator(const typename MDGrid::LevelSubDomainInterfaceIterator it) const {
    assert(_subDomain == it->subDomain1() || _subDomain == it->subDomain2());
    if (_subDomain == it.domain1())
      return LevelIntersectionIteratorWrapper<const GridImp>(*this,it->firstCell()->level(),it->firstMultiDomainIntersectionIterator());
    else
      return LevelIntersectionIteratorWrapper<const GridImp>(*this,it->secondCell()->level(),it->secondMultiDomainIntersectionIterator());
  }

  template<typename Intersection>
  IntersectionType intersectionType(const Intersection& intersection) const {
    return SubDomainGrid::getRealImplementation(intersection).intersectionType();
  }

private:

  MDGrid& _grid;
  SubDomainIndex _subDomain;
  GlobalIdSetImp _globalIdSet;
  LocalIdSetImp _localIdSet;
  LeafIndexSetImp _leafIndexSet;
  mutable std::vector<Dune::shared_ptr<LevelIndexSetImp> > _levelIndexSets;

  SubDomainGrid(MDGrid& grid, SubDomainIndex subDomain) :
    _grid(grid),
    _subDomain(subDomain),
    _globalIdSet(*this,grid._hostGrid.globalIdSet()),
    _localIdSet(*this,grid._hostGrid.localIdSet()),
    _leafIndexSet(*this,grid.leafIndexSet())
  {
    update();
  }

  template<typename EntityType>
  bool containsMultiDomainEntity(const EntityType& e) const {
    if (_grid.supportLevelIndexSets())
      return levelIndexSet(e.level()).containsMultiDomainEntity(e);
    else
      return leafIndexSet().containsMultiDomainEntity(e);
  }

  template<typename EntityType>
  bool containsHostEntity(const EntityType& e) const {
    if (_grid.supportLevelIndexSets())
      return levelIndexSet(e.level()).containsHostEntity(e);
    else
      return leafIndexSet().containsHostEntity(e);
  }

  SubDomainGrid(const SubDomainGrid& rv);
  SubDomainGrid& operator=(const SubDomainGrid& rv);


  template<typename Impl>
  struct DataHandleWrapper
    : public Dune::CommDataHandleIF<DataHandleWrapper<Impl>,
                                    typename Impl::DataType
                                    >
  {

    bool contains(int dim, int codim) const
    {
      return _impl.contains(dim,codim); // TODO: check if codim supported
    }

    bool fixedsize(int dim, int codim) const
    {
      //_impl.fixedsize(dim,codim); // TODO: warning if true?
      return false;
    }

    template<typename Entity>
    std::size_t size(const Entity& e) const
    {
      if (_grid.containsHostEntity(e))
        return _impl.size(*_grid.subDomainEntityPointer(*_grid._grid.wrapHostEntity(e)));
      else
        return 0;
    }

    template<typename MessageBufferImp, typename Entity>
    void gather(MessageBufferImp& buf, const Entity& e) const
    {
      if (_grid.containsHostEntity(e))
        _impl.gather(buf,*_grid.subDomainEntityPointer(*_grid._grid.wrapHostEntity(e)));
    }

    template<typename MessageBufferImp, typename Entity>
    void scatter(MessageBufferImp& buf, const Entity& e, std::size_t n)
    {
      if (_grid.containsHostEntity(e))
        _impl.scatter(buf,*_grid.subDomainEntityPointer(*_grid._grid.wrapHostEntity(e)),n);
    }

    DataHandleWrapper(Impl& impl, const SubDomainGrid<MDGrid>& grid)
      : _impl(impl)
      , _grid(grid)
    {}

    Impl& _impl;
    const SubDomainGrid<MDGrid>& _grid;

  };

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_MULTIDOMAINGRID_HH
