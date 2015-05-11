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
#include <dune/grid/multidomaingrid/subdomaingrid/iterator.hh>
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

private:

  static const int dim  = MDGrid::dimension;
  static const int dimw = MDGrid::dimensionworld;

public:

  struct Traits
  {
    /** \brief The type that is implementing the grid. */
    using Grid = SubDomainGrid<MDGrid>;


    using LeafIntersection = Dune::Intersection<
      const Grid,
      IntersectionWrapper<
        const Grid,
        IndexSetWrapper<
          const Grid,
          typename MDGrid::LeafIndexSetImp
          >,
        typename MDGrid::LeafGridView::Intersection
        >
      >;

    using LevelIntersection = Dune::Intersection<
      const Grid,
      IntersectionWrapper<
        const Grid,
        IndexSetWrapper<
          const Grid,
          typename MDGrid::LevelIndexSetImp
          >,
        typename MDGrid::LevelGridView::Intersection
        >
      >;

    using LeafIntersectionIterator = Dune::IntersectionIterator<
      const Grid,
      IntersectionIteratorWrapper<
        const Grid,
        IndexSetWrapper<
          const Grid,
          typename MDGrid::LeafIndexSetImp
          >,
        typename MDGrid::LeafGridView::IntersectionIterator
        >,
      IntersectionWrapper<
        const Grid,
        IndexSetWrapper<
          const Grid,
          typename MDGrid::LeafIndexSetImp
          >,
        typename MDGrid::LeafGridView::Intersection
        >
      >;

    using LevelIntersectionIterator = Dune::IntersectionIterator<
      const Grid,
      IntersectionIteratorWrapper<
        const Grid,
        IndexSetWrapper<
          const Grid,
          typename MDGrid::LevelIndexSetImp
          >,
        typename MDGrid::LevelGridView::IntersectionIterator
        >,
      IntersectionWrapper<
        const Grid,
        IndexSetWrapper<
          const Grid,
          typename MDGrid::LevelIndexSetImp
          >,
        typename MDGrid::LevelGridView::Intersection
        >
      >;

    using HierarchicIterator = Dune::EntityIterator< 0, const Grid, HierarchicIteratorWrapper< const Grid > >;


    template <PartitionIteratorType pitype>
    struct Partition
    {

      using LevelGridView = Dune::GridView<LevelGridViewTraits<const Grid,pitype> >;

      using LeafGridView = Dune::GridView<LeafGridViewTraits<const Grid,pitype> >;

    };

    template <int cd>
    struct Codim
    {

      using Geometry      = Dune::Geometry<dim-cd, dimw, const Grid, GeometryWrapper>;
      using LocalGeometry = Dune::Geometry<dim-cd, dim, const Grid, LocalGeometryWrapper>;

      using Entity        = Dune::Entity<cd, dim, const Grid, EntityWrapper>;
      using EntityPointer = Dune::EntityPointer<const Grid, EntityPointerWrapper<cd,const Grid> >;

      using EntitySeed    = EntitySeedWrapper<typename MDGrid::HostGrid::template Codim<cd>::EntitySeed>;

      template <PartitionIteratorType pitype>
      struct Partition
      {

        using LevelIterator = Dune::EntityIterator<
          cd,
          const Grid,
          IteratorWrapper<
            typename Traits::template Partition<pitype>::LevelGridView,
            typename MDGrid::LevelGridView::template Codim<cd>::template Partition<pitype>::Iterator,
            cd,
            pitype,
            const Grid
            >
          >;

        using LeafIterator = Dune::EntityIterator<
          cd,
          const Grid,
          IteratorWrapper<
            typename Traits::template Partition<pitype>::LeafGridView,
            typename MDGrid::LeafGridView::template Codim<cd>::template Partition<pitype>::Iterator,
            cd,
            pitype,
            const Grid
            >
          >;
      };

      using LeafIterator  = typename Partition< All_Partition >::LeafIterator;
      using LevelIterator = typename Partition< All_Partition >::LevelIterator;

    private:
      friend class Dune::Entity<cd, dim, const Grid, EntityWrapper>;

    };

    using LevelIndexSet = IndexSetWrapper<const Grid, typename MDGrid::LevelIndexSetImp>;
    using LeafIndexSet  = IndexSetWrapper<const Grid, typename MDGrid::LeafIndexSetImp>;

    using GlobalIdSet = IdSet<
      const Grid,
      IdSetWrapper<
        const Grid,
        typename MDGrid::HostGrid::Traits::GlobalIdSet
        >,
      typename MDGrid::HostGrid::Traits::GlobalIdSet::IdType
      >;

    using LocalIdSet = IdSet<
      const Grid,
      IdSetWrapper<
        const Grid,
        typename MDGrid::HostGrid::Traits::LocalIdSet
        >,
      typename MDGrid::HostGrid::Traits::LocalIdSet::IdType
      >;

    using CollectiveCommunication = typename MDGrid::CollectiveCommunication;
  };

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
  friend class EntityWrapperBase;

  template<int codim, int dim, typename GridImp>
  friend class EntityWrapper;

  template<int codim, typename GridImp>
  friend class EntityPointerWrapper;

  template<typename, typename, int codim, PartitionIteratorType pitype, typename GridImp>
  friend class IteratorWrapper;

  template<typename GridImp>
  friend class HierarchicIteratorWrapper;

  template<int mydim, int coorddim, typename GridImp>
  friend class GeometryWrapper;

  template<int mydim, int coorddim, typename GridImp>
  friend class LocalGeometryWrapper;

  template<typename GridImp, typename WrappedIndexSet>
  friend class IndexSetWrapper;

  template<typename GridImp, typename WrappedIdSet>
  friend class IdSetWrapper;

  template<typename GridImp>
  friend struct ::Dune::mdgrid::detail::HostGridAccessor;

  template<typename,typename,typename>
  friend class IntersectionIteratorWrapper;

  template<typename,typename,typename>
  friend class IntersectionWrapper;

  template<typename, PartitionIteratorType>
  friend class LevelGridView;

  template<typename, PartitionIteratorType>
  friend class LeafGridView;

  typedef GridDefaultImplementation<MDGrid::dimension,
                                    MDGrid::dimensionworld,
                                    typename MDGrid::ctype,
                                    SubDomainGridFamily<MDGrid>
                                    > BaseT;

  using GridImp = SubDomainGrid<MDGrid>;

public:

  using MultiDomainGrid = MDGrid;
  using MDGridType DUNE_DEPRECATED_MSG("Deprecated in 2.4, Use MultiDomainGrid instead") = MultiDomainGrid;

  using HostGrid = typename MDGrid::HostGrid;
  using HostGridType DUNE_DEPRECATED_MSG("Deprecated in 2.4, Use HostGrid instead") = HostGrid;

private:

  typedef IndexSetWrapper<const SubDomainGrid<MDGrid>, typename MDGrid::LevelIndexSetImp> LevelIndexSetImp;

  typedef IndexSetWrapper<const SubDomainGrid<MDGrid>, typename MDGrid::LeafIndexSetImp> LeafIndexSetImp;

  typedef IdSetWrapper<const SubDomainGrid<MDGrid>, typename HostGrid::Traits::GlobalIdSet> GlobalIdSetImp;

  typedef IdSetWrapper<const SubDomainGrid<MDGrid>, typename HostGrid::Traits::LocalIdSet> LocalIdSetImp;

public:

  typedef SubDomainGridFamily<MDGrid> GridFamily;
  typedef typename GridFamily::Traits Traits;

  /** \brief The type used for coordinates */
  typedef typename MDGrid::ctype ctype;

  /** \brief The type used for subdomain numbers */
  typedef typename MDGrid::SubDomainIndex SubDomainIndex;

  enum IntersectionType { neighbor, foreign, boundary, processor };

  using BaseT::dimension;
  using BaseT::dimensionworld;

  /** @name Dune grid interface methods */
  /*@{*/

  //! Reconstruct EntityPointer from EntitySeed
  template<typename EntitySeed>
  typename Traits::template Codim<EntitySeed::codimension>::EntityPointer
  DUNE_DEPRECATED_MSG("entityPointer() is deprecated and will be removed after the release of dune-grid 2.4. Use entity() instead to directly obtain an Entity object.")
  entityPointer(const EntitySeed& entitySeed) const
  {
    return
      EntityPointerWrapper<EntitySeed::codimension,const GridImp>(
        this,
        typename MDGrid::template Codim<EntitySeed::codimension>::EntityPointer(
          _grid.entityPointer(entitySeed)
        )
      );
  }

  template<typename EntitySeed>
  typename Traits::template Codim<EntitySeed::codimension>::Entity
  entity(const EntitySeed& entitySeed) const
  {
    return
      EntityWrapper<EntitySeed::codimension,dimension,const GridImp>(
        this,
        typename MDGrid::template Codim<EntitySeed::codimension>::Entity(
          _grid.entity(entitySeed)
        )
      );
  }

  int maxLevel() const {
    return _grid.maxLevel();
  }

  template<int codim>
  typename Traits::template Codim<codim>::LevelIterator lbegin(int level) const {
    return IteratorWrapper<
      typename Traits::template Partition<All_Partition>::LevelGridView,
      typename MultiDomainGrid::LevelGridView::template Codim<codim>::template Partition<All_Partition>::Iterator,
      codim,
      All_Partition,
      const GridImp>(
                     this,
                     &this->levelGridView(level).indexSet(),
                     _grid.levelGridView(level).template begin<codim>(),
                     _grid.levelGridView(level).template end<codim>()
                     );
  }

  template<int codim>
  typename Traits::template Codim<codim>::LevelIterator lend(int level) const {
    return IteratorWrapper<
      typename Traits::template Partition<All_Partition>::LevelGridView,
      typename MultiDomainGrid::LevelGridView::template Codim<codim>::template Partition<All_Partition>::Iterator,
      codim,
      All_Partition,
      const GridImp>(
                     this,
                     &this->levelGridView(level).indexSet(),
                     _grid.levelGridView(level).template end<codim>(),
                     _grid.levelGridView(level).template end<codim>()
                     );
  }

  template<int codim, PartitionIteratorType pitype>
  typename Traits::template Codim<codim>::template Partition<pitype>::LevelIterator lbegin(int level) const {
    return IteratorWrapper<
      typename Traits::template Partition<pitype>::LevelGridView,
      typename MultiDomainGrid::LevelGridView::template Codim<codim>::template Partition<pitype>::Iterator,
      codim,
      pitype,
      const GridImp>(
                     this,
                     &this->levelGridView(level).indexSet(),
                     _grid.levelGridView(level).template begin<codim,pitype>(),
                     _grid.levelGridView(level).template end<codim,pitype>()
                     );
  }

  template<int codim, PartitionIteratorType pitype>
  typename Traits::template Codim<codim>::template Partition<pitype>::LevelIterator lend(int level) const {
    return IteratorWrapper<
      typename Traits::template Partition<pitype>::LevelGridView,
      typename MultiDomainGrid::LevelGridView::template Codim<codim>::template Partition<pitype>::Iterator,
      codim,
      pitype,
      const GridImp>(
                     this,
                     &this->levelGridView(level).indexSet(),
                     _grid.levelGridView(level).template end<codim,pitype>(),
                     _grid.levelGridView(level).template end<codim,pitype>()
                     );
  }

  template<int codim>
  typename Traits::template Codim<codim>::LeafIterator leafbegin() const {
    return IteratorWrapper<
      typename Traits::template Partition<All_Partition>::LeafGridView,
      typename MultiDomainGrid::LeafGridView::template Codim<codim>::template Partition<All_Partition>::Iterator,
      codim,
      All_Partition,
      const GridImp>(
                     this,
                     &this->leafGridView().indexSet(),
                     _grid.leafGridView().template begin<codim>(),
                     _grid.leafGridView().template end<codim>()
                     );
  }

  template<int codim>
  typename Traits::template Codim<codim>::LeafIterator leafend() const {
    return IteratorWrapper<
      typename Traits::template Partition<All_Partition>::LeafGridView,
      typename MultiDomainGrid::LeafGridView::template Codim<codim>::template Partition<All_Partition>::Iterator,
      codim,
      All_Partition,
      const GridImp>(
                     this,
                     &this->leafGridView().indexSet(),
                     _grid.leafGridView().template end<codim>(),
                     _grid.leafGridView().template end<codim>()
                     );
  }

  template<int codim, PartitionIteratorType pitype>
  typename Traits::template Codim<codim>::template Partition<pitype>::LeafIterator leafbegin() const {
    return IteratorWrapper<
      typename Traits::template Partition<pitype>::LeafGridView,
      typename MultiDomainGrid::LeafGridView::template Codim<codim>::template Partition<pitype>::Iterator,
      codim,
      pitype,
      const GridImp>(
                     this,
                     &this->leafGridView().indexSet(),
                     _grid.leafGridView().template begin<codim,pitype>(),
                     _grid.leafGridView().template end<codim,pitype>()
                     );
  }

  template<int codim, PartitionIteratorType pitype>
  typename Traits::template Codim<codim>::template Partition<pitype>::LeafIterator leafend() const {
    return IteratorWrapper<
      typename Traits::template Partition<pitype>::LeafGridView,
      typename MultiDomainGrid::LeafGridView::template Codim<codim>::template Partition<pitype>::Iterator,
      codim,
      pitype,
      const GridImp>(
                     this,
                     &this->leafGridView().indexSet(),
                     _grid.leafGridView().template end<codim,pitype>(),
                     _grid.leafGridView().template end<codim,pitype>()
                     );
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
        _levelIndexSets.push_back(std::make_shared<LevelIndexSetImp>(*this,_grid.levelIndexSet(_levelIndexSets.size())));
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
  typename Traits::template Codim<EntityType::codimension>::Entity subDomainEntity(const EntityType& mdEntity) const {
    return EntityWrapper<EntityType::codimension,dimension,const GridImp>(this,mdEntity);
  }

  template<typename EntityType>
  static const typename MDGrid::template MultiDomainEntity<EntityType>::type& multiDomainEntity(const EntityType& e) {
    return SubDomainGrid::getRealImplementation(e).multiDomainEntity();
  }

  template<typename EntityType>
  static typename MDGrid::template MultiDomainEntityPointer<EntityType>::type multiDomainEntityPointer(const EntityType& e) {
    return SubDomainGrid::getRealImplementation(e).multiDomainEntity();
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
      return IntersectionIteratorWrapper<
        const GridImp,
        typename GridImp::LeafGridView::IndexSet,
        typename MDGrid::LeafGridView::IntersectionIterator
        >(
          &this->leafGridView().indexSet(),
          it->firstMultiDomainIntersectionIterator()
          );
    else
      return IntersectionIteratorWrapper<
        const GridImp,
        typename GridImp::LeafGridView::IndexSet,
        typename MDGrid::LeafGridView::IntersectionIterator
        >(
          &this->leafGridView().indexSet(),
          it->secondMultiDomainIntersectionIterator()
          );
  }

  typename Traits::LevelIntersectionIterator subDomainIntersectionIterator(const typename MDGrid::LevelSubDomainInterfaceIterator it) const {
    assert(_subDomain == it->subDomain1() || _subDomain == it->subDomain2());
    if (_subDomain == it->subDomain1())
      return IntersectionIteratorWrapper<
        const GridImp,
        typename GridImp::LevelGridView::IndexSet,
        typename MDGrid::LevelGridView::IntersectionIterator
        >(
          &this->levelGridView(it->firstMultiDomainIntersectionIterator()->inside().level()).indexSet(),
          it->firstMultiDomainIntersectionIterator()
          );
    else
      return IntersectionIteratorWrapper<
        const GridImp,
        typename GridImp::LevelGridView::IndexSet,
        typename MDGrid::LevelGridView::IntersectionIterator
        >(
          &this->levelGridView(it->secondMultiDomainIntersectionIterator()->inside().level()).indexSet(),
          it->secondMultiDomainIntersectionIterator()
          );
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
        return _impl.size(_grid.subDomainEntity(_grid._grid.wrapHostEntity(e)));
      else
        return 0;
    }

    template<typename MessageBufferImp, typename Entity>
    void gather(MessageBufferImp& buf, const Entity& e) const
    {
      if (_grid.containsHostEntity(e))
        _impl.gather(buf,_grid.subDomainEntity(_grid._grid.wrapHostEntity(e)));
    }

    template<typename MessageBufferImp, typename Entity>
    void scatter(MessageBufferImp& buf, const Entity& e, std::size_t n)
    {
      if (_grid.containsHostEntity(e))
        _impl.scatter(buf,_grid.subDomainEntity(_grid._grid.wrapHostEntity(e)),n);
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
