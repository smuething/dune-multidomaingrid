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

#include <dune/grid/common/grid.hh>

#include <dune/grid/multidomaingrid/subdomainset.hh>

#include <dune/grid/multidomaingrid/subdomaingrid/subdomaingrid.hh>

#include <dune/grid/multidomaingrid/geometry.hh>
#include <dune/grid/multidomaingrid/entity.hh>
#include <dune/grid/multidomaingrid/entitypointer.hh>
#include <dune/grid/multidomaingrid/leafiterator.hh>
#include <dune/grid/multidomaingrid/leveliterator.hh>
#include <dune/grid/multidomaingrid/hierarchiciterator.hh>
#include <dune/grid/multidomaingrid/intersection.hh>
#include <dune/grid/multidomaingrid/intersectioniterator.hh>
#include <dune/grid/multidomaingrid/idsets.hh>
#include <dune/grid/multidomaingrid/indexsets.hh>
#include <dune/grid/multidomaingrid/mdgridtraits.hh>

#include <dune/grid/multidomaingrid/subdomaintosubdomaininterfaceiterator.hh>
#include <dune/grid/multidomaingrid/allsubdomaininterfacesiterator.hh>

namespace Dune {

namespace mdgrid {

template<typename T>
boost::shared_ptr<T> make_shared_ptr(T* ptr) {
  return boost::shared_ptr<T>(ptr);
}

template<typename HostGrid, typename MDGridTraits>
class MultiDomainGrid;

template<typename HostGrid, typename MDGridTraits>
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
            template<class> class LeafSubDomainInterfaceIteratorImp,
            template<class> class LevelSubDomainInterfaceIteratorImp,
            template<class> class LeafAllSubDomainInterfacesIteratorImp,
            template<class> class LevelAllSubDomainInterfacesIteratorImp,
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

    typedef LeafSubDomainInterfaceIteratorImp<const GridImp> LeafSubDomainInterfaceIterator;
    typedef LevelSubDomainInterfaceIteratorImp<const GridImp> LevelSubDomainInterfaceIterator;

    typedef LeafAllSubDomainInterfacesIteratorImp<const GridImp> LeafAllSubDomainInterfacesIterator;
    typedef LevelAllSubDomainInterfacesIteratorImp<const GridImp> LevelAllSubDomainInterfacesIterator;
  };

  typedef MultiDomainGridTraits<
    HostGrid::dimension,
    HostGrid::dimensionworld,
    MultiDomainGrid<HostGrid,MDGridTraits>,
    GeometryWrapper,
    EntityWrapper,
    EntityPointerWrapper,
    LevelIteratorWrapper,
    LeafIntersectionWrapper, // leaf intersection
    LevelIntersectionWrapper, // level intersection
    LeafIntersectionIteratorWrapper, // leaf intersection iterator
    LevelIntersectionIteratorWrapper, // level intersection iterator
    HierarchicIteratorWrapper,
    LeafIteratorWrapper,
    IndexSetWrapper<const MultiDomainGrid<HostGrid,MDGridTraits>, typename HostGrid::LevelGridView>,
    IndexSetWrapper<const MultiDomainGrid<HostGrid,MDGridTraits>, typename HostGrid::LeafGridView>,
    IdSetWrapper<const MultiDomainGrid<HostGrid,MDGridTraits>, typename HostGrid::Traits::GlobalIdSet>,
    typename HostGrid::Traits::GlobalIdSet::IdType,
    IdSetWrapper<const MultiDomainGrid<HostGrid,MDGridTraits>, typename HostGrid::Traits::LocalIdSet>,
    typename HostGrid::Traits::LocalIdSet::IdType,
    typename HostGrid::CollectiveCommunication,
    LeafSubDomainInterfaceIterator,
    LevelSubDomainInterfaceIterator,
    LeafAllSubDomainInterfacesIterator,
    LevelAllSubDomainInterfacesIterator
    > Traits;

};


//! A meta grid for dividing an existing DUNE grid into subdomains that can be accessed as a grid in their own right.
/**
 * \tparam HostGrid             The type of the underlying grid implementation.
 * \tparam MDGridTraitsType     A traits type for customizing how the MultiDomainGrid manages the partitioning information.
 */

template<typename HostGrid, typename MDGridTraitsType>
class MultiDomainGrid :
    public GridDefaultImplementation<HostGrid::dimension,
				     HostGrid::dimensionworld,
				     typename HostGrid::ctype,
				     MultiDomainGridFamily<HostGrid,MDGridTraitsType> > {


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

  template<typename,typename,typename,typename>
  friend class IntersectionIteratorWrapper;

  template<typename GridImp>
  friend class LeafIntersectionIteratorWrapper;

  template<typename GridImp>
  friend class LeafIntersectionWrapper;

  template<typename GridImp>
  friend class LevelIntersectionIteratorWrapper;

  template<typename GridImp>
  friend class LevelIntersectionWrapper;

  template<typename>
  friend class subdomain::SubDomainGrid;

  template<typename>
  friend struct subdomain::SubDomainGridFamily;

  template<int,int,typename>
  friend class subdomain::EntityWrapper;

  template <typename>
  friend class LeafSubDomainInterfaceIterator;

  template <typename>
  friend class LevelSubDomainInterfaceIterator;

  template <typename>
  friend class LeafAllSubDomainInterfacesIterator;

  template <typename>
  friend class LevelAllSubDomainInterfacesIterator;

  template<typename,typename,typename,typename>
  friend class SubDomainInterface;

  template<typename,typename,typename,typename,typename>
  friend class subdomain::IntersectionIteratorWrapper;

  template<typename,typename,typename,typename>
  friend class subdomain::IntersectionWrapper;

  typedef MultiDomainGrid<HostGrid,MDGridTraitsType> GridImp;
  typedef HostGrid HostGridType;

  typedef IndexSetWrapper<const GridImp, typename HostGrid::LevelGridView> LevelIndexSetImp;

  typedef IndexSetWrapper<const GridImp, typename HostGrid::LeafGridView> LeafIndexSetImp;

  typedef IdSetWrapper<const GridImp, typename HostGrid::Traits::GlobalIdSet> GlobalIdSetImp;

  typedef IdSetWrapper<const GridImp, typename HostGrid::Traits::LocalIdSet> LocalIdSetImp;

  enum State { stateFixed, stateMarking, statePreUpdate, statePostUpdate, statePreAdapt, statePostAdapt };

  typedef GridImp ThisType;

public:

  typedef MultiDomainGridFamily<HostGrid,MDGridTraitsType> GridFamily;
  typedef typename GridFamily::Traits Traits;
  typedef MDGridTraitsType MDGridTraits;
  typedef typename HostGrid::ctype ctype;

private:

  typedef std::map<typename Traits::LocalIdSet::IdType,typename MDGridTraits::template Codim<0>::SubDomainSet> AdaptationStateMap;

  // typedefs for extracting the host entity types from our own entities

  template<typename Entity>
  struct HostEntity {
    typedef typename HostGrid::Traits::template Codim<Entity::codimension>::Entity type;
  };

  template<typename Entity>
  struct HostEntityPointer {
    typedef typename HostGrid::Traits::template Codim<Entity::codimension>::EntityPointer type;
  };

  // typedefs for extracting the multidomain entity types from subdomain entities

  template<typename Entity>
  struct MultiDomainEntity {
    typedef typename Traits::template Codim<Entity::codimension>::Entity type;
  };

  template<typename Entity>
  struct MultiDomainEntityPointer {
    typedef typename Traits::template Codim<Entity::codimension>::EntityPointer type;
  };


public:

  //! The (integer) type used to identify subdomains.
  typedef typename MDGridTraits::SubDomainIndexType SubDomainIndexType;
  typedef SubDomainIndexType SubDomainType DUNE_DEPRECATED;

  //! The largest number of subdomains any given grid cell may belong to.
  static const std::size_t maxNumberOfSubDomains = MDGridTraits::maxSubDomainsPerCell;
  //! The largest allowed index for a subdomain.
  /**
   * \note As subdomain indices always start at 0, this constant also determines the maximum
   * number of possible subdomains.
   */
  static const SubDomainIndexType maxSubDomainIndex = MDGridTraits::maxSubDomainIndex;

  //! The type used for representing the grid of a subdomain, always a specialization of Dune::mdgrid::subdomain::SubDomainGrid.
  typedef subdomain::SubDomainGrid<ThisType> SubDomainGrid;


  //! The type of the iterators over the codim 1 interface between two subdomains on the leaf view.
  typedef typename Traits::LeafSubDomainInterfaceIterator LeafSubDomainInterfaceIterator;

  //! The type of the iterators over the codim 1 interface between two subdomains on a level view.
  typedef typename Traits::LevelSubDomainInterfaceIterator LevelSubDomainInterfaceIterator;

  //! The type of the iterators over the codim 1 interfaces between all subdomains on the leaf view.
  typedef typename Traits::LeafAllSubDomainInterfacesIterator LeafAllSubDomainInterfacesIterator;

  //! The type of the iterators over the codim 1 interfaces between all subdomains on a level view.
  typedef typename Traits::LevelAllSubDomainInterfacesIterator LevelAllSubDomainInterfacesIterator;

  /** @name Constructors */
  /*@{*/
  //! Constructs a new MultiDomainGrid from the given host grid.
  /**
   *
   * \param hostGrid                the host grid that will be wrapped by the MultiDomainGrid
   * \param supportLevelIndexSets   flag indicating support for level index sets on subdomains
   */
  explicit MultiDomainGrid(HostGrid& hostGrid, bool supportLevelIndexSets = true) :
    _hostGrid(hostGrid),
    _leafIndexSet(*this,hostGrid.leafView()),
    _globalIdSet(*this),
    _localIdSet(*this),
    _state(stateFixed),
    _adaptState(stateFixed),
    _supportLevelIndexSets(supportLevelIndexSets),
    _maxAssignedSubDomainIndex(0)
  {
    updateIndexSets();
  }
  /*@}*/

  /** @name Dune grid interface methods */
  /*@{*/

  //! The current maximum level of the grid.
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
  /*@}*/

  /** @name Methods for iterating over subdomain interfaces */
  /*@{*/
  //! Returns an iterator over the leaf interface of two subdomains.
  /**
   * The resulting iterator will visit all cell intersections that are part of both subdomains.
   *
   * \attention The iterator assumes the two subdomains to be non-overlapping! If there is an overlap,
   * some intersections will be iterated over twice!
   *
   * \param subDomain1 the first subdomain
   * \param subDomain2 the second subdomain
   */
  LeafSubDomainInterfaceIterator leafSubDomainInterfaceBegin(SubDomainIndexType subDomain1, SubDomainIndexType subDomain2) const {
    return LeafSubDomainInterfaceIterator(*this,subDomain1,subDomain2);
  }

  //! Returns the corresponding end iterator for leafSubDomainInterfaceBegin().
  LeafSubDomainInterfaceIterator leafSubDomainInterfaceEnd(SubDomainIndexType subDomain1, SubDomainIndexType subDomain2) const {
    return LeafSubDomainInterfaceIterator(*this,subDomain1,subDomain2,true);
  }

  //! Returns an iterator over the interface of two subdomains at the given level.
  /**
   * The resulting iterator will visit all cell intersections that are part of both subdomains.
   *
   * \attention The iterator assumes the two subdomains to be non-overlapping! If there is an overlap,
   * some intersections will be iterated over twice!
   *
   * \param subDomain1 the first subdomain
   * \param subDomain2 the second subdomain
   * \param level      the grid level over which to iterate
   */
  LevelSubDomainInterfaceIterator levelSubDomainInterfaceBegin(SubDomainIndexType subDomain1, SubDomainIndexType subDomain2, int level) const {
    return LevelSubDomainInterfaceIterator(*this,subDomain1,subDomain2,level);
  }

  //! Returns the corresponding end iterator for levelSubDomainInterfaceBegin().
  /**
   * \param level the grid level to be iterated over.
   */
  LevelSubDomainInterfaceIterator levelSubDomainInterfaceEnd(SubDomainIndexType subDomain1, SubDomainIndexType subDomain2, int level) const {
    return LevelSubDomainInterfaceIterator(*this,subDomain1,subDomain2,level,true);
  }

  //! Returns an iterator over all subdomain interfaces on the leaf view.
  /**
   * This method returns an iterator that will visit all pairwise surface interfaces
   * between subdomains on the leaf view. In particular, given to adjacent cells
   * \f$e_1\f$ and \f$e_2\f$ and two subdomains \f$s_1\f$ and \f$s_2\f$, this iterator
   * will visit the intersection between \f$e_1\f$ and \f$e_2\f$ iff all of the following hold:
   * \f{eqnarray*} e_1 \in s_1,\ e_1 \not\in s_2,\\ e_2 \not\in s_1,\ e_1 \in s_2.\f}
   * In essence, the two subdomains have to be locally disjoint on \f$e_1\f$ and \f$e_2\f$.
   *
   * The iterator will only traverse the host grid once for visiting all subdomain interfaces.
   * Incrementing the iterator might thus result in an iterator pointing to the same grid
   * intersection, but to a different pair of subdomains. The subdomains pointed to by ther
   * iterator can be retrieved by calling LeafAllSubDomainInterfacesIterator::subDomain1()
   * and LeafAllSubDomainInterfacesIterator::subDomain2(), respectively.
   */
  LeafAllSubDomainInterfacesIterator leafAllSubDomainInterfacesBegin() const {
    return LeafAllSubDomainInterfacesIterator(*this);
  }

  //! Returns the corresponding end iterator for leafAllSubDomainInterfacesBegin().
  LeafAllSubDomainInterfacesIterator leafAllSubDomainInterfacesEnd() const {
    return LeafAllSubDomainInterfacesIterator(*this,true);
  }

  //! Returns an iterator over all subdomain interfaces on the requested level view.
  /**
   * This method returns an iterator that will visit all pairwise surface interfaces
   * between subdomains on the requested level view. In particular, given to adjacent cells
   * \f$e_1\f$ and \f$e_2\f$ and two subdomains \f$s_1\f$ and \f$s_2\f$, this iterator
   * will visit the intersection between \f$e_1\f$ and \f$e_2\f$ iff all of the following hold:
   * \f{eqnarray*} e_1 \in s_1,\ e_1 \not\in s_2,\\ e_2 \not\in s_1,\ e_1 \in s_2.\f}
   * In essence, the two subdomains have to be locally disjoint on \f$e_1\f$ and \f$e_2\f$.
   *
   * The iterator will only traverse the host grid once for visiting all subdomain interfaces.
   * Incrementing the iterator might thus result in an iterator pointing to the same grid
   * intersection, but to a different pair of subdomains. The subdomains pointed to by ther
   * iterator can be retrieved by calling LevelAllSubDomainInterfacesIterator::subDomain1()
   * and LevelAllSubDomainInterfacesIterator::subDomain2(), respectively.
   *
   * \param level the grid level to be iterated over.
   */
  LevelAllSubDomainInterfacesIterator levelAllSubDomainInterfacesBegin(int level) const {
    return LevelAllSubDomainInterfacesIterator(*this,level);
  }

  //! Returns the corresponding end iterator for levelAllSubDomainInterfacesBegin().
  /**
   * \param level the grid level to be iterated over.
   */
  LevelAllSubDomainInterfacesIterator levelAllSubDomainInterfacesEnd(int level) const {
    return LevelAllSubDomainInterfacesIterator(*this,level,true);
  }
  /*@}*/

  /** @name Dune grid interface methods */
  /*@{*/
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

  const typename Traits::LevelIndexSet& levelIndexSet(unsigned int level) const {
    if (!_supportLevelIndexSets) {
      DUNE_THROW(GridError,"level index set support not enabled for this grid");
    }
    assert(level <= maxLevel());
    return *_levelIndexSets[level];
  }

  const typename Traits::LeafIndexSet& leafIndexSet() const {
    return _leafIndexSet;
  }

  void globalRefine(int refCount) {
    saveMultiDomainState();
    _hostGrid.globalRefine(refCount);
    updateIndexSets();
    restoreMultiDomainState();
  }

  bool mark(int refCount, const typename Traits::template Codim<0>::Entity& e) {
    assert(_state == stateFixed);
    return _hostGrid.mark(refCount, hostEntity(e));
  }

  int getMark(const typename Traits::template Codim<0>::Entity& e) {
    assert(_state == stateFixed);
    return _hostGrid.getMark(hostEntity(e));
  }

  bool preAdapt() {
    assert(_state == stateFixed && _adaptState == stateFixed);
    _adaptState = statePreAdapt;
    bool result = _hostGrid.preAdapt();
    return result;
  }

  bool adapt() {
    assert(_state == stateFixed && _adaptState == statePreAdapt);
    _adaptState = statePostAdapt;
    saveMultiDomainState();
    bool result = _hostGrid.adapt();
    updateIndexSets();
    restoreMultiDomainState();
    return result;
  }

  void postAdapt() {
    assert(_state == stateFixed && _adaptState == statePostAdapt);
    _adaptState = stateFixed;
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

  template<typename DataHandleImp, typename DataTypeImp>
  void communicate (CommDataHandleIF<DataHandleImp,DataTypeImp> &data,
                    InterfaceType iftype,
                    CommunicationDirection dir,
                    int level) const
  {
    DataHandleWrapper<CommDataHandleIF<DataHandleImp,DataTypeImp> > datahandle(data,*this);
    _hostGrid.communicate(datahandle,iftype,dir,level);
  }

  template<typename DataHandleImp, typename DataTypeImp>
  void communicate (CommDataHandleIF<DataHandleImp,DataTypeImp> &data,
                    InterfaceType iftype,
                    CommunicationDirection dir) const
  {
    DataHandleWrapper<CommDataHandleIF<DataHandleImp,DataTypeImp> > datahandle(data,*this);
    _hostGrid.communicate(datahandle,iftype,dir);
  }

  size_t numBoundarySegments() const
  {
    return _hostGrid.numBoundarySegments();
  }

  /*@}*/

  /** @name Subdomain creation- and adaptation methods */
  /*@{*/
  //! Prepares the grid for (re-)assigning cells to subdomains.
  /**
   * After calling this method, it becomes possible to invoke the various methods
   * for cell assignment to subdomains. When you are done marking, call preUpdateSubDomains().
   *
   * IMPORTANT: Reassigning subdomains and grid adaptation are mutually exclusive,
   * it is not possibly to do both at the same time. This restriction is enforced
   * by the grid.
   */
  void startSubDomainMarking() {
    assert(_state == stateFixed && _adaptState == stateFixed);
    _tmpLeafIndexSet.reset(new LeafIndexSetImp(_leafIndexSet));
    _tmpLeafIndexSet->reset(false);
    _state = stateMarking;
  }

  //! Calculates the new subdomain layout, but does not update the current subdomains yet.
  /**
   * After calling this method, you can query the grid for the changes that will occur when the
   * new subdomain layout becomes active. This includes the possibility to obtain the new indices
   * entities will be assigned in the modified subdomains.
   *
   * To switch the grid over to the new layout, call updateSubDomains().
   */
  void preUpdateSubDomains() {
    assert(_state == stateMarking && _adaptState == stateFixed);
    if (_supportLevelIndexSets) {
      for (unsigned int l = 0; l <= maxLevel(); ++l) {
        _tmpLevelIndexSets.push_back(make_shared_ptr(new LevelIndexSetImp(*this,_hostGrid.levelView(l))));
      }
    }
    _tmpLeafIndexSet->update(_tmpLevelIndexSets,true);
    _state = statePreUpdate;
  }

  //! Switches the subdomain layout over to the new layout.
  /**
   *
   */
  void updateSubDomains() {
    assert(_state == statePreUpdate && _adaptState == stateFixed);
    _leafIndexSet.swap(*_tmpLeafIndexSet);
    if (_supportLevelIndexSets) {
      for (unsigned int l = 0; l <= maxLevel(); ++l) {
        _levelIndexSets[l]->swap(*_tmpLevelIndexSets[l]);
      }
    }
    _state = statePostUpdate;
  }

  //! clears the saved state of the subdomain layout that was active before the last call to updateSubDomains().
  void postUpdateSubDomains() {
    assert(_state == statePostUpdate && _adaptState == stateFixed);
    _tmpLevelIndexSets.clear();
    _tmpLeafIndexSet.reset(NULL);
    _state = stateFixed;
  }

  //! Adds the given leaf entity to the specified subdomain.
  void addToSubDomain(SubDomainIndexType subDomain, const typename Traits::template Codim<0>::Entity& e) {
    assert(_state == stateMarking);
    assert(e.isLeaf());
    _maxAssignedSubDomainIndex = std::max(_maxAssignedSubDomainIndex,subDomain);
    _tmpLeafIndexSet->addToSubDomain(subDomain,e);
  }

  //! Removes the given leaf entity from the specified subdomain.
  void removeFromSubDomain(SubDomainIndexType subDomain, const typename Traits::template Codim<0>::Entity& e) {
    assert(_state == stateMarking);
    assert(e.isLeaf());
    _tmpLeafIndexSet->removeFromSubDomain(subDomain,e);
  }

  //! Assigns the given leaf entity to the specified subdomain, clearing any previous subdomain assignments.
  void assignToSubDomain(SubDomainIndexType subDomain, const typename Traits::template Codim<0>::Entity& e) {
    assert(_state == stateMarking);
    assert(e.isLeaf());
    _maxAssignedSubDomainIndex = std::max(_maxAssignedSubDomainIndex,subDomain);
    _tmpLeafIndexSet->assignToSubDomain(subDomain,e);
  }

  //! Removes the given leaf entity from all subdomains it currently belongs to.
  void removeFromAllSubDomains(const typename Traits::template Codim<0>::Entity& e) {
    assert(_state == stateMarking);
    assert(e.isLeaf());
    _tmpLeafIndexSet->removeFromAllSubDomains(e);
  }
  /*@}*/

  /** @name Access to the subdomain grids */
  /*@{*/
  //! Returns a reference to the SubDomainGrid associated with the given subdomain.
  const SubDomainGrid& subDomain(SubDomainIndexType subDomain) const {
    boost::shared_ptr<SubDomainGrid>& subGridPointer = _subDomainGrids[subDomain];
    if (!subGridPointer) {
      subGridPointer.reset(new SubDomainGrid(const_cast<MultiDomainGrid&>(*this),subDomain));
      // subGridPointer->update();
    }
    return *subGridPointer;
  }

  //! Returns a reference to the SubDomainGrid associated with the given subdomain.
  SubDomainGrid& subDomain(SubDomainIndexType subDomain) {
    boost::shared_ptr<SubDomainGrid>& subGridPointer = _subDomainGrids[subDomain];
    if (!subGridPointer) {
      subGridPointer.reset(new SubDomainGrid(*this,subDomain));
      // subGridPointer->update();
    }
    return *subGridPointer;
  }

  //! Returns the largest subdomain index that was ever assigned to a cell in this grid.
  /**
   * This method returns the largest subdomain index that was passed to addToSubDomain() or
   * assignToSubDomain() since this MultiDomainGrid was created. Keep in mind that the subdomain
   * belonging to that index might not contain any entities anymore if all entities have been
   * removed from it at a later point.
   */
  SubDomainIndexType maxAssignedSubDomainIndex() const
  {
    return _maxAssignedSubDomainIndex;
  }

  //! Indicates whether this MultiDomainGrid instance supports level index sets on its SubDomainGrids.
  bool supportLevelIndexSets() const {
    return _supportLevelIndexSets;
  }
  /*@}*/

  /** @name Entity conversion methods */
  /*@{*/
  //! Returns a reference to the corresponding host entity.
  /**
   * \warning The returned reference will only be valid as long as the passed-in reference to the
   * MultiDomainGrid entity! If you need a persistent reference to the host entity, use
   * hostEntityPointer() instead.
   */
  template<typename EntityType>
  const typename HostEntity<EntityType>::type& hostEntity(const EntityType& e) const {
    return *(getRealImplementation(e).hostEntityPointer());
  }

  //! Returns an EntityPointer to the corresponding host entity.
  template<typename EntityType>
  const typename HostEntityPointer<EntityType>::type hostEntityPointer(const EntityType& e) const {
    return getRealImplementation(e).hostEntityPointer();
  }

  template<typename EntityType>
  const typename MultiDomainEntity<EntityType>::type& multiDomainEntity(const EntityType& e) const {
    return *(SubDomainGrid::getRealImplementation(e).multiDomainEntityPointer());
  }

  //! Returns an EntityPointer to the corresponding MultiDomain entity.
  template<typename EntityType>
  const typename MultiDomainEntityPointer<EntityType>::type multiDomainEntityPointer(const EntityType& e) const {
    return SubDomainGrid::getRealImplementation(e).multiDomainEntityPointer();
  }
 /*@}*/

private:

  HostGrid& _hostGrid;

  std::vector<boost::shared_ptr<LevelIndexSetImp> > _levelIndexSets;
  LeafIndexSetImp _leafIndexSet;

  std::vector<boost::shared_ptr<LevelIndexSetImp> > _tmpLevelIndexSets;
  boost::scoped_ptr<LeafIndexSetImp> _tmpLeafIndexSet;

  GlobalIdSetImp _globalIdSet;
  LocalIdSetImp _localIdSet;

  State _state;
  State _adaptState;
  const bool _supportLevelIndexSets;

  mutable std::map<SubDomainIndexType,boost::shared_ptr<SubDomainGrid> > _subDomainGrids;
  SubDomainIndexType _maxAssignedSubDomainIndex;

  AdaptationStateMap _adaptationStateMap;

  void updateIndexSets() {
    // make sure we have enough LevelIndexSets
    if (_supportLevelIndexSets) {
      while (_levelIndexSets.size() <= maxLevel()) {
        _levelIndexSets.push_back(make_shared_ptr(new LevelIndexSetImp(*this,_hostGrid.levelView(_levelIndexSets.size()))));
      }
      // and make sure we don't have too many...
      if (_levelIndexSets.size() > maxLevel() + 1)
        {
          _levelIndexSets.resize(maxLevel() + 1);
        }
    }

    _leafIndexSet.reset(true);
    _leafIndexSet.update(_levelIndexSets,true);

    _globalIdSet.update(_hostGrid.globalIdSet());
    _localIdSet.update(_hostGrid.localIdSet());
  }

  void saveMultiDomainState() {
    typedef typename ThisType::LeafGridView GV;
    GV gv = this->leafView();
    typedef typename GV::template Codim<0>::Iterator Iterator;
    typedef typename GV::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GV::template Codim<0>::Entity Entity;
    typedef typename MDGridTraits::template Codim<0>::SubDomainSet SubDomainSet;
    for (Iterator it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
      const Entity& e = *it;
      const SubDomainSet& subDomains = gv.indexSet().subDomains(e);
      _adaptationStateMap[localIdSet().id(e)] = subDomains;
      EntityPointer ep(e);
      while (ep->mightVanish()) {
        ep = ep->father();
        typename Traits::LocalIdSet::IdType id = localIdSet().id(*ep);
        // if the entity has not been added to the set, create a new entry
        // (entity returns false and does not change the map if there is already an entry for "id")
        if (!_adaptationStateMap.insert(typename AdaptationStateMap::value_type(id,subDomains)).second) {
          // otherwise add the leaf entity's subdomains to the existing set of subdomains
          _adaptationStateMap[id].addAll(subDomains);
        }
      }
    }
  }

  void restoreMultiDomainState() {
    typedef typename ThisType::LeafGridView GV;
    GV gv = this->leafView();
    typedef typename GV::template Codim<0>::Iterator Iterator;
    typedef typename GV::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GV::template Codim<0>::Entity Entity;
    typedef typename MDGridTraits::template Codim<0>::SubDomainSet SubDomainSet;
    for (Iterator it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
      EntityPointer ep(it);
      // First try to exploit the information in the underlying grid
      while (ep->isNew()) {
        ep = ep->father();
      }
      // This might not work, as there are no isNew() marks for globalrefine()
      // We thus have to look up the former leaf entity in our adaptation map
      typename AdaptationStateMap::iterator asmit = _adaptationStateMap.find(localIdSet().id(*ep));
      while(asmit == _adaptationStateMap.end()) {
        ep = ep->father();
        asmit = _adaptationStateMap.find(localIdSet().id(*ep));
      }
      _leafIndexSet.addToSubDomains(asmit->second, *it);
    }
    _leafIndexSet.update(_levelIndexSets,false);
    _adaptationStateMap.clear();
  }

  template<typename GridView, typename HostGridView>
  typename GridView::IntersectionIterator multiDomainIntersectionIterator(typename HostGridView::IntersectionIterator iit) const {
    return typename GridView::IntersectionIterator::ImplementationType(iit);
  }

  template<typename Entity>
  typename Traits::template Codim<Entity::codimension>::EntityPointer wrapHostEntity(const Entity& e) const {
    return wrapHostEntity<Entity::codimension>(e);
  }

  template<int codim>
  typename Traits::template Codim<codim>::EntityPointer wrapHostEntity(const typename HostGrid::template Codim<codim>::Entity& e)
  {
    return EntityPointerWrapper<codim,GridImp>(e);
  }


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
      return _impl.fixedsize(dim,codim);
    }

    template<typename Entity>
    std::size_t size(const Entity& e) const
    {
      return _impl.size(*_grid.wrapHostEntity(e));
    }

    template<typename MessageBufferImp, typename Entity>
    void gather(MessageBufferImp& buf, const Entity& e) const
    {
      _impl.gather(buf,*_grid.wrapHostEntity(e));
    }

    template<typename MessageBufferImp, typename Entity>
    void scatter(MessageBufferImp& buf, const Entity& e, std::size_t n)
    {
      _impl.scatter(buf,*_grid.wrapHostEntity(e),n);
    }

    DataHandleWrapper(Impl& impl, const MultiDomainGrid<HostGrid,MDGridTraitsType>& grid)
      : _impl(impl)
      , _grid(grid)
    {}

    Impl& _impl;
    const MultiDomainGrid<HostGrid,MDGridTraitsType>& _grid;

  };

};

enum MultiDomainGridType { multiDomainGrid, subDomainGrid, other };

template<typename T>
struct GridType {
  static const MultiDomainGridType v = other;
};

template<class HostGrid, typename MDGridTraits>
struct GridType<MultiDomainGrid<HostGrid,MDGridTraits> > {
  static const MultiDomainGridType v = multiDomainGrid;
};

template<class MDGrid>
struct GridType<subdomain::SubDomainGrid<MDGrid> > {
  static const MultiDomainGridType v = subDomainGrid;
};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_MULTIDOMAINGRID_HH
