#ifndef DUNE_MULTIDOMAINGRID_MULTIDOMAINGRID_HH
#define DUNE_MULTIDOMAINGRID_MULTIDOMAINGRID_HH

#include <string>
#include <boost/shared_ptr.hpp>

#include <dune/grid/multidomaingrid/subdomainset.hh>

#include <dune/grid/multidomaingrid/geometry.hh>
#include <dune/grid/multidomaingrid/entity.hh>
#include <dune/grid/multidomaingrid/entitypointer.hh>
#include <dune/grid/multidomaingrid/leafiterator.hh>
#include <dune/grid/multidomaingrid/leveliterator.hh>
#include <dune/grid/multidomaingrid/hierarchiciterator.hh>
#include <dune/grid/multidomaingrid/intersectioniterator.hh>
#include <dune/grid/multidomaingrid/idsets.hh>
#include <dune/grid/multidomaingrid/indexsets.hh>


namespace Dune {

namespace mdgrid {

template<typename HostGrid>
class MultiDomainGrid;

template<typename HostGrid>
struct MultiDomainGridFamily {

  typedef GridTraits<
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
  friend struct detail::HostGridTraits;

  template<typename GridImp>
  friend class LeafIntersectionIteratorWrapper;

  template<typename GridImp>
  friend class LevelIntersectionIteratorWrapper;

  typedef MultiDomainGrid<HostGrid> GridImp;
  typedef HostGrid HostGridType;

  typedef IndexSetWrapper<const MultiDomainGrid<HostGrid>, typename HostGrid::LevelGridView> LevelIndexSetImp;

  typedef IndexSetWrapper<const MultiDomainGrid<HostGrid>, typename HostGrid::LeafGridView> LeafIndexSetImp;

  typedef IdSetWrapper<const MultiDomainGrid<HostGrid>, typename HostGrid::Traits::GlobalIdSet> GlobalIdSetImp;

  typedef IdSetWrapper<const MultiDomainGrid<HostGrid>, typename HostGrid::Traits::LocalIdSet> LocalIdSetImp;

  enum State { fixed, marking, preUpdate, postUpdate };

public:

  typedef MultiDomainGridFamily<HostGrid> GridFamily;
  typedef typename GridFamily::Traits Traits;
  typedef typename HostGrid::ctype ctype;
  typedef IntegralTypeSubDomainSet<3> SubDomainSet;
  typedef typename SubDomainSet::DomainType SubDomainType;

  explicit MultiDomainGrid(HostGrid& hostGrid) :
    _hostGrid(hostGrid),
    _leafIndexSet(*this,hostGrid.leafView()),
    _globalIdSet(*this),
    _localIdSet(*this)
  {
    updateIndexSets();
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
    _tmpLeafIndexSet.update(_tmpLevelIndexSets,false);
    _state = preUpdate;
  }

  void updateSubDomains() {
    assert(_state == preUpdate);
    _leafIndexSet.swap(_tmpLeafIndexSet);
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
    _tmpLeafIndexSet.addToSubDomain(subDomain,e);
  }


  void removeFromSubDomain(SubDomainType subDomain, const typename Traits::template Codim<0>::Entity& e) {
    assert(_state == marking);
    assert(e.isLeaf());
    _tmpLeafIndexSet.removeFromSubDomain(subDomain,e);
  }


  void assignToSubDomain(SubDomainType subDomain, const typename Traits::template Codim<0>::Entity& e) {
    assert(_state == marking);
    assert(e.isLeaf());
    _tmpLeafIndexSet.assignToSubDomain(subDomain,e);
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

  template<typename Entity>
  struct HostEntity {
    typedef typename HostGrid::Traits::template Codim<Entity::codimension>::Entity type;
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

  void updateIndexSets() {
    // make sure we have enough LevelIndexSets
    while (_levelIndexSets.size() <= maxLevel()) {
      _levelIndexSets.push_back(new LevelIndexSetImp(*this,_hostGrid.levelView(_levelIndexSets.size())));
    }

    _leafIndexSet.update(_levelIndexSets,true);

    _globalIdSet.update(_hostGrid.globalIdSet());
    _localIdSet.update(_hostGrid.localIdSet());
  }

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_MULTIDOMAINGRID_HH
