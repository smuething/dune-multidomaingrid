#ifndef DUNE_MULTIDOMAINGRID_MULTIDOMAINGRID_HH
#define DUNE_MULTIDOMAINGRID_MULTIDOMAINGRID_HH

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
    IndexSetWrapper<const MultiDomainGrid<HostGrid>, typename HostGrid::Traits::LevelIndexSet>,
    IndexSetWrapper<const MultiDomainGrid<HostGrid>, typename HostGrid::Traits::LeafIndexSet>,
    IdSetWrapper<const MultiDomainGrid<HostGrid>, typename HostGrid::Traits::GlobalIdSet>,
    typename HostGrid::Traits::GlobalIdSet::IdType,
    IdSetWrapper<const MultiDomainGrid<HostGrid>, typename HostGrid::Traits::LocalIdSet>,
    typename HostGrid::Traits::LocalIdSet::IdType,
    CollectiveCommunication<MultiDomainGrid<HostGrid> >
    > Traits;

};

template<typename HostGrid>
class MultiDomainGrid :
    public GridDefaultImplementation<HostGrid::dimension,
				     HostGrid::dimensionworld,
				     HostGrid::ctype,
				     MultiDomainGridFamily<HostGrid> > {

public:

  typedef MultiDomainGridFamily<HostGrid> GridFamily;
  typedef GridFamily::Traits Traits;
  typedef typename HostGrid::ctype ctype;

  explicit MultiDomainGrid(HostGrid& hostGrid) :
    _hostGrid(hostGrid)
  {}

  std::string name() const {
    return "MultiDomainGrid";
  }

  std::size_t maxLevel() const {
    return _hostGrid->maxLevel();
  }

  template<int codim>
  typename Traits::template Codim<codim>::LevelIterator lbegin(int level) const {
    return ...;
  }

  template<int codim>
  typename Traits::template Codim<codim>::LevelIterator lend(int level) const {
    return ...;
  }

  template<int codim, PartitionIteratorType PiType>
  typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin(int level) const {
    return ...;
  }

  template<int codim, PartitionIteratorType PiType>
  typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend(int level) const {
    return ...;
  }

  template<int codim>
  typename Traits::template Codim<codim>::LeafIterator leafbegin() const {
    return ...;
  }

  template<int codim>
  typename Traits::template Codim<codim>::LeafIterator leafEnd() const {
    return ...;
  }

  template<int codim, PartitionIteratorType PiType>
  typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const {
    return ...;
  }

  template<int codim, PartitionIteratorType PiType>
  typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const {
    return ...;
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
    return ...;
  }

  const typename Traits::LocalIdSet& localIdSet() const {
    return ...;
  }

  const typename Traits::LevelIndexSet& levelIndexSet(int level) const {
    return ...;
  }

  const typename Traits::LeafIndexSet& leafIndexSet() const {
    return ...;
  }

  void globalRefine(int refCount) {
    _hostGrid.globalRefine(refCount);
  }

  bool mark(int refCount, const typename Traits::template Codim<0>::Entity& e) {
    return _hostGrid.mark(refCount, getHostEntity(e));
  }

  int getMark(const typename Traits::template Codim<0>::Entity& e) {
    return _hostGrid.getMark(getHostEntity(e));
  }

  bool preAdapt() {
    return _hostGrid.preAdapt();
  }

  bool adapt() {
    return _hostGrid.adapt();
  }

  void postAdapt() {
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

private:

  HostGrid& _hostGrid;

  template<typename Entity>
  struct HostEntity {
    typedef typename HostGrid::Traits::template Codim<Entity::codimension>::Entity type;
  };

  template<typename EntityType>
  typename HostEntity<EntityType>::type& hostEntity(EntityType& e) const {
    return e.wrappedEntity();
  }

  template<typename EntityType>
  typename const HostEntity<EntityType>::type& hostEntity(const EntityType& e) const {
    return e.wrappedEntity();
  }

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_MULTIDOMAINGRID_HH
