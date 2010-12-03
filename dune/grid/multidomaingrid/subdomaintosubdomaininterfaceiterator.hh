#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINTOSUBDOMAININTERFACEITERATOR_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINTOSUBDOMAININTERFACEITERATOR_HH

#include <dune/grid/multidomaingrid/subdomaininterfaceiterator.hh>

namespace Dune {

namespace mdgrid {

template<typename GridImp>
class LeafSubDomainInterfaceIterator;

template<typename GridImp>
class LevelSubDomainInterfaceIterator;

template<typename SubDomainIndexType>
class SubDomainToSubDomainController
{

  template<typename GridImp,
           typename GridView,
           typename HostGridView,
           typename IntersectionController
           >
  friend class SubDomainInterface;

  template<typename GridImp,
           typename GridView,
           typename HostGridView,
           typename IntersectionController
           >
  friend class SubDomainInterfaceIterator;

  template<typename GridImp>
  friend class LeafSubDomainInterfaceIterator;

  template<typename GridImp>
  friend class LevelSubDomainInterfaceIterator;

  template<typename Iterator>
  bool incrementToNextValidEntity(Iterator& it) {
    while (it._hostIterator != it._hostEnd) {
      if (it._gridView.indexSet().containsForSubDomain(_subDomain1,*it._hostIterator)) {
        it._hostIntersectionIterator = it._hostGridView.ibegin(*it._hostIterator);
        it._hostIntersectionEnd = it._hostGridView.iend(*it._hostIterator);
        return true;
      }
      ++it._hostIterator;
    }
    return false;
  }

  template<typename Iterator>
  void incrementToNextValidPosition(Iterator& it) {
    for (;;) {
      while(it._hostIntersectionIterator != it._hostIntersectionEnd) {
        if (it._hostIntersectionIterator->neighbor() && it._gridView.indexSet().containsForSubDomain(_subDomain2,*it._hostIntersectionIterator->outside())) {
          return;
        }
        ++it._hostIntersectionIterator;
      }
      ++it._hostIterator;
      if (!incrementToNextValidEntity(it)) {
        return;
      }
    }
  }

  template<typename Iterator>
  void increment(Iterator& it) {
    ++it._hostIntersectionIterator;
    incrementToNextValidPosition(it);
  }

  template<typename Iterator>
  void incrementToStartPosition(Iterator& it)
  {
    if (incrementToNextValidEntity(it)) {
      incrementToNextValidPosition(it);
    }
  }

  SubDomainToSubDomainController(SubDomainIndexType subDomain1, SubDomainIndexType subDomain2)
    : _subDomain1(subDomain1)
    , _subDomain2(subDomain2)
  {}

  SubDomainIndexType subDomain1() const
  {
    return _subDomain1;
  }

  SubDomainIndexType subDomain2() const
  {
    return _subDomain2;
  }

  const SubDomainIndexType _subDomain1;
  const SubDomainIndexType _subDomain2;
};


template<typename GridImp>
class LeafSubDomainInterfaceIterator :
    public SubDomainInterfaceIterator<GridImp,
                                      typename GridImp::LeafGridView,
                                      typename detail::HostGridAccessor<GridImp>::Type::LeafGridView,
                                      SubDomainToSubDomainController<typename GridImp::SubDomainIndexType>
                                      >
{

  template<typename, typename, typename, typename>
  friend class SubDomainInterfaceIterator;

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename,typename>
  friend class MultiDomainGrid;

  typedef SubDomainToSubDomainController<typename GridImp::SubDomainIndexType> Controller;

  typedef SubDomainInterfaceIterator<GridImp,
                                     typename GridImp::LeafGridView,
                                     typename GridImp::HostGridType::LeafGridView,
                                     Controller
                                     > Base;

  typedef typename Base::Intersection::SubDomainIndexType SubDomainIndexType;

  LeafSubDomainInterfaceIterator(const GridImp& grid, SubDomainIndexType subDomain1, SubDomainIndexType subDomain2, bool end=false) :
    Base(grid.leafView(),grid._hostGrid.leafView(),Controller(subDomain1,subDomain2),end)
  {}

};


template<typename GridImp>
class LevelSubDomainInterfaceIterator :
    public SubDomainInterfaceIterator<GridImp,
                                      typename GridImp::LevelGridView,
                                      typename detail::HostGridAccessor<GridImp>::Type::LevelGridView,
                                      SubDomainToSubDomainController<typename GridImp::SubDomainIndexType>
                                      >
{

  template<typename, typename, typename, typename>
  friend class SubDomainInterfaceIterator;

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename,typename>
  friend class MultiDomainGrid;

  typedef SubDomainToSubDomainController<typename GridImp::SubDomainIndexType> Controller;

  typedef SubDomainInterfaceIterator<GridImp,
                                     typename GridImp::LevelGridView,
                                     typename GridImp::HostGridType::LevelGridView,
                                     Controller
                                     > Base;

  typedef typename Base::Intersection::SubDomainIndexType SubDomainIndexType;

  LevelSubDomainInterfaceIterator(const GridImp& grid, SubDomainIndexType subDomain1, SubDomainIndexType subDomain2, int level, bool end=false) :
    Base(grid.levelView(level),grid._hostGrid.levelView(level),Controller(subDomain1,subDomain2),end)
  {}

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINTOSUBDOMAININTERFACEITERATOR_HH
