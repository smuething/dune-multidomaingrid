#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINTOSUBDOMAININTERFACEITERATOR_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINTOSUBDOMAININTERFACEITERATOR_HH

#include <dune/grid/multidomaingrid/hostgridaccessor.hh>
#include <dune/grid/multidomaingrid/subdomaininterfaceiterator.hh>

namespace Dune {

namespace mdgrid {

template<typename GridImp>
class LeafSubDomainInterfaceIterator;

template<typename GridImp>
class LevelSubDomainInterfaceIterator;

template<typename SubDomainIndex>
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
        const auto& hostIntersection = *it._hostIntersectionIterator;
        if (hostIntersection.neighbor() && it._gridView.indexSet().containsForSubDomain(_subDomain2,hostIntersection.outside())) {
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

  SubDomainToSubDomainController(SubDomainIndex subDomain1, SubDomainIndex subDomain2)
    : _subDomain1(subDomain1)
    , _subDomain2(subDomain2)
  {}

  SubDomainIndex subDomain1() const
  {
    return _subDomain1;
  }

  SubDomainIndex subDomain2() const
  {
    return _subDomain2;
  }

  const SubDomainIndex _subDomain1;
  const SubDomainIndex _subDomain2;
};


template<typename GridImp>
class LeafSubDomainInterfaceIterator :
    public SubDomainInterfaceIterator<GridImp,
                                      typename GridImp::LeafGridView,
                                      typename detail::HostGridAccessor<GridImp>::Type::LeafGridView,
                                      SubDomainToSubDomainController<typename GridImp::SubDomainIndex>
                                      >
{

  template<typename, typename, typename, typename>
  friend class SubDomainInterfaceIterator;

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename,typename>
  friend class MultiDomainGrid;

  typedef SubDomainToSubDomainController<typename GridImp::SubDomainIndex> Controller;

  typedef SubDomainInterfaceIterator<GridImp,
                                     typename GridImp::LeafGridView,
                                     typename GridImp::HostGrid::LeafGridView,
                                     Controller
                                     > Base;

  typedef typename Base::Intersection::SubDomainIndex SubDomainIndex;

  LeafSubDomainInterfaceIterator(const GridImp& grid, SubDomainIndex subDomain1, SubDomainIndex subDomain2, bool end=false) :
    Base(grid.leafGridView(),grid._hostGrid.leafGridView(),Controller(subDomain1,subDomain2),end)
  {}

};


template<typename GridImp>
class LevelSubDomainInterfaceIterator :
    public SubDomainInterfaceIterator<GridImp,
                                      typename GridImp::LevelGridView,
                                      typename detail::HostGridAccessor<GridImp>::Type::LevelGridView,
                                      SubDomainToSubDomainController<typename GridImp::SubDomainIndex>
                                      >
{

  template<typename, typename, typename, typename>
  friend class SubDomainInterfaceIterator;

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename,typename>
  friend class MultiDomainGrid;

  typedef SubDomainToSubDomainController<typename GridImp::SubDomainIndex> Controller;

  typedef SubDomainInterfaceIterator<GridImp,
                                     typename GridImp::LevelGridView,
                                     typename GridImp::HostGrid::LevelGridView,
                                     Controller
                                     > Base;

  typedef typename Base::Intersection::SubDomainIndex SubDomainIndex;

  LevelSubDomainInterfaceIterator(const GridImp& grid, SubDomainIndex subDomain1, SubDomainIndex subDomain2, int level, bool end=false) :
    Base(grid.levelGridView(level),grid._hostGrid.levelGridView(level),Controller(subDomain1,subDomain2),end)
  {}

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINTOSUBDOMAININTERFACEITERATOR_HH
