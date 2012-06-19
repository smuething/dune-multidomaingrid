#ifndef DUNE_MULTIDOMAINGRID_ALLSUBDOMAININTERFACESITERATOR_HH
#define DUNE_MULTIDOMAINGRID_ALLSUBDOMAININTERFACESITERATOR_HH

#include <dune/common/fvector.hh>
#include <dune/grid/multidomaingrid/subdomaininterfaceiterator.hh>

namespace Dune {

namespace mdgrid {


template<typename GridImp>
class LeafAllSubDomainInterfacesIterator;

template<typename GridImp>
class LevelAllSubDomainInterfacesIterator;

template<typename SubDomainSet>
class AllInterfacesController
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
  friend class LeafAllSubDomainInterfacesIterator;

  template<typename GridImp>
  friend class LevelAllSubDomainInterfacesIterator;

  typedef typename SubDomainSet::SubDomainIndexType SubDomainIndexType;
  typedef typename SubDomainSet::Iterator SubDomainIterator;

  template<typename Iterator>
  bool calculateInterfacingSubDomains(Iterator& it)
  {
    const typename Iterator::HostIntersectionIterator::Intersection::EntityPointer outside = it._hostIntersectionIterator->outside();
    const SubDomainSet& subDomains2 = it._gridView.indexSet().subDomainsForHostEntity(*outside);
    _interfacingSubDomains1.difference(*_subDomains1,subDomains2);
    _interfacingSubDomains2.difference(subDomains2,*_subDomains1);
    _subDomain1Iterator = _interfacingSubDomains1.begin();
    _subDomain1End = _interfacingSubDomains1.end();
    _subDomain2Iterator = _interfacingSubDomains2.begin();
    _subDomain2End = _interfacingSubDomains2.end();
    return !(_interfacingSubDomains1.empty() || _interfacingSubDomains2.empty());
  }

  template<typename Iterator>
  void incrementToNextValidPosition(Iterator& it) {
    if (_subDomain2Iterator != _subDomain2End) {
      return;
    }
    ++_subDomain1Iterator;
    if (_subDomain1Iterator != _subDomain1End) {
      _subDomain2Iterator = _interfacingSubDomains2.begin();
      return;
    }
    incrementToNextValidIntersection(it);
  }

  template<typename Iterator>
  bool incrementToNextValidIntersection(Iterator& it)
  {
    for (;;) {
      do {
        ++it._hostIntersectionIterator;
        if (it._hostIntersectionIterator->neighbor() &&
            calculateInterfacingSubDomains(it))
          return true;
      } while (it._hostIntersectionIterator != it._hostIntersectionEnd);
      ++it._hostIterator;
      if (it._hostIterator == it._hostEnd)
        return false;
      _subDomains1 = &it._gridView.indexSet().subDomainsForHostEntity(*it._hostIterator);
      it._hostIntersectionIterator = it._hostGridView.ibegin(*it._hostIterator);
      it._hostIntersectionEnd = it._hostGridView.iend(*it._hostIterator);
      if (it._hostIntersectionIterator->neighbor() &&
          calculateInterfacingSubDomains(it))
        return true;
    }
  }

  template<typename Iterator>
  void increment(Iterator& it) {
    ++_subDomain2Iterator;
    incrementToNextValidPosition(it);
  }

  template<typename Iterator>
  void incrementToStartPosition(Iterator& it)
  {
    if (it._hostIterator != it._hostEnd) {
      _subDomains1 = &it._gridView.indexSet().subDomainsForHostEntity(*it._hostIterator);
      if (!it._hostIntersectionIterator->neighbor() || !calculateInterfacingSubDomains(it))
        incrementToNextValidIntersection(it);
    }
  }

  AllInterfacesController()
    : _subDomain1Iterator(_interfacingSubDomains1.end())
    , _subDomain2Iterator(_interfacingSubDomains2.end())
    , _subDomain1End(_subDomain1Iterator)
    , _subDomain2End(_subDomain2Iterator)
  {}

  SubDomainIndexType subDomain1() const
  {
    return *_subDomain1Iterator;
  }

  SubDomainIndexType subDomain2() const
  {
    return *_subDomain2Iterator;
  }

  const SubDomainSet* _subDomains1;
  SubDomainSet _interfacingSubDomains1;
  SubDomainSet _interfacingSubDomains2;
  SubDomainIterator _subDomain1Iterator;
  SubDomainIterator _subDomain2Iterator;
  SubDomainIterator _subDomain1End;
  SubDomainIterator _subDomain2End;
};


template<typename GridImp>
class LeafAllSubDomainInterfacesIterator :
    public SubDomainInterfaceIterator<GridImp,
                                      typename GridImp::LeafGridView,
                                      typename detail::HostGridAccessor<GridImp>::Type::LeafGridView,
                                      AllInterfacesController<typename GridImp::MDGridTraits::template Codim<0>::SubDomainSet>
                                      >
{

  template<typename, typename, typename, typename>
  friend class SubDomainInterfaceIterator;

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename,typename>
  friend class MultiDomainGrid;

  typedef AllInterfacesController<typename GridImp::MDGridTraits::template Codim<0>::SubDomainSet> Controller;

  typedef SubDomainInterfaceIterator<GridImp,
                                     typename GridImp::LeafGridView,
                                     typename GridImp::HostGridType::LeafGridView,
                                     Controller
                                     > Base;

  LeafAllSubDomainInterfacesIterator(const GridImp& grid, bool end=false) :
    Base(grid.leafView(),grid._hostGrid.leafView(),Controller(),end)
  {}

};


template<typename GridImp>
class LevelAllSubDomainInterfacesIterator :
    public SubDomainInterfaceIterator<GridImp,
                                      typename GridImp::LevelGridView,
                                      typename detail::HostGridAccessor<GridImp>::Type::LevelGridView,
                                      AllInterfacesController<typename GridImp::MDGridTraits::template Codim<0>::SubDomainSet>
                                      >
{

  template<typename, typename, typename, typename>
  friend class SubDomainInterfaceIterator;

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename,typename>
  friend class MultiDomainGrid;

  typedef AllInterfacesController<typename GridImp::MDGridTraits::template Codim<0>::SubDomainSet> Controller;

  typedef SubDomainInterfaceIterator<GridImp,
                                     typename GridImp::LevelGridView,
                                     typename GridImp::HostGridType::LevelGridView,
                                     Controller
                                     > Base;

  LevelAllSubDomainInterfacesIterator(const GridImp& grid, int level, bool end=false) :
    Base(grid.levelView(level),grid._hostGrid.levelView(level),Controller(),end)
  {}

};


} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_ALLSUBDOMAININTERFACESITERATOR_HH
