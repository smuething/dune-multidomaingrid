#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAININTERFACEITERATOR_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAININTERFACEITERATOR_HH

#include <dune/common/iteratorfacades.hh>

namespace Dune {

namespace mdgrid {

template<typename GridImp,
         typename WrapperImp,
         typename GridView,
         typename HostGridView,
         typename IntersectionType>
class SubDomainInterfaceIterator : public ForwardIteratorFacade<SubDomainInterfaceIterator<
                                                                  GridImp,
                                                                  WrapperImp,
                                                                  GridView,
                                                                  HostGridView,
                                                                  IntersectionType
                                                                  >,
                                                                IntersectionType>
{

  typedef typename HostGridView::template Codim<0>::Iterator HostIterator;
  typedef typename HostGridView::IntersectionIterator HostIntersectionIterator;
  typedef IntersectionType Intersection;
  typedef typename GridImp::SubDomainType SubDomainType;

  typedef typename GridImp::Traits::template Codim<0>::EntityPointer EntityPointer;
  typedef typename GridImp::Traits::template Codim<0>::Entity Entity;
  typedef typename GridImp::Traits::template Codim<1>::Geometry Geometry;
  typedef typename GridImp::Traits::template Codim<1>::LocalGeometry LocalGeometry;

  typedef typename GridImp::ctype ctype;
  static const int dimension = GridImp::dimension;
  static const int dimensionworld = GridImp::dimensionworld;

  typedef FieldVector<ctype,dimensionworld> GlobalCoords;
  typedef FieldVector<ctype,dimension - 1> LocalCoords;

protected:

  SubDomainInterfaceIterator(const GridView& gridView,
                             const HostGridView& hostGridView,
                             SubDomainType domain1,
                             SubDomainType domain2,
                             bool end) :
    _gridView(gridView),
    _hostGridView(hostGridView),
    _domain1(domain1),
    _domain2(domain2),
    _hostIterator(end ? hostGridView.template end<0>() : hostGridView.template begin<0>()),
    _hostEnd(hostGridView.template end<0>()),
    _hostIntersectionIterator(hostGridView.ibegin(*hostGridView.template begin<0>())),
    _hostIntersectionEnd(hostGridView.iend(*hostGridView.template begin<0>()))
  {
    if (incrementToNextValidEntity()) {
      incrementToNextValidPosition();
    }
  }

public:

  bool incrementToNextValidEntity() {
    while (_hostIterator != _hostEnd) {
      if (_gridView.indexSet().containsForSubDomain(_domain1,*_hostIterator)) {
        _hostIntersectionIterator = _hostGridView.ibegin(*_hostIterator);
        _hostIntersectionEnd = _hostGridView.iend(*_hostIterator);
        return true;
      }
      ++_hostIterator;
    }
    return false;
  }

  void incrementToNextValidPosition() {
    for (;;) {
      while(_hostIntersectionIterator != _hostIntersectionEnd) {
        if (_hostIntersectionIterator->neighbor() && _gridView.indexSet().containsForSubDomain(_domain2,*_hostIntersectionIterator->outside())) {
          return;
        }
        ++_hostIntersectionIterator;
      }
      ++_hostIterator;
      if (!incrementToNextValidEntity()) {
        return;
      }
    }
  }

  bool equals(const WrapperImp& rhs) const {
    assert(_domain1 == rhs._domain1 && _domain2 == rhs._domain2);
    return _hostIterator == rhs._hostIterator && (_hostIterator == _hostEnd || _hostIntersectionIterator == rhs._hostIntersectionIterator); //TODO: domains?
  }

  bool equals(const SubDomainInterfaceIterator& rhs) const {
    assert(_domain1 == rhs._domain1 && _domain2 == rhs._domain2);
    return _hostIterator == rhs._hostIterator && (_hostIterator == _hostEnd || _hostIntersectionIterator == rhs._hostIntersectionIterator); //TODO: domains?
  }

  void increment() {
    ++_hostIntersectionIterator;
    incrementToNextValidPosition();
    _geometry.clear();
    _geometryInInside.clear();
    _geometryInOutside.clear();
  }

  const Intersection& dereference() const {
    return reinterpret_cast<const Intersection&>(*this);
  }

  const Intersection& operator*() const {
    return dereference();
  }

  const Intersection* operator->() const {
    return &(dereference());
  }

  EntityPointer firstEntity() const {
    return EntityPointerWrapper<0,GridImp>(_hostIntersectionIterator->inside());
  }

  EntityPointer secondEntity() const {
    return EntityPointerWrapper<0,GridImp>(_hostIntersectionIterator->outside());
  }

  bool conforming() const {
    return _hostIntersectionIterator->conforming();
  }

  const LocalGeometry& geometryInFirstEntity() const {
    if (!_geometryInInside.isSet()) {
      _geometryInInside.reset(_hostIntersectionIterator->geometryInInside());
    }
    return _geometryInInside;
  }

  const LocalGeometry& geometryInSecondEntity() const {
    if (!_geometryInOutside.isSet()) {
      _geometryInOutside.reset(_hostIntersectionIterator->geometryInOutside());
    }
    return _geometryInOutside;
  }

  const Geometry& geometry() const {
    if (!_geometry.isSet()) {
      _geometry.reset(_hostIntersectionIterator->geometry());
    }
    return _geometry;
  }

  GeometryType type() const {
    return _hostIntersectionIterator->type();
  }

  int indexInFirstEntity() const {
    return _hostIntersectionIterator->indexInInside();
  }

  int indexInSecondEntity() const {
    return _hostIntersectionIterator->indexInOutside();
  }

  GlobalCoords firstOuterNormal(const LocalCoords& local) const {
    return _hostIntersectionIterator->outerNormal(local);
  }

  GlobalCoords firstIntegrationOuterNormal(const LocalCoords& local) const {
    return _hostIntersectionIterator->integrationOuterNormal(local);
  }

  GlobalCoords firstUnitOuterNormal(const LocalCoords& local) const {
    return _hostIntersectionIterator->unitOuterNormal(local);
  }

  GlobalCoords secondOuterNormal(const LocalCoords& local) const {
    return -_hostIntersectionIterator->outerNormal(local);
  }

  GlobalCoords secondIntegrationOuterNormal(const LocalCoords& local) const {
    return -_hostIntersectionIterator->integrationOuterNormal(local);
  }

  GlobalCoords secondUnitOuterNormal(const LocalCoords& local) const {
    return -_hostIntersectionIterator->unitOuterNormal(local);
  }

private:

  GridView _gridView;
  HostGridView _hostGridView;

  SubDomainType _domain1;
  SubDomainType _domain2;

  HostIterator _hostIterator;
  HostIterator _hostEnd;

  HostIntersectionIterator _hostIntersectionIterator;
  HostIntersectionIterator _hostIntersectionEnd;

  MakeableGeometryWrapper<LocalGeometry::mydimension,LocalGeometry::coorddimension,GridImp> _geometryInInside, _geometryInOutside;
  MakeableGeometryWrapper<Geometry::mydimension,Geometry::coorddimension,GridImp> _geometry;

};

template<typename GridImp>
class LeafSubDomainInterfaceIterator :
    public SubDomainInterfaceIterator<GridImp,
                                      LeafSubDomainInterfaceIterator<GridImp>,
                                      typename GridImp::LeafGridView,
                                      typename detail::HostGridAccessor<GridImp>::Type::LeafGridView,
                                      LeafSubDomainInterfaceIterator<GridImp> >
{

  template<typename, typename, typename, typename, typename>
  friend class SubDomainInterfaceIterator;

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename>
  friend class MultiDomainGrid;

  typedef SubDomainInterfaceIterator<GridImp,
                                     LeafSubDomainInterfaceIterator<GridImp>,
                                     typename GridImp::LeafGridView,
                                     typename GridImp::HostGridType::LeafGridView,
                                     LeafSubDomainInterfaceIterator<GridImp> > Base;

  typedef LeafSubDomainInterfaceIterator<GridImp> Intersection;
  typedef typename GridImp::SubDomainType SubDomainType;

  LeafSubDomainInterfaceIterator(const GridImp& grid, SubDomainType subDomain1, SubDomainType subDomain2, bool end=false) :
    Base(grid.leafView(),grid._hostGrid.leafView(),subDomain1,subDomain2,end)
  {}

};


template<typename GridImp>
class LevelSubDomainInterfaceIterator :
    public SubDomainInterfaceIterator<GridImp,
                                      LevelSubDomainInterfaceIterator<GridImp>,
                                      typename GridImp::LevelGridView,
                                      typename detail::HostGridAccessor<GridImp>::Type::LevelGridView,
                                      LevelSubDomainInterfaceIterator<GridImp> >
{

  template<typename, typename, typename, typename, typename>
  friend class SubDomainInterfaceIterator;

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename>
  friend class MultiDomainGrid;

  typedef SubDomainInterfaceIterator<GridImp,
                                     LevelSubDomainInterfaceIterator<GridImp>,
                                     typename GridImp::LevelGridView,
                                     typename GridImp::HostGridType::LevelGridView,
                                     LevelSubDomainInterfaceIterator<GridImp> > Base;

  typedef LevelSubDomainInterfaceIterator<GridImp> Intersection;
  typedef typename GridImp::SubDomainType SubDomainType;

  LevelSubDomainInterfaceIterator(const GridImp& grid, SubDomainType subDomain1, SubDomainType subDomain2, int level, bool end=false) :
    Base(grid.levelView(level),grid._hostGrid.levelView(level),subDomain1,subDomain2,end)
  {}

};


} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAININTERFACEITERATOR_HH
