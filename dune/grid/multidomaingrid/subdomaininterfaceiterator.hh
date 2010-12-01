#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAININTERFACEITERATOR_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAININTERFACEITERATOR_HH

#include <dune/common/iteratorfacades.hh>

namespace Dune {

namespace mdgrid {

template<typename SubDomainIndexType>
class SubDomainToSubDomainController;

template<typename GridImp,
         typename WrapperImp,
         typename GridView,
         typename HostGridView,
         typename IntersectionType,
         typename IterationController>
class SubDomainInterfaceIterator : public ForwardIteratorFacade<SubDomainInterfaceIterator<
                                                                  GridImp,
                                                                  WrapperImp,
                                                                  GridView,
                                                                  HostGridView,
                                                                  IntersectionType,
                                                                  IterationController
                                                                  >,
                                                                IntersectionType>
{

  typedef typename HostGridView::template Codim<0>::Iterator HostIterator;
  typedef typename HostGridView::IntersectionIterator HostIntersectionIterator;
  typedef typename GridView::IntersectionIterator MultiDomainIntersectionIterator;
  typedef IntersectionType Intersection;
  typedef typename GridImp::SubDomainIndexType SubDomainIndexType;
  typedef SubDomainIndexType SubDomainType DUNE_DEPRECATED;

  typedef typename GridImp::Traits::template Codim<0>::EntityPointer EntityPointer;
  typedef typename GridImp::Traits::template Codim<0>::Entity Entity;
  typedef typename GridImp::Traits::template Codim<1>::Geometry Geometry;
  typedef typename GridImp::Traits::template Codim<1>::LocalGeometry LocalGeometry;

  typedef typename GridImp::SubDomainGrid::Traits::template Codim<0>::Entity SubDomainEntity;

  typedef typename GridImp::ctype ctype;
  static const int dimension = GridImp::dimension;
  static const int dimensionworld = GridImp::dimensionworld;

  typedef FieldVector<ctype,dimensionworld> GlobalCoords;
  typedef FieldVector<ctype,dimension - 1> LocalCoords;

  template<typename, typename, typename, typename>
  friend class ForwardIteratorFacade;

  template<typename>
  friend class SubDomainToSubDomainController;

protected:

  SubDomainInterfaceIterator(const GridView& gridView,
                             const HostGridView& hostGridView,
                             IterationController controller,
                             bool end) :
    _gridView(gridView),
    _hostGridView(hostGridView),
    _controller(controller),
    _hostIterator(end ? hostGridView.template end<0>() : hostGridView.template begin<0>()),
    _hostEnd(hostGridView.template end<0>()),
    _hostIntersectionIterator(hostGridView.ibegin(*hostGridView.template begin<0>())),
    _hostIntersectionEnd(hostGridView.iend(*hostGridView.template begin<0>())),
    _inverseHostIntersection(hostGridView.ibegin(*hostGridView.template begin<0>())),
    _inverseHostIntersectionValid(false)
  {
    _controller.incrementToStartPosition(*this);
  }

public:

  // The following two methods have to be public because of non-member comparison operators in the
  // IteratorFacade framework!

  bool equals(const WrapperImp& rhs) const {
    return (_hostIterator == rhs._hostIterator &&
            (_hostIterator == _hostEnd ||
             (_hostIntersectionIterator == rhs._hostIntersectionIterator &&
             _controller.subDomain1() == rhs._controller.subDomain1() &&
             _controller.subDomain2() == rhs._controller.subDomain2()
              )
             )
            );
  }

  bool equals(const SubDomainInterfaceIterator& rhs) const {
    return (_hostIterator == rhs._hostIterator &&
            (_hostIterator == _hostEnd ||
             (_hostIntersectionIterator == rhs._hostIntersectionIterator &&
             _controller.subDomain1() == rhs._controller.subDomain1() &&
             _controller.subDomain2() == rhs._controller.subDomain2()
              )
             )
            );
  }


private:

  void increment() {
    _controller.increment(*this);
    _inverseHostIntersectionValid = false;
    _geometry.clear();
    _geometryInInside.clear();
    _geometryInOutside.clear();
  }

  void findInverseHostIntersection() const {
    assert(_hostIntersectionIterator->neighbor());
    _inverseHostIntersection = _hostGridView.ibegin(*_hostIntersectionIterator->outside());
    while (_hostIntersectionIterator->inside() != _inverseHostIntersection->outside()) {
      ++_inverseHostIntersection;
    }
    _inverseHostIntersectionValid = true;
  }

  HostIntersectionIterator secondHostIntersectionIterator() const {
    if (!_inverseHostIntersectionValid)
      findInverseHostIntersection();
    return _inverseHostIntersection;
  }

  HostIntersectionIterator firstHostIntersectionIterator() const {
    return _hostIntersectionIterator;
  }

public:

  MultiDomainIntersectionIterator secondMultiDomainIntersectionIterator() const {
    return _gridView.grid().template multiDomainIntersectionIterator<GridView,HostGridView>(secondHostIntersectionIterator());
  }

  MultiDomainIntersectionIterator firstMultiDomainIntersectionIterator() const {
    return _gridView.grid().template multiDomainIntersectionIterator<GridView,HostGridView>(firstHostIntersectionIterator());
  }

private:

  const Intersection& dereference() const {
    return reinterpret_cast<const Intersection&>(*this);
  }


public:

  const Intersection& operator*() const {
    return dereference();
  }

  const Intersection* operator->() const {
    return &(dereference());
  }

  //! Returns an EntityPointer to the corresponding cell in the first subdomain.
  EntityPointer firstCell() const {
    return EntityPointerWrapper<0,GridImp>(_hostIntersectionIterator->inside());
  }

  //! Returns an EntityPointer to the corresponding cell in the second subdomain.
  EntityPointer secondCell() const {
    return EntityPointerWrapper<0,GridImp>(_hostIntersectionIterator->outside());
  }

  //! Returns true if this intersection is conforming.
  bool conforming() const {
    return _hostIntersectionIterator->conforming();
  }

  //! Returns the local geometry in the corresponding cell of the first subdomain.
  const LocalGeometry& geometryInFirstCell() const {
    if (!_geometryInInside.isSet()) {
      _geometryInInside.reset(_hostIntersectionIterator->geometryInInside());
    }
    return _geometryInInside;
  }

  //! Returns the local geometry in the corresponding cell of the second subdomain.
  const LocalGeometry& geometryInSecondCell() const {
    if (!_geometryInOutside.isSet()) {
      _geometryInOutside.reset(_hostIntersectionIterator->geometryInOutside());
    }
    return _geometryInOutside;
  }

  //! Returns the global geometry of the intersection.
  const Geometry& geometry() const {
    if (!_geometry.isSet()) {
      _geometry.reset(_hostIntersectionIterator->geometry());
    }
    return _geometry;
  }

  //! Returns the GeometryType of this intersection.
  GeometryType type() const {
    return _hostIntersectionIterator->type();
  }

  //! Returns the subindex of the corresponding face in the cell belonging to the
  //! first subdomain.
  int indexInFirstCell() const {
    return _hostIntersectionIterator->indexInInside();
  }

  //! Returns the subindex of the corresponding face in the cell belonging to the
  //! second subdomain.
  int indexInSecondCell() const {
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

  SubDomainIndexType subDomain1() const {
    return _controller.subDomain1();
  }

  SubDomainIndexType subDomain2() const {
    return _controller.subDomain2();
  }

  SubDomainIndexType domain1() const DUNE_DEPRECATED {
    return this->subDomain1();
  }

  SubDomainIndexType domain2() const DUNE_DEPRECATED {
    return this->subDomain2();
  }

private:

  GridView _gridView;
  HostGridView _hostGridView;

  IterationController _controller;

  HostIterator _hostIterator;
  HostIterator _hostEnd;

  HostIntersectionIterator _hostIntersectionIterator;
  HostIntersectionIterator _hostIntersectionEnd;

  mutable HostIntersectionIterator _inverseHostIntersection;
  mutable bool _inverseHostIntersectionValid;

  MakeableGeometryWrapper<LocalGeometry::mydimension,LocalGeometry::coorddimension,GridImp> _geometryInInside, _geometryInOutside;
  MakeableGeometryWrapper<Geometry::mydimension,Geometry::coorddimension,GridImp> _geometry;

};


template<typename SubDomainIndexType>
class SubDomainToSubDomainController
{

  template<typename GridImp,
           typename WrapperImp,
           typename GridView,
           typename HostGridView,
           typename IntersectionType,
           typename IntersectionController
           >
  friend class SubDomainInterfaceIterator;

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
                                      LeafSubDomainInterfaceIterator<GridImp>,
                                      typename GridImp::LeafGridView,
                                      typename detail::HostGridAccessor<GridImp>::Type::LeafGridView,
                                      LeafSubDomainInterfaceIterator<GridImp>,
                                      SubDomainToSubDomainController<typename GridImp::SubDomainIndexType>
                                      >
{

  template<typename, typename, typename, typename, typename, typename>
  friend class SubDomainInterfaceIterator;

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename,typename>
  friend class MultiDomainGrid;

  typedef SubDomainToSubDomainController<typename GridImp::SubDomainIndexType> Controller;

  typedef SubDomainInterfaceIterator<GridImp,
                                     LeafSubDomainInterfaceIterator<GridImp>,
                                     typename GridImp::LeafGridView,
                                     typename GridImp::HostGridType::LeafGridView,
                                     LeafSubDomainInterfaceIterator<GridImp>,
                                     Controller
                                     > Base;

  typedef LeafSubDomainInterfaceIterator<GridImp> Intersection;
  typedef typename GridImp::SubDomainIndexType SubDomainIndexType;
  typedef SubDomainIndexType SubDomainType DUNE_DEPRECATED;

  LeafSubDomainInterfaceIterator(const GridImp& grid, SubDomainIndexType subDomain1, SubDomainIndexType subDomain2, bool end=false) :
    Base(grid.leafView(),grid._hostGrid.leafView(),Controller(subDomain1,subDomain2),end)
  {}

};


template<typename GridImp>
class LevelSubDomainInterfaceIterator :
    public SubDomainInterfaceIterator<GridImp,
                                      LevelSubDomainInterfaceIterator<GridImp>,
                                      typename GridImp::LevelGridView,
                                      typename detail::HostGridAccessor<GridImp>::Type::LevelGridView,
                                      LevelSubDomainInterfaceIterator<GridImp>,
                                      SubDomainToSubDomainController<typename GridImp::SubDomainIndexType>
                                      >
{

  template<typename, typename, typename, typename, typename, typename>
  friend class SubDomainInterfaceIterator;

  template<int, int, typename>
  friend class EntityWrapper;

  template<typename,typename>
  friend class MultiDomainGrid;

  typedef SubDomainToSubDomainController<typename GridImp::SubDomainIndexType> Controller;

  typedef SubDomainInterfaceIterator<GridImp,
                                     LevelSubDomainInterfaceIterator<GridImp>,
                                     typename GridImp::LevelGridView,
                                     typename GridImp::HostGridType::LevelGridView,
                                     LevelSubDomainInterfaceIterator<GridImp>,
                                     Controller
                                     > Base;

  typedef LevelSubDomainInterfaceIterator<GridImp> Intersection;
  typedef typename GridImp::SubDomainIndexType SubDomainIndexType;
  typedef SubDomainIndexType SubDomainType DUNE_DEPRECATED;

  LevelSubDomainInterfaceIterator(const GridImp& grid, SubDomainIndexType subDomain1, SubDomainIndexType subDomain2, int level, bool end=false) :
    Base(grid.levelView(level),grid._hostGrid.levelView(level),Controller(subDomain1,subDomain2),end)
  {}

};


} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAININTERFACEITERATOR_HH
