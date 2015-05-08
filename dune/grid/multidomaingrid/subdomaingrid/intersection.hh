#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_INTERSECTION_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_INTERSECTION_HH

#include <dune/grid/common/intersection.hh>

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<int codim, int dim, typename GridImp>
class EntityWrapper;

template<typename GridImp,
         typename IndexSet,
         typename MultiDomainIntersection_
         >
class IntersectionWrapper {

  template<class, class, class>
  friend class Dune::IntersectionIterator;

  template<class,class,class>
  friend class IntersectionIteratorWrapper;

  template<class, class>
  friend class Dune::Intersection;

  template<typename MDGrid>
  friend class SubDomainGrid;

  using MultiDomainIntersection = MultiDomainIntersection_;
  using EntityPointer           = typename GridImp::Traits::template Codim<0>::EntityPointer;
  using Entity                  = typename GridImp::Traits::template Codim<0>::Entity;
  using Geometry                = typename GridImp::Traits::template Codim<1>::Geometry;
  using LocalGeometry           = typename GridImp::Traits::template Codim<1>::LocalGeometry;
  using EntityWrapper           = Dune::mdgrid::subdomain::EntityWrapper<0,GridImp::dimension,GridImp>;

  static const int dimension = GridImp::dimension;
  static const int dimensionworld = GridImp::dimensionworld;

  using ctype            = typename GridImp::ctype;
  using GlobalCoords     = FieldVector<ctype,dimensionworld>;
  using LocalCoords      = FieldVector<ctype,dimension - 1>;

  IntersectionWrapper()
    : _indexSet(nullptr)
    , _intersectionTypeTested(false)
  {}

  IntersectionWrapper(const IndexSet* indexSet, const MultiDomainIntersection& multiDomainIntersection)
    : _indexSet(indexSet)
    , _multiDomainIntersection(multiDomainIntersection)
    , _intersectionTypeTested(false)
  {}

private:

  const typename GridImp
    ::MultiDomainGrid
    ::template ReturnImplementationType<MultiDomainIntersection>
    ::ImplementationType
    ::HostIntersection&
  hostIntersection() const {
    return GridImp::MultiDomainGrid::getRealImplementation(_multiDomainIntersection).hostIntersection();
  }


  bool equals(const IntersectionWrapper& rhs) const {
    return _indexSet == rhs._indexSet && _multiDomainIntersection == rhs._multiDomainIntersection;
  }

  void checkIntersectionType() const {
    if (!_intersectionTypeTested) {
      if (_multiDomainIntersection.boundary()) {
        _intersectionType = GridImp::boundary;
        _intersectionTypeTested = true;
        return;
      }
      if (!_multiDomainIntersection.neighbor()) {
        _intersectionType = GridImp::processor;
        _intersectionTypeTested = true;
        return;
      }
      if (_indexSet->containsMultiDomainEntity(_multiDomainIntersection.outside())) {
        _intersectionType = GridImp::neighbor;
        _intersectionTypeTested = true;
        return;
      } else {
        _intersectionType = GridImp::foreign;
        _intersectionTypeTested = true;
        return;
      }
      assert(false && "Should never get here - invalid intersection type!");
    }
  }

  bool boundary() const {
    checkIntersectionType();
    return
      _intersectionType == GridImp::boundary ||
      _intersectionType == GridImp::foreign;
  }

  int boundaryId() const {
    return _multiDomainIntersection.boundaryId();
  }

  std::size_t boundarySegmentIndex() const {
    checkIntersectionType();
    // FIXME: We need to do something about subdomain boundaries in the interior of
    //        the MultiDomainGrid
    return _intersectionType == GridImp::boundary ? _multiDomainIntersection.boundarySegmentIndex() : 0;
  }

  bool neighbor() const {
    checkIntersectionType();
    return _intersectionType == GridImp::neighbor;
  }

  Entity inside() const {
    return {EntityWrapper(&_indexSet->_grid,_multiDomainIntersection.inside())};
  }

  Entity outside() const {
    checkIntersectionType();
    assert(_intersectionType == GridImp::neighbor);
    return {EntityWrapper(&_indexSet->_grid,_multiDomainIntersection.outside())};
  }

  bool conforming() const {
    return _multiDomainIntersection.conforming();
  }

  LocalGeometry geometryInInside() const {
    return LocalGeometry(hostIntersection().geometryInInside());
  }

  LocalGeometry geometryInOutside() const {
    checkIntersectionType();
    assert(_intersectionType == GridImp::neighbor);
    return LocalGeometry(hostIntersection().geometryInOutside());
  }

  Geometry geometry() const {
    return Geometry(hostIntersection().geometry());
  }

  GeometryType type() const {
    return _multiDomainIntersection.type();
  }

  int indexInInside() const {
    return _multiDomainIntersection.indexInInside();
  }

  int indexInOutside() const {
    checkIntersectionType();
    assert(_intersectionType == GridImp::neighbor);
    return _multiDomainIntersection.indexInOutside();
  }

  GlobalCoords outerNormal(const LocalCoords& local) const {
    return _multiDomainIntersection.outerNormal(local);
  }

  GlobalCoords integrationOuterNormal(const LocalCoords& local) const {
    return _multiDomainIntersection.integrationOuterNormal(local);
  }

  GlobalCoords unitOuterNormal(const LocalCoords& local) const {
    return _multiDomainIntersection.unitOuterNormal(local);
  }

  GlobalCoords centerUnitOuterNormal() const {
    return _multiDomainIntersection.centerUnitOuterNormal();
  }

  typename GridImp::IntersectionType intersectionType() const {
    checkIntersectionType();
    return _intersectionType;
  }

  const MultiDomainIntersection& multiDomainIntersection() const
  {
    return _multiDomainIntersection;
  }

private:

  const IndexSet* _indexSet;
  MultiDomainIntersection _multiDomainIntersection;
  mutable bool _intersectionTypeTested;
  mutable typename GridImp::IntersectionType _intersectionType;

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_INTERSECTION_HH
