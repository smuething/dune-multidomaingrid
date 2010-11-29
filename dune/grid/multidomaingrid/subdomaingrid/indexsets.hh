#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_INDEXSETS_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_INDEXSETS_HH

#include <unordered_map>
#include <vector>
#include <array>
#include <algorithm>
#include <type_traits>
#include <tuple>
#include <boost/scoped_ptr.hpp>
#include <boost/bind.hpp>
#include <dune/grid/multidomaingrid/utility.hh>

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<typename MDGrid>
class SubDomainGrid;


template<typename GridImp, typename MDIndexSet>
class IndexSetWrapper :
    public Dune::IndexSet<GridImp,IndexSetWrapper<GridImp,MDIndexSet>,
                          typename MDIndexSet::IndexType>
{

  template<typename, typename>
  friend class IndexSetWrapper;

  template<typename MDGrid>
  friend class SubDomainGrid;

  template<int, PartitionIteratorType, typename>
  friend class LeafIteratorWrapper;

  template<int, PartitionIteratorType, typename>
  friend class LevelIteratorWrapper;

  template<typename, typename, typename, typename>
  friend class IntersectionWrapper;

  typedef IndexSetWrapper<GridImp,MDIndexSet> ThisType;

  typedef typename remove_const<GridImp>::type::HostGridType HostGrid;
  typedef typename remove_const<GridImp>::type::MDGridType MDGrid;

  typedef typename remove_const<GridImp>::type::ctype ctype;

public:

  //typedef typename remove_const<GridImp>::type::SubDomainSet SubDomainSet;
  typedef typename remove_const<GridImp>::type::SubDomainIndexType SubDomainIndexType;
  typedef SubDomainIndexType SubDomainType DUNE_DEPRECATED;
  typedef typename MDIndexSet::IndexType IndexType;
  static const int dimension = remove_const<GridImp>::type::dimension;
  //typedef typename SubDomainSet::DomainType DomainType;
  //static const std::size_t maxSubDomains = MDGrid::MDGridTraits::template Codim

private:

  typedef typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity Codim0Entity;

public:

  template<int codim>
  IndexType index(const typename remove_const<GridImp>::type::Traits::template Codim<codim>::Entity& e) const {
    return _mdIndexSet.template index<codim>(_grid.domain(),_grid.multiDomainEntity(e));
  }

  template<typename Entity>
  IndexType index(const Entity& e) const {
    return index<Entity::codimension>(e);
  }

  template<int codim>
  IndexType subIndex(const Codim0Entity& e, int i) const {
    return _mdIndexSet.subIndex(_grid.domain(),_grid.multiDomainEntity(e),i,codim);
  }

  IndexType subIndex(const Codim0Entity& e, int i, unsigned int codim) const {
    return _mdIndexSet.subIndex(_grid.domain(),_grid.multiDomainEntity(e),i,codim);
  }

  const std::vector<GeometryType>& geomTypes(int codim) const {
    return _mdIndexSet.geomTypesForSubDomain(_grid.domain(),codim);
  }

  IndexType size(GeometryType type) const {
    return _mdIndexSet.sizeForSubDomain(_grid.domain(),type);
  }

  IndexType size(int codim) const {
    return _mdIndexSet.sizeForSubDomain(_grid.domain(),codim);
  }

  template<typename EntityType>
  bool contains(const EntityType& e) const {
    return _mdIndexSet.contains(_grid.domain(),_grid.multiDomainEntity(e));
  }

  bool operator==(const IndexSetWrapper& rhs) const {
    return (_grid == rhs._grid && &_mdIndexSet == &rhs._mdIndexSet);
  }

private:

  const GridImp& _grid;
  const MDIndexSet& _mdIndexSet;

  IndexSetWrapper(const GridImp& grid, const MDIndexSet& mdIndexSet) :
    _grid(grid),
    _mdIndexSet(mdIndexSet)
  {}

  template<typename EntityType>
  bool containsHostEntity(const EntityType& he) const {
    return _mdIndexSet.containsForSubDomain(_grid.domain(),he);
  }

  template<typename EntityType>
  bool containsMultiDomainEntity(const EntityType& mde) const {
    return _mdIndexSet.contains(_grid.domain(),mde);
  }

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_INDEXSETS_HH
