#ifndef DUNE_MULTIDOMAINGRID_INDEXSETS_HH
#define DUNE_MULTIDOMAINGRID_INDEXSETS_HH

#include <unordered_map>
#include <vector>
#include <array>
#include <algorithm>
#include <type_traits>
#include <tuple>
#include <boost/bind.hpp>
#include <dune/grid/multidomaingrid/utility.hh>

namespace Dune {

namespace mdgrid {

template<typename HostGrid>
class MultiDomainGrid;


template<typename GridImp, typename HostIndexSet>
class IndexSetWrapper :
    public Dune::IndexSet<GridImp,IndexSetWrapper<GridImp,HostIndexSet>,
		       typename HostIndexSet::IndexType>
{

  template<typename HostGrid>
  friend class MultiDomainGrid;

  typedef typename remove_const<GridImp>::type::HostGridType HostGrid;
  typedef typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity Codim0Entity;

public:

  typedef typename HostIndexSet::IndexType IndexType;
  static const int dimension = remove_const<GridImp>::type::dimension;

  template<int codim>
  IndexType index(const typename remove_const<GridImp>::type::Traits::template Codim<codim>::Entity& e) const {
    return _hostIndexSet->index(_grid.hostEntity(e));
  }

  template<typename Entity>
  IndexType index(const Entity& e) const {
    return _hostIndexSet->index(_grid.hostEntity(e));
  }

  template<int codim>
  IndexType subIndex(const Codim0Entity& e, int i) const {
    return _hostIndexSet->subIndex(_grid.hostEntity(e),i,codim);
  }

  IndexType subIndex(const Codim0Entity& e, int i, unsigned int codim) const {
    return _hostIndexSet->subIndex(_grid.hostEntity(e),i,codim);
  }

  const std::vector<GeometryType>& geomTypes(int codim) const {
    return _hostIndexSet->geomTypes(codim);
  }

  IndexType size(GeometryType type) const {
    return _hostIndexSet->size(type);
  }

  IndexType size(int codim) const {
    return _hostIndexSet->size(codim);
  }

  template<typename EntityType>
  bool contains(const EntityType& e) const {
    return _hostIndexSet->contains(_grid.hostEntity(e));
  }

private:

  const GridImp& _grid;
  const HostIndexSet* _hostIndexSet;

  IndexSetWrapper(const GridImp& grid) :
    _grid(grid),
    _hostIndexSet(NULL)
  {}

  void update(const HostIndexSet& indexSet) {
    _hostIndexSet = &indexSet;
  }

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_INDEXSETS_HH
