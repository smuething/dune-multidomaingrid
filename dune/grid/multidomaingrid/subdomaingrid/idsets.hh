#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_IDSETS_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_IDSETS_HH

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<typename MDGrid>
class SubDomainGrid;

template<typename GridImp, typename HostIdSet>
class IdSetWrapper :
    public Dune::IdSet<GridImp,IdSetWrapper<GridImp,HostIdSet>,
		       typename HostIdSet::IdType>
{

  template<typename>
  friend class SubDomainGrid;

  typedef typename remove_const<GridImp>::type::HostGridType HostGrid;
  typedef typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity Codim0Entity;

public:

  typedef typename HostIdSet::IdType IdType;

  template<int codim>
  IdType id(const typename remove_const<GridImp>::type::Traits::template Codim<codim>::Entity& e) const {
    assert(_grid.containsHostEntity(_grid.hostEntity(e)));
    return _hostIdSet.id(_grid.hostEntity(e));
  }

  template<typename Entity>
  IdType id(const Entity& e) const {
    assert(_grid.containsHostEntity(_grid.hostEntity(e)));
    return _hostIdSet.id(_grid.hostEntity(e));
  }

  template<int codim>
  IdType subId(const Codim0Entity& e, int i) const {
    assert(_grid.containsHostEntity(_grid.hostEntity(e)));
    return _hostIdSet.subId(_grid.hostEntity(e),i,codim);
  }

  IdType subId(const Codim0Entity& e, int i, unsigned int codim) const {
    assert(_grid.containsHostEntity(_grid.hostEntity(e)));
    return _hostIdSet.subId(_grid.hostEntity(e),i,codim);
  }

private:

  const GridImp& _grid;
  const HostIdSet& _hostIdSet;

  IdSetWrapper(const GridImp& grid, const HostIdSet& hostIdSet) :
    _grid(grid),
    _hostIdSet(hostIdSet)
  {}

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_IDSETS_HH
