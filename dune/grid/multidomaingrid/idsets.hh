#ifndef DUNE_MULTIDOMAINGRID_IDSETS_HH
#define DUNE_MULTIDOMAINGRID_IDSETS_HH

namespace Dune {

namespace mdgrid {

template<typename GridImp, typename WrappedIdSet>
class IdSetWrapper :
    public Dune::IdSet<GridImp,IdSetWrapper<GridImp,WrappedIdSet>,
		       typename WrappedIdSet::IdType>
{

  typedef typename remove_const<GridImp>::type::HostGridType HostGrid;
  typedef typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity Codim0Entity;

public:

  typedef WrappedIdSet::IdType IdType;

  template<int codim>
  IdType id(const typename remove_const<GridImp>::type::Traits::template Codim<codim>::Entity& e) const {
    return _wrappedIdSet.id(_grid.hostEntity(e));
  }

  template<typename Entity>
  IdType id(const Entity& e) const {
    return _wrappedIdSet.id(_grid.hostEntity(e));
  }

  template<int codim>
  IdType subId(const Codim0Entity& e, int i) const {
    return _wrappedIdSet.subId(_grid.hostEntity(e),i,codim);
  }

  IdType subId(const Codim0Entity& e, int i, unsigned int codim) const {
    return _wrappedIdSet.subId(_grid.hostEntity(e),i,codim);
  }

private:

  const GridImp& _grid;
  const WrappedIdSet& _wrappedIdSet;

  IdSetWrapper(const GridImp& grid, const WrappedIdSet& wrappedIdSet) :
    _grid(grid),
    _wrappedIdSet(wrappedIdSet)
  {}

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_IDSETS_HH
