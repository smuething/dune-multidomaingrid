#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_GRIDVIEW_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_GRIDVIEW_HH

#include <dune/grid/common/gridview.hh>
#include <dune/grid/common/defaultgridview.hh>

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<typename>
class LevelIntersectionIteratorWrapper;

template<typename>
class LeafIntersectionIteratorWrapper;

template<typename GridImp, PartitionIteratorType pitype>
class LevelGridView
  : public DefaultLevelGridView<GridImp,pitype>
{

  typedef DefaultLevelGridView<GridImp,pitype> BaseT;

public:

  typedef typename BaseT::IntersectionIterator IntersectionIterator;

  LevelGridView(const GridImp& grid, int level)
    : BaseT(grid,level)
  {}

  IntersectionIterator ibegin(const typename BaseT::template Codim<0>::Entity& entity) const
  {
    return LevelIntersectionIteratorWrapper<const GridImp>(this->grid(),entity.level(),this->grid()._grid.levelView(entity.level()).ibegin(*GridImp::getRealImplementation(entity).multiDomainEntityPointer()));
  }

  IntersectionIterator iend(const typename BaseT::template Codim<0>::Entity& entity) const
  {
    return LevelIntersectionIteratorWrapper<const GridImp>(this->grid(),entity.level(),this->grid()._grid.levelView(entity.level()).iend(*GridImp::getRealImplementation(entity).multiDomainEntityPointer()));
  }

};

template<typename GridImp, PartitionIteratorType pitype>
struct LevelGridViewTraits
  : public DefaultLevelGridViewTraits<GridImp,pitype>
{
  typedef LevelGridView<GridImp,pitype> GridViewImp;
};



template<typename GridImp, PartitionIteratorType pitype>
class LeafGridView
  : public DefaultLeafGridView<GridImp,pitype>
{

  typedef DefaultLeafGridView<GridImp,pitype> BaseT;

public:

  typedef typename BaseT::IntersectionIterator IntersectionIterator;

  LeafGridView(const GridImp& grid)
    : BaseT(grid)
  {}

  IntersectionIterator ibegin(const typename BaseT::template Codim<0>::Entity& entity) const
  {
    return LeafIntersectionIteratorWrapper<const GridImp>(this->grid(),this->grid()._grid.leafView().ibegin(*GridImp::getRealImplementation(entity).multiDomainEntityPointer()));
  }

  IntersectionIterator iend(const typename BaseT::template Codim<0>::Entity& entity) const
  {
    return LeafIntersectionIteratorWrapper<const GridImp>(this->grid(),this->grid()._grid.leafView().iend(*GridImp::getRealImplementation(entity).multiDomainEntityPointer()));
  }

};

template<typename GridImp, PartitionIteratorType pitype>
struct LeafGridViewTraits
  : public DefaultLeafGridViewTraits<GridImp,pitype>
{
  typedef LeafGridView<GridImp,pitype> GridViewImp;
};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_GRIDVIEW_HH
