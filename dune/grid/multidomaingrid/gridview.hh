#ifndef DUNE_MULTIDOMAINGRID_GRIDVIEW_HH
#define DUNE_MULTIDOMAINGRID_GRIDVIEW_HH

#include <dune/grid/common/gridview.hh>
#include <dune/grid/common/defaultgridview.hh>

namespace Dune {

namespace mdgrid {

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
    return LevelIntersectionIteratorWrapper<const GridImp>(this->grid()._hostGrid.levelGridView(entity.level()).ibegin(*GridImp::getRealImplementation(entity).hostEntityPointer()));
  }

  IntersectionIterator iend(const typename BaseT::template Codim<0>::Entity& entity) const
  {
    return LevelIntersectionIteratorWrapper<const GridImp>(this->grid()._hostGrid.levelGridView(entity.level()).iend(*GridImp::getRealImplementation(entity).hostEntityPointer()));
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

  IntersectionIterator ibegin(const typename BaseT:: template Codim<0>::Entity& entity) const
  {
    return LeafIntersectionIteratorWrapper<const GridImp>(this->grid()._hostGrid.leafGridView().ibegin(*GridImp::getRealImplementation(entity).hostEntityPointer()));
  }

  IntersectionIterator iend(const typename BaseT:: template Codim<0>::Entity& entity) const
  {
    return LeafIntersectionIteratorWrapper<const GridImp>(this->grid()._hostGrid.leafGridView().iend(*GridImp::getRealImplementation(entity).hostEntityPointer()));
  }

};

template<typename GridImp, PartitionIteratorType pitype>
struct LeafGridViewTraits
  : public DefaultLeafGridViewTraits<GridImp,pitype>
{
  typedef LeafGridView<GridImp,pitype> GridViewImp;
};



} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_GRIDVIEW_HH
