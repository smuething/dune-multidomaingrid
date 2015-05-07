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

  using BaseT = DefaultLevelGridView<GridImp,pitype>;

public:

  using typename BaseT::IntersectionIterator;

  LevelGridView(const GridImp& grid, int level)
    : BaseT(grid,level)
  {}

  IntersectionIterator ibegin(const typename BaseT::template Codim<0>::Entity& entity) const
  {
    return IntersectionIteratorWrapper<
      const GridImp,
      typename BaseT::IndexSet,
      typename GridImp::MultiDomainGrid::LevelGridView::IntersectionIterator
      >(
        &this->indexSet(),
        this->grid()._grid.levelGridView(entity.level()).ibegin(GridImp::getRealImplementation(entity).multiDomainEntity())
        );
  }

  IntersectionIterator iend(const typename BaseT::template Codim<0>::Entity& entity) const
  {
    return IntersectionIteratorWrapper<
      const GridImp,
      typename BaseT::IndexSet,
      typename GridImp::MultiDomainGrid::LevelGridView::IntersectionIterator
      >(
        &this->indexSet(),
        this->grid()._grid.levelGridView(entity.level()).iend(GridImp::getRealImplementation(entity).multiDomainEntity())
        );
  }

};

template<typename GridImp, PartitionIteratorType pitype>
struct LevelGridViewTraits
  : public DefaultLevelGridViewTraits<GridImp,pitype>
{
  using GridViewImp = LevelGridView<GridImp,pitype>;
};



template<typename GridImp, PartitionIteratorType pitype>
class LeafGridView
  : public DefaultLeafGridView<GridImp,pitype>
{

  using BaseT = DefaultLeafGridView<GridImp,pitype>;

public:

  using typename BaseT::IntersectionIterator;

  LeafGridView(const GridImp& grid)
    : BaseT(grid)
  {}

  IntersectionIterator ibegin(const typename BaseT::template Codim<0>::Entity& entity) const
  {
    return IntersectionIteratorWrapper<
      const GridImp,
      typename BaseT::IndexSet,
      typename GridImp::MultiDomainGrid::LeafGridView::IntersectionIterator
      >(
        &this->indexSet(),
        this->grid()._grid.leafGridView().ibegin(GridImp::getRealImplementation(entity).multiDomainEntity())
        );
  }

  IntersectionIterator iend(const typename BaseT::template Codim<0>::Entity& entity) const
  {
    return IntersectionIteratorWrapper<
      const GridImp,
      typename BaseT::IndexSet,
      typename GridImp::MultiDomainGrid::LeafGridView::IntersectionIterator
      >(
        &this->indexSet(),
        this->grid()._grid.leafGridView().iend(GridImp::getRealImplementation(entity).multiDomainEntity())
        );
  }

};

template<typename GridImp, PartitionIteratorType pitype>
struct LeafGridViewTraits
  : public DefaultLeafGridViewTraits<GridImp,pitype>
{
  using GridViewImp = LeafGridView<GridImp,pitype>;
};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_GRIDVIEW_HH
