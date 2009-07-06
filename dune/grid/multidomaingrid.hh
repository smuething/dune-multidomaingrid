#ifndef DUNE_MULTIDOMAINGRID_MULTIDOMAINGRID_HH
#define DUNE_MULTIDOMAINGRID_MULTIDOMAINGRID_HH

#include <vector>
#include <string>

#include <dune/common/collectivecommunication.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>

#include <dune/common/scopedptr.hh>

namespace Dune {

template <typename HostGrid>
class MultiDomainGrid;

template<int dim, typename HostGrid>
struct MultiDomainGridFamily
{
  typedef GridTraits<
    dim,
    HostGrid::dimensionworld,
    MultiDomainGrid<HostGrid>,
    MultiDomainGridGeometry,
    MultiDomainGridEntity,
    MultiDomainGridEntityPointer,
    MultiDomainGridLevelIterator,
    MultiDomainGridLeafIntersection,
    MultiDomainGridLevelIntersection,
    MultiDomainGridLeafIntersectionIterator,
    MultiDomainGridLevelIntersectionIterator,
    MultiDomainGridHierarchicIterator,
    MultiDomainGridLeafIterator,
    MultiDomainGridLevelIndexSet<const MultiDomainGrid<HostGrid> >,
    MultiDomainGridLeafIndexSet<const MultiDomainGrid<HostGrid> >,
    MultiDomainGridGlobalIdSet<const MultiDomainGrid<HostGrid> >,
    typename HostGrid::Traits::GlobalIdSet::IdType,
    MultiDomainGridLocalIdSet<const MultiDomainGrid<HostGrid> >,
    typename HostGrid::Traits::LocalIdSet::IdType,
    CollectiveCommunication<MultiDomainGrid<HostGrid> >
    > Traits;
};

template<typename HostGrid>
class MultiDomainGrid :
    public GridDefaultImplementation<HostGrid::dimension,
				     HostGrid::dimensionworld,
				     HostGrid::ctype,
				     MultiDomainGridFamily<
				       HostGrid::dimension,
				       HostGrid
				       >
				     > {

public:

  typedef MultiDomainGridFamily<HostGrid::dimension,HostGrid> GridFamily;

  typedef typename GridFamily::Traits Traits;

  typedef typename HostGrid::ctype ctype;

  explicit MultiDomainGrid(HostGrid& hostgrid) :
    _hostgrid(hostgrid),
    _leafIndexSet(*this),
    _globalIdSet(*this),
    _localIdSet(*this)
  {

  }

private:

};

} // namespace multidomaingrid

} // namespace Dune


#endif // DUNE_MULTIDOMAINGRID_MULTIDOMAINGRID_HH
