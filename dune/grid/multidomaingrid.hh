#ifndef DUNE_MULTIDOMAINGRID_HH
#define DUNE_MULTIDOMAINGRID_HH

#include <dune/grid/multidomaingrid/multidomaingrid.hh>

namespace Dune {

using mdgrid::MultiDomainGrid;

namespace Capabilities {

  template<class HostGrid, int codim>
  struct hasEntity<MultiDomainGrid<HostGrid>, codim>
  {
    static const bool v = hasEntity<HostGrid,codim>::v;
  };
  

  template<class HostGrid>
  struct isParallel<MultiDomainGrid<HostGrid> >
  {
    static const bool v = isParallel<HostGrid>::v;
  };


  template<class HostGrid>
  struct hasHangingNodes<MultiDomainGrid<HostGrid> >
  {
    static const bool v = hasHangingNodes<HostGrid>::v;
  };


  template<class HostGrid>
  struct isLevelwiseConforming<MultiDomainGrid<HostGrid> >
  {
    static const bool v = isLevelwiseConforming<HostGrid>::v;
  };


  template<class HostGrid>
  struct isLeafwiseConforming<MultiDomainGrid<HostGrid> >
  {
    static const bool v = isLeafwiseConforming<HostGrid>::v;
  };


  template<class HostGrid>
  struct isUnstructured<MultiDomainGrid<HostGrid> >
  {
    static const bool v = isUnstructured<HostGrid>::v;
  };


} // namespace Capabilities

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_HH
