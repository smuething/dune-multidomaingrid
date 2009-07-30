#ifndef DUNE_MULTIDOMAINGRID_HH
#define DUNE_MULTIDOMAINGRID_HH

#include <dune/common/collectivecommunication.hh>
#include <dune/grid/multidomaingrid/multidomaingrid.hh>
#include <dune/grid/multidomaingrid/multidomainmcmgmapper.hh>

namespace Dune {

using mdgrid::MultiDomainGrid;
using mdgrid::MultiDomainMCMGMapper;


namespace Capabilities {

  template<typename Grid>
  struct isAdaptable
  {
    static const bool v = true;
  };

}


// MultiDomainGrid capabilities

namespace Capabilities {

  template<class HostGrid, typename MDGridTraits, int codim>
  struct hasEntity<MultiDomainGrid<HostGrid,MDGridTraits>, codim>
  {
    static const bool v = hasEntity<HostGrid,codim>::v;
  };


  template<class HostGrid, typename MDGridTraits>
  struct isParallel<MultiDomainGrid<HostGrid,MDGridTraits> >
  {
    static const bool v = isParallel<HostGrid>::v;
  };


  template<class HostGrid, typename MDGridTraits>
  struct hasHangingNodes<MultiDomainGrid<HostGrid,MDGridTraits> >
  {
    static const bool v = hasHangingNodes<HostGrid>::v;
  };


  template<class HostGrid, typename MDGridTraits>
  struct isLevelwiseConforming<MultiDomainGrid<HostGrid,MDGridTraits> >
  {
    static const bool v = isLevelwiseConforming<HostGrid>::v;
  };


  template<class HostGrid, typename MDGridTraits>
  struct isLeafwiseConforming<MultiDomainGrid<HostGrid,MDGridTraits> >
  {
    static const bool v = isLeafwiseConforming<HostGrid>::v;
  };


  template<class HostGrid, typename MDGridTraits>
  struct IsUnstructured<MultiDomainGrid<HostGrid,MDGridTraits> >
  {
    static const bool v = IsUnstructured<HostGrid>::v;
  };


  template<typename HostGrid, typename MDGridTraits>
  struct isAdaptable<MultiDomainGrid<HostGrid,MDGridTraits> >
  {
    static const bool v = isAdaptable<HostGrid>::v;
  };

} // namespace Capabilities



// SubDomainGrid capabilities

namespace Capabilities {

  template<class MDGrid, int codim>
  struct hasEntity<Dune::mdgrid::subdomain::SubDomainGrid<MDGrid>, codim>
  {
    static const bool v = hasEntity<MDGrid,codim>::v;
  };


  template<class MDGrid>
  struct isParallel<Dune::mdgrid::subdomain::SubDomainGrid<MDGrid> >
  {
    static const bool v = isParallel<MDGrid>::v;
  };


  template<class MDGrid>
  struct hasHangingNodes<Dune::mdgrid::subdomain::SubDomainGrid<MDGrid> >
  {
    static const bool v = hasHangingNodes<MDGrid>::v;
  };


  template<class MDGrid>
  struct isLevelwiseConforming<Dune::mdgrid::subdomain::SubDomainGrid<MDGrid> >
  {
    static const bool v = isLevelwiseConforming<MDGrid>::v;
  };


  template<class MDGrid>
  struct isLeafwiseConforming<Dune::mdgrid::subdomain::SubDomainGrid<MDGrid> >
  {
    static const bool v = isLeafwiseConforming<MDGrid>::v;
  };


  template<class MDGrid>
  struct IsUnstructured<Dune::mdgrid::subdomain::SubDomainGrid<MDGrid> >
  {
    static const bool v = IsUnstructured<MDGrid>::v;
  };


  template<typename MDGrid>
  struct isAdaptable<Dune::mdgrid::subdomain::SubDomainGrid<MDGrid> >
  {
    static const bool v = false;
  };


} // namespace Capabilities

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_HH
