#ifndef DUNE_MULTIDOMAINGRID_HH
#define DUNE_MULTIDOMAINGRID_HH

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
  struct IsUnstructured<MultiDomainGrid<HostGrid> >
  {
    static const bool v = IsUnstructured<HostGrid>::v;
  };


  template<typename HostGrid>
  struct isAdaptable<MultiDomainGrid<HostGrid> >
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
