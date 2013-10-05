#ifndef DUNE_MULTIDOMAINGRID_HOSTGRIDACCESSOR_HH
#define DUNE_MULTIDOMAINGRID_HOSTGRIDACCESSOR_HH

namespace Dune {
namespace mdgrid {
namespace detail {

  template<typename GridImp>
  struct HostGridAccessor {
    typedef typename GridImp::HostGridType::Traits Traits;
    typedef typename GridImp::HostGridType Type;
  };

} // end namespace detail
} // end namespace mdgrid
} // end namespace Dune

#endif // DUNE_MULTIDOMAINGRID_HOSTGRIDACCESSOR_HH
