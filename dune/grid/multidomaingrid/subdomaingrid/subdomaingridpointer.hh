#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_SUBDOMAINGRIDPOINTER_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_SUBDOMAINGRIDPOINTER_HH

namespace Dune {

namespace mdgrid {

namespace subdomain {

template<typename>
class SubDomainGrid;

template<typename SubDomainGridType>
class SubDomainGridPointer {

  template<typename>
  friend class SubDomainGrid;

public:

  typedef SubDomainGridType SubDomainGrid;

  const SubDomainGrid& operator*() const {
    return _grid;
  }

  const SubDomainGrid* operator->() const {
    return &_grid;
  }

  SubDomainGridPointer(const SubDomainGridPointer& rhs) {
    _grid.reset(*rhs);
  }

  const SubDomainGridPointer& operator=(const SubDomainGridPointer& rhs) {
    _grid.reset(*rhs);
  }

private:
  SubDomainGrid _grid;

  SubDomainGridPointer(const typename SubDomainGridType::MDGridType& multiDomainGrid, typename SubDomainGridType::SubDomainType subDomain) :
    _grid(multiDomainGrid,subDomain)
  {}

};


} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_SUBDOMAINGRIDPOINTER_HH
