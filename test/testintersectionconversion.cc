#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/multidomaingrid.hh>
#include <iostream>
#include <cassert>

int main(int argc, char** argv)
{

  Dune::FieldVector<double,2> L(1.0);
  Dune::FieldVector<int,2> s(2);
  Dune::FieldVector<bool,2> p(false);
  int overlap = 0;

  typedef Dune::YaspGrid<2> HostGrid;
  HostGrid hostgrid(L,s,p,overlap);

  //typedef Dune::MultiDomainGrid<HostGrid,Dune::mdgrid::FewSubDomainsTraits<2,8> > MDGrid;
  typedef Dune::MultiDomainGrid<HostGrid,Dune::mdgrid::ArrayBasedTraits<2,8,8> > MDGrid;
  MDGrid mdgrid(hostgrid,true);

  typedef MDGrid::LeafGridView MDGV;
  typedef MDGV::Codim<0>::Iterator Iterator;

  MDGV mdgv = mdgrid.leafGridView();

  mdgrid.startSubDomainMarking();

  for(Iterator it = mdgv.begin<0>(); it != mdgv.end<0>(); ++it)
      mdgrid.addToSubDomain(0,*it);

  mdgrid.preUpdateSubDomains();
  mdgrid.updateSubDomains();
  mdgrid.postUpdateSubDomains();

  typedef MDGrid::SubDomainGrid SDGrid;
  const SDGrid& sdgrid = mdgrid.subDomain(0);
  SDGrid::LeafGridView sdgv = sdgrid.leafGridView();

  SDGrid::LeafGridView::Codim<0>::Iterator it = sdgv.begin<0>();
  SDGrid::LeafGridView::IntersectionIterator iit = sdgv.ibegin(*it);

  const MDGrid::LeafGridView::Intersection& is1 = sdgrid.multiDomainIntersection(*iit);
  const MDGrid::LeafGridView::Intersection& is2 = mdgrid.multiDomainIntersection(*iit);

  assert(is1.geometry().center() == is2.geometry().center());
}
