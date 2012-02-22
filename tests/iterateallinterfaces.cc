#include "config.h"
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/multidomaingrid.hh>
#include <iostream>

template<typename ST, typename E>
void assembleLocalInterfaceTerm(ST s1, const E& e1, ST s2, const E& e2)
{
  typename E::Geometry geo1 = e1.geometry();
  std::cout << s1 << " -> " << s2 << " (" << geo1.center() << " -> " << e2.geometry().center() << ")" << std::endl;
}

template<typename MDGV>
void iterate(MDGV mdgv)
{
  typedef typename MDGV::template Codim<0>::Iterator Iterator;
  typedef typename MDGV::template Codim<0>::EntityPointer EP;
  typedef typename MDGV::IntersectionIterator IIterator;
  typedef typename MDGV::IndexSet::SubDomainIndexType SubDomainIndexType;
  const typename MDGV::IndexSet& is = mdgv.indexSet();

  for (Iterator it = mdgv.template begin<0>(); it != mdgv.template end<0>(); ++it)
    {
      // this assumes that subdomains form a partition of the host grid
      const SubDomainIndexType insideSubDomain = *is.subDomains(*it).begin();
      for (IIterator iit = mdgv.ibegin(*it); iit != mdgv.iend(*it); ++iit)
        {
          if (!iit->neighbor())
            continue;
          const EP outside = iit->outside();
          const SubDomainIndexType outsideSubDomain = *is.subDomains(*outside).begin();
          if (insideSubDomain != outsideSubDomain)
            assembleLocalInterfaceTerm(insideSubDomain,*it,outsideSubDomain,*outside);
        }
    }
}

template<typename Grid>
void iterate2(const Grid& grid)
{
  typedef typename Grid::LeafAllSubDomainInterfacesIterator Iterator;
  typedef typename Grid::template Codim<0>::EntityPointer EP;

  for (Iterator it = grid.leafAllSubDomainInterfacesBegin(); it != grid.leafAllSubDomainInterfacesEnd(); ++it)
    {
      const EP ep1 = it->firstCell();
      const EP ep2 = it->secondCell();
      std::cout << it->subDomain1() << " -> " << it->subDomain2() << " (" << ep1->geometry().center() << " -> " << ep2->geometry().center() << ") " << it->geometry().center() << std::endl;
    }
  for (auto it = grid.leafSubDomainInterfaceBegin(0,1); it != grid.leafSubDomainInterfaceEnd(0,1); ++it)
    {
      const EP ep1 = it->firstCell();
      const EP ep2 = it->secondCell();
      std::cout << it->subDomain1() << " -> " << it->subDomain2() << " (" << ep1->geometry().center() << " -> " << ep2->geometry().center() << ") " << it->geometry().center() << std::endl;
    }
}


int main(int argc, char** argv)
{
  if (argc < 2)
    {
      std::cerr << "Usage: " << argv[0] << " <refinement level>" << std::endl;
      exit(1);
    }

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
  typedef MDGrid::SubDomainIndexType SubDomainIndexType;

  MDGV mdgv = mdgrid.leafView();

  mdgrid.startSubDomainMarking();

  for(Iterator it = mdgv.begin<0>(); it != mdgv.end<0>(); ++it)
    {
      SubDomainIndexType subdomain = 0;
      if (it->geometry().center()[0] > 0.5)
        subdomain += 1;
      if (it->geometry().center()[1] > 0.5)
        subdomain += 2;
      mdgrid.addToSubDomain(subdomain,*it);
      if (it->geometry().center()[1] > 0.5)
        mdgrid.addToSubDomain(4,*it);
      if (it->geometry().center()[0] > 0.5)
        mdgrid.addToSubDomain(5,*it);
      if (it->geometry().center()[0] > 0.5 && it->geometry().center()[1] > 0.5)
        mdgrid.addToSubDomain(6,*it);
      if (it->geometry().center()[0] > 0.5 && it->geometry().center()[1] < 0.5)
        mdgrid.addToSubDomain(7,*it);
    }

  mdgrid.preUpdateSubDomains();
  mdgrid.updateSubDomains();
  mdgrid.postUpdateSubDomains();

  mdgrid.globalRefine(atoi(argv[1]));

  //iterate(mdgv);
  iterate2(mdgrid);
}
