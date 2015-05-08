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
  typedef typename MDGV::template Codim<0>::Entity Entity;
  typedef typename MDGV::IndexSet::SubDomainIndex SubDomainIndex;
  const typename MDGV::IndexSet& is = mdgv.indexSet();

  for (const auto& cell : elements(mdgv))
    {
      // this assumes that subdomains form a partition of the host grid
      const SubDomainIndex insideSubDomain = *is.subDomains(cell).begin();
      for (const auto& intersection : intersections(mdgv,cell))
        {
          if (!intersection.neighbor())
            continue;
          const Entity outside = intersection.outside();
          const SubDomainIndex outsideSubDomain = intersection.subDomains(outside).begin();
          if (insideSubDomain != outsideSubDomain)
            assembleLocalInterfaceTerm(insideSubDomain,cell,outsideSubDomain,outside);
        }
    }
}

template<typename Grid>
void iterate2(const Grid& grid)
{
  typedef typename Grid::LeafAllSubDomainInterfacesIterator Iterator;
  typedef typename Grid::template Codim<0>::Entity Entity;

  for (Iterator it = grid.leafAllSubDomainInterfacesBegin(); it != grid.leafAllSubDomainInterfacesEnd(); ++it)
    {
      const Entity e1 = it->firstCell();
      const Entity e2 = it->secondCell();
      std::cout << it->subDomain1() << " -> " << it->subDomain2() << " (" << e1.geometry().center() << " -> " << e2.geometry().center() << ") " << it->geometry().center() << std::endl;
    }
  for (auto it = grid.leafSubDomainInterfaceBegin(0,1); it != grid.leafSubDomainInterfaceEnd(0,1); ++it)
    {
      const Entity e1 = it->firstCell();
      const Entity e2 = it->secondCell();
      std::cout << it->subDomain1() << " -> " << it->subDomain2() << " (" << e1.geometry().center() << " -> " << e2.geometry().center() << ") " << it->geometry().center() << std::endl;
    }
}


template<typename MDGrid>
void driver(MDGrid& mdgrid)
{
  typedef typename MDGrid::LeafGridView MDGV;
  typedef typename MDGrid::SubDomainIndex SubDomainIndex;

  MDGV mdgv = mdgrid.leafGridView();

  mdgrid.startSubDomainMarking();

  for(const auto& cell : elements(mdgv))
    {
      SubDomainIndex subdomain = 0;
      if (cell.geometry().center()[0] > 0.5)
        subdomain += 1;
      if (cell.geometry().center()[1] > 0.5)
        subdomain += 2;
      mdgrid.addToSubDomain(subdomain,cell);
      if (cell.geometry().center()[1] > 0.5)
        mdgrid.addToSubDomain(4,cell);
      if (cell.geometry().center()[0] > 0.5)
        mdgrid.addToSubDomain(5,cell);
      if (cell.geometry().center()[0] > 0.5 && cell.geometry().center()[1] > 0.5)
        mdgrid.addToSubDomain(6,cell);
      if (cell.geometry().center()[0] > 0.5 && cell.geometry().center()[1] < 0.5)
        mdgrid.addToSubDomain(7,cell);
    }

  mdgrid.preUpdateSubDomains();
  mdgrid.updateSubDomains();
  mdgrid.postUpdateSubDomains();

  //iterate(mdgv);
  iterate2(mdgrid);
}


int main(int argc, char** argv)
{

  try {

    Dune::MPIHelper::instance(argc,argv);

    int refinement = 5;

    if (argc < 2)
      {
        std::cerr << "Usage: " << argv[0] << " <refinement level>" << std::endl;
        std::cerr << std::endl << "Defaulting to refinement level = " << refinement << std::endl;
      }
    else
      refinement = atoi(argv[1]);

    Dune::FieldVector<double,2> L(1.0);
    Dune::array<int,2> s = {{2, 2}};

    typedef Dune::YaspGrid<2> HostGrid;
    HostGrid hostgrid(L,s);
    hostgrid.globalRefine(refinement);

    {
      typedef Dune::MultiDomainGrid<HostGrid,Dune::mdgrid::FewSubDomainsTraits<2,8> > MDGrid;
      MDGrid mdgrid(hostgrid,true);
      driver(mdgrid);
    }

    {
      typedef Dune::MultiDomainGrid<HostGrid,Dune::mdgrid::ArrayBasedTraits<2,8,8> > MDGrid;
      MDGrid mdgrid(hostgrid,true);
      driver(mdgrid);
    }

    {
      typedef Dune::mdgrid::DynamicSubDomainCountTraits<2,8> Traits;
      typedef Dune::MultiDomainGrid<HostGrid,Dune::mdgrid::DynamicSubDomainCountTraits<2,8> > MDGrid;
      MDGrid mdgrid(hostgrid,Traits(8),true);
      driver(mdgrid);
    }

    return 0;

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

}
