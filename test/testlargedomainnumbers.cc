#include "config.h"
#include <dune/common/mpihelper.hh>
#include <dune/grid/multidomaingrid.hh>
#include "output.hh"
#include <cstdlib>


int main(int argc, char** argv) {
  try {
    Dune::MPIHelper::instance(argc,argv);
    typedef Dune::YaspGrid<2> GridType;
    Dune::GridPtr<GridType> gridPtr("/Users/muethisn/Documents/dune/ws/dune-grid-howto/grids/unitcube2.dgf");
    GridType& wgrid = *gridPtr;
    typedef Dune::MultiDomainGrid<GridType,Dune::mdgrid::ArrayBasedTraits<GridType::dimension,1,65536> > Grid;
    Grid grid(wgrid,false);
    grid.globalRefine(8);
    typedef Grid::LeafGridView GridView;
    GridView gv = grid.leafView();
    typedef GridView::Codim<0>::Iterator Iterator;
    typedef GridView::Codim<2>::Iterator VIterator;
    typedef GridView::Codim<0>::Entity Entity;
    typedef GridView::Codim<0>::Geometry Geometry;
    grid.startSubDomainMarking();
    GridView::IndexSet::IndexType sd=0;
    int c = 0;
    const int sdsize = atoi(argv[1]);
    for (Iterator it = gv.begin<0>(); it != gv.end<0>(); ++it, ++c) {
     if (c == sdsize) {
        ++sd;
        c = 0;
      }
      const Entity& e = *it;
      grid.addToSubDomain(sd,e);
    }
    std::cout << "Number of subdomains: " << sd + 1 << std::endl;
    grid.preUpdateSubDomains();
    grid.updateSubDomains();
    grid.postUpdateSubDomains();

    const Grid::SubDomainGrid& sd0 = grid.subDomain(0);
    vtkOut(gv,"largedomainnumbers_leafView",grid.leafSubDomainInterfaceBegin(0,1),grid.leafSubDomainInterfaceEnd(0,1));

    //vtkOut2(grid.subDomain(0).leafView(),"subdomain0");
    //vtkOut2(grid.subDomain(1).leafView(),"subdomain1");

    /*for (int i = 0; i <= 2; ++i) {
      std::cout << "codim " << i << ":";
      for (int s = 0; s <= 2; ++s) {
	std::cout << " " << is.size(i,s);
      }
      std::cout << std::endl;
      }*/
  } catch (Dune::Exception& e) {
    std::cout << e << std::endl;
  }
  return 0;
}