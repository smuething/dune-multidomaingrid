#include "config.h"
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/multidomaingrid.hh>
#include "output.hh"
#include <cstdlib>


int main(int argc, char** argv) {
  try {
    Dune::MPIHelper::instance(argc,argv);
    typedef Dune::YaspGrid<2> GridType;
    Dune::GridPtr<GridType> gridPtr("/Users/smuething/Documents/dune/ws/dune-grid-howto/grids/unitcube2.dgf");
    GridType& wgrid = *gridPtr;
    typedef Dune::MultiDomainGrid<GridType,Dune::mdgrid::ArrayBasedTraits<GridType::dimension,1,65536> > Grid;
    Grid grid(wgrid,false);
    //typedef Dune::mdgrid::DynamicSubDomainCountTraits<GridType::dimension,1> MDGridTraits;
    //typedef Dune::MultiDomainGrid<GridType,MDGridTraits> Grid;
    //MDGridTraits md_grid_traits(65536);
    //Grid grid(wgrid,md_grid_traits,false);
    grid.globalRefine(8);
    typedef Grid::LeafGridView GridView;
    GridView gv = grid.leafGridView();
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
    vtkOut(gv,"largedomainnumbers_leafGridView",grid.leafSubDomainInterfaceBegin(0,1),grid.leafSubDomainInterfaceEnd(0,1));

    //vtkOut2(grid.subDomain(0).leafGridView(),"subdomain0");
    //vtkOut2(grid.subDomain(1).leafGridView(),"subdomain1");

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
