#include "config.h"

#include <dune/common/unused.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/multidomaingrid.hh>
#include "output.hh"


int main(int argc, char** argv) {
  try {
    Dune::MPIHelper::instance(argc,argv);
    typedef Dune::YaspGrid<2> GridType;
    Dune::FieldVector<double,2> L(1.0);
    Dune::array<int,2> N = { {1,1} };
    GridType wgrid(L,N);
    typedef Dune::MultiDomainGrid<GridType,Dune::mdgrid::FewSubDomainsTraits<GridType::dimension,4> > Grid;
    Grid grid(wgrid);
    grid.globalRefine(5);
    typedef Grid::LeafGridView GridView;
    GridView gv = grid.leafGridView();
    grid.startSubDomainMarking();
    for (const auto& cell : elements(gv)) {
      auto c = cell.geometry().global(
        Dune::ReferenceElements<GridType::ctype,2>::general(cell.type()).position(0,0));
      double x = c[0];
      double y = c[1];
      if (x > 0.2) {
        if (y > 0.3 && y < 0.7) {
          if (x < 0.8)
            grid.addToSubDomain(1,cell);
          else // if (x > 0.6)
            grid.addToSubDomain(0,cell);
        } else {
            grid.addToSubDomain(0,cell);
        }
      }
    }
    grid.preUpdateSubDomains();
    grid.updateSubDomains();
    grid.postUpdateSubDomains();

    Grid::SubDomainGrid::LeafGridView::Codim<0>::Entity se0(*grid.subDomain(0).leafGridView().begin<0>());
    const auto& e0 DUNE_UNUSED = grid.subDomain(0).multiDomainEntity(se0);

    printStatus(grid,"partitioning1");

    grid.startSubDomainMarking();
    for (const auto& cell : elements(gv)) {
      const Grid::MDGridTraits::Codim<0>::SubDomainSet& domains = gv.indexSet().subDomains(cell);
      if (domains.contains(0) && domains.contains(1)) {
        grid.removeFromSubDomain(0,cell);
        grid.removeFromSubDomain(1,cell);
      } else if (domains.contains(0)) {
        grid.assignToSubDomain(1,cell);
      } else if (domains.contains(1)) {
        grid.assignToSubDomain(0,cell);
      } else {
        grid.addToSubDomain(0,cell);
        grid.addToSubDomain(1,cell);
      }
    }

    grid.preUpdateSubDomains();
    grid.updateSubDomains();
    grid.postUpdateSubDomains();

    printStatus(grid,"partitioning2");

    /*for (int i = 0; i <= 2; ++i) {
      std::cout << "codim " << i << ":";
      for (int s = 0; s <= 2; ++s) {
        std::cout << " " << is.size(i,s);
      }
      std::cout << std::endl;
      }*/
    return 0;
  } catch (Dune::Exception& e) {
    std::cout << e << std::endl;
    return 1;
  } catch (...) {
    std::cout << "generic exception" << std::endl;
    return 2;
  }
}
