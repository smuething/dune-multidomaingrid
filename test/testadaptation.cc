#include "config.h"
#define GRIDDIM 2
#define ALUGRID_SIMPLEX

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/multidomaingrid.hh>
#include "output.hh"

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

int main(int argc, char** argv) {
#if HAVE_UG
  try {
    Dune::MPIHelper::instance(argc,argv);
    typedef Dune::UGGrid<2> AdaptableGrid;
    Dune::FieldVector<double,2> lower_left(0.0);
    Dune::FieldVector<double,2> upper_right(1.0);
    std::array<unsigned int,2> N = {{16,16}};
    std::shared_ptr<AdaptableGrid> gridPtr = Dune::StructuredGridFactory<AdaptableGrid>::createCubeGrid(
      lower_left,
      upper_right,
      N
      );
    AdaptableGrid& wgrid = *gridPtr;
    typedef Dune::MultiDomainGrid<
      AdaptableGrid,
      Dune::mdgrid::FewSubDomainsTraits<
        AdaptableGrid::dimension,
        4
        >
      > Grid;
    Grid grid(wgrid);
    typedef Grid::LeafGridView GridView;
    GridView gv = grid.leafGridView();

    grid.startSubDomainMarking();
    for (const auto& cell : elements(gv)) {
      Dune::FieldVector<AdaptableGrid::ctype,2> c = cell.geometry().global(
        Dune::ReferenceElements<AdaptableGrid::ctype,2>::general(cell.type()).position(0,0));
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

    int counter = 0;

    printStatus(grid,"adaptation",counter++);

    grid.globalRefine(1);

    printStatus(grid,"adaptation",counter++);

    for (const auto& cell : elements(gv)) {
      Dune::FieldVector<AdaptableGrid::ctype,2> c = cell.geometry().global(
        Dune::ReferenceElements<AdaptableGrid::ctype,2>::general(cell.type()).position(0,0));
      double y = c[1];
      if (y > 0.5) {
        grid.mark(1,cell);
      } else {
        grid.mark(-1,cell);
      }
    }
    grid.preAdapt();
    grid.adapt();
    grid.postAdapt();

    printStatus(grid,"adaptation",counter++);

    for (const auto& cell : elements(gv)) {
      Dune::FieldVector<AdaptableGrid::ctype,2> c = cell.geometry().global(
        Dune::ReferenceElements<AdaptableGrid::ctype,2>::general(cell.type()).position(0,0));
      double y = c[1];
      if (y > 0.5) {
        grid.mark(-1,cell);
      } else {
        grid.mark(1,cell);
      }
    }
    grid.preAdapt();
    grid.adapt();
    grid.postAdapt();

    printStatus(grid,"adaptation",counter++);

    return 0;

  } catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (std::exception &e){
    std::cerr << "std reported error: " << e.what() << std::endl;
    return 2;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 3;
  }
#else
  std::cerr << "You need UGGrid to run this test." << std::endl;
  return 77;
#endif // HAVE_UG
}
