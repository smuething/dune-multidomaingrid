#include "config.h"
#define GRIDDIM 2
#define ALUGRID_SIMPLEX

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/multidomaingrid.hh>
#include "output.hh"

int main(int argc, char** argv) {
  try {
    Dune::MPIHelper::instance(argc,argv);
    typedef Dune::ALUSimplexGrid<2,2> AdaptableGridType;
    Dune::GridPtr<GridType> gridPtr("../../dune-grid-howto/grids/unitcube2.dgf");
    GridType& wgrid = *gridPtr;
    typedef Dune::MultiDomainGrid<GridType,Dune::mdgrid::FewSubDomainsTraits<GridType::dimension,4> > Grid;
    Grid grid(wgrid);
    grid.globalRefine(5);
    typedef Grid::LeafGridView GridView;
    GridView gv = grid.leafGridView();
    typedef GridView::Codim<0>::Iterator Iterator;
    typedef GridView::Codim<2>::Iterator VIterator;
    typedef GridView::Codim<0>::Entity Entity;
    typedef GridView::Codim<0>::Geometry Geometry;
    grid.startSubDomainMarking();
    for (Iterator it = gv.begin<0>(); it != gv.end<0>(); ++it) {
      const Entity& e = *it;
      //IndexSet::SubDomainSet& sds = is.subDomainSet(e);
      Dune::FieldVector<GridType::ctype,2> c = e.geometry().global(
        Dune::ReferenceElements<GridType::ctype,2>::general(e.type()).position(0,0));
      double x = c[0];
      double y = c[1];
      if (x > 0.2) {
        if (y > 0.3 && y < 0.7) {
          if (x < 0.8)
            grid.addToSubDomain(1,e);
          else // if (x > 0.6)
            grid.addToSubDomain(0,e);
        } else {
            grid.addToSubDomain(0,e);
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

    for (Iterator it = gv.begin<0>(); it != gv.end<0>(); ++it) {
      const Entity& e = *it;
      Dune::FieldVector<GridType::ctype,2> c = e.geometry().global(
        Dune::ReferenceElements<GridType::ctype,2>::general(e.type()).position(0,0));
      double x = c[0];
      double y = c[1];
      if (y > 0.5) {
        grid.mark(1,e);
      } else {
        grid.mark(-1,e);
      }
    }
    grid.preAdapt();
    grid.adapt();
    grid.postAdapt();

    printStatus(grid,"adaptation",counter++);

    for (Iterator it = gv.begin<0>(); it != gv.end<0>(); ++it) {
      const Entity& e = *it;
      Dune::FieldVector<GridType::ctype,2> c = e.geometry().global(
        Dune::ReferenceElements<GridType::ctype,2>::general(e.type()).position(0,0));
      double x = c[0];
      double y = c[1];
      if (y > 0.5) {
        grid.mark(-1,e);
      } else {
        grid.mark(1,e);
      }
    }
    grid.preAdapt();
    grid.adapt();
    grid.postAdapt();

    printStatus(grid,"adaptation",counter++);

  } catch (Dune::Exception& e) {
    std::cout << e << std::endl;
  }
  return 0;
}
