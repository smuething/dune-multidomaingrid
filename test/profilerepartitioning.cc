#include "config.h"

#include <dune/common/mpihelper.hh>
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#include <dune/grid/multidomaingrid.hh>
#include "output.hh"


int main(int argc, char** argv) {
  try {
    Dune::MPIHelper::instance(argc,argv);
    //typedef Dune::YaspGrid<2> GridType;
    Dune::GridPtr<GridType> gridPtr("/Users/muethisn/Documents/dune/ws/dune-grid-howto/grids/unitcube2.dgf");
    GridType& wgrid = *gridPtr;
    typedef Dune::MultiDomainGrid<GridType,Dune::mdgrid::FewSubDomainsTraits<GridType::dimension,4> > Grid;
    Grid grid(wgrid);
    grid.globalRefine(8);
    typedef Grid::LeafGridView GridView;
    GridView gv = grid.leafView();
    typedef GridView::Codim<0>::Iterator Iterator;
    typedef GridView::Codim<2>::Iterator VIterator;
    typedef GridView::Codim<0>::Entity Entity;
    typedef GridView::Codim<2>::Entity Vertex;
    typedef GridView::Codim<0>::Geometry Geometry;
    for (int i = 0; i < 1000; ++i) {
      if (i % 10 == 0)
        std::cout << i << std::endl;
      grid.startSubDomainMarking();
      for (Iterator it = gv.begin<0>(); it != gv.end<0>(); ++it) {
        const Entity& e = *it;
        //IndexSet::SubDomainSet& sds = is.subDomainSet(e);
        Dune::FieldVector<GridType::ctype,2> c = e.geometry().global(Dune::GenericReferenceElements<GridType::ctype,2>::general(e.type()).position(0,0));
        grid.removeFromAllSubDomains(e);
        double x = c[0];
        double y = c[1];
        if (x > 0.2) {
          if (y > 0.3 && y < 0.7) {
            if (x < 0.8)
              grid.addToSubDomain(1,e);
            if (x > 0.6)
              grid.addToSubDomain(0,e);
          } else {
            grid.addToSubDomain(0,e);
          }
        }
      }
      grid.preUpdateSubDomains();
      grid.updateSubDomains();
      grid.postUpdateSubDomains();

      grid.startSubDomainMarking();
      for (Iterator it = gv.begin<0>(); it != gv.end<0>(); ++it) {
        const Entity& e = *it;
        const Grid::MDGridTraits::Codim<0>::SubDomainSet& domains = gv.indexSet().subDomains(e);
        if (domains.contains(0) && domains.contains(1)) {
          grid.removeFromSubDomain(0,e);
          grid.removeFromSubDomain(1,e);
        } else if (domains.contains(0)) {
          grid.assignToSubDomain(1,e);
        } else if (domains.contains(1)) {
          grid.assignToSubDomain(0,e);
        } else {
          grid.addToSubDomain(0,e);
          grid.addToSubDomain(1,e);
        }
      }

      grid.preUpdateSubDomains();
      grid.updateSubDomains();
      grid.postUpdateSubDomains();
    }

  } catch (Dune::Exception& e) {
    std::cout << e << std::endl;
  }
  return 0;
}
