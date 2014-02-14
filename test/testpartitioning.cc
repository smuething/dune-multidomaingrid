#include "config.h"

#include <dune/common/parallel/mpihelper.hh>
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
    grid.globalRefine(5);
    typedef Grid::LeafGridView GridView;
    GridView gv = grid.leafView();
    typedef GridView::Codim<0>::Iterator Iterator;
    typedef GridView::Codim<2>::Iterator VIterator;
    typedef GridView::Codim<0>::Entity Entity;
    typedef GridView::Codim<2>::Entity Vertex;
    typedef GridView::Codim<0>::Geometry Geometry;
    grid.startSubDomainMarking();
    for (Iterator it = gv.begin<0>(); it != gv.end<0>(); ++it) {
      const Entity& e = *it;
      //IndexSet::SubDomainSet& sds = is.subDomainSet(e);
      Dune::FieldVector<GridType::ctype,2> c = e.geometry().global(Dune::GenericReferenceElements<GridType::ctype,2>::general(e.type()).position(0,0));
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

    Grid::SubDomainGrid::LeafGridView::Codim<0>::EntityPointer p0 = grid.subDomain(0).leafView().begin<0>();
    const Entity& e0 = grid.subDomain(0).multiDomainEntity(*p0);
    const Entity& e02 = grid.multiDomainEntity(*p0);
    Grid::SubDomainGrid::LeafGridView::Codim<0>::EntityPointer p02 = grid.subDomain(0).subDomainEntityPointer(e0);

    Grid::SubDomainGrid::LeafGridView::Codim<2>::EntityPointer p1 = grid.subDomain(0).leafView().begin<2>();
    const Vertex& e1 = grid.subDomain(0).multiDomainEntity(*p1);
    const Vertex& e12 = grid.multiDomainEntity(*p1);

    printStatus(grid,"partitioning1");

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

    printStatus(grid,"partitioning2");

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
