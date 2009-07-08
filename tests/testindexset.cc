#include "config.h"
#include <iostream>
#include <dune/common/mpihelper.hh>
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/multidomaingrid/indexset.hh>
#include <dune/grid/multidomaingrid/subdomainset.hh>

int main(int argc, char** argv) {
  try {
    Dune::MPIHelper::instance(argc,argv);
    //typedef Dune::YaspGrid<2> GridType;
    Dune::GridPtr<GridType> gridPtr("/Users/muethisn/Documents/dune/ws/dune-grid-howto/grids/unitcube2.dgf");
    GridType& grid = *gridPtr;
    grid.globalRefine(5);
    typedef GridType::LeafGridView GridView;
    GridView gv = grid.leafView();
    typedef Dune::multidomaingrid::IndexSet<GridView,Dune::multidomaingrid::IntegralTypeSubDomainSet<7> > IndexSet;
    IndexSet is(gv);
    is.update(true);
    typedef GridView::Codim<0>::Iterator Iterator;
    typedef GridView::Codim<2>::Iterator VIterator;
    typedef GridView::Codim<0>::Entity Entity;
    typedef GridView::Codim<0>::Geometry Geometry;
    for (Iterator it = gv.begin<0>(); it != gv.end<0>(); ++it) {
      const Entity& e = *it;
      IndexSet::SubDomainSet& sds = is.subDomainSet(e);
      Dune::FieldVector<GridType::ctype,2> c = e.geometry().global(Dune::GenericReferenceElements<GridType::ctype,2>::general(e.type()).position(0,0));
      double x = c[0];
      double y = c[1];
      if (x > 0.2) {
	if (y > 0.3 && y < 0.7) {
	  if (x < 0.8)
	    sds.add(1);
	  if (x > 0.6)
	    sds.add(0);
	} else {
	  sds.add(0);
	}
      }
    }
    is.update(false);
    std::vector<int> hcid(gv.indexSet().size(0),0);
    std::vector<int> sdc0(gv.indexSet().size(0),0);
    std::vector<int> sdc1(gv.indexSet().size(0),0);
    for (Iterator it = gv.begin<0>(); it != gv.end<0>(); ++it) {
      GridView::IndexSet::IndexType idx = gv.indexSet().index(*it);
      const IndexSet::SubDomainSet& sds = is.subDomainSet(*it);
      hcid[idx] = idx;
      sdc0[idx] = sds.contains(0) ? is.indexForSubDomain(0,*it) : -1;
      sdc1[idx] = sds.contains(1) ? is.indexForSubDomain(1,*it) : -1;
    }
    std::vector<int> hvid(gv.indexSet().size(2),0);
    std::vector<int> sdv0(gv.indexSet().size(2),0);
    std::vector<int> sdv1(gv.indexSet().size(2),0);
    for (VIterator it = gv.begin<2>(); it != gv.end<2>(); ++it) {
      GridView::IndexSet::IndexType idx = gv.indexSet().index(*it);
      const IndexSet::SubDomainSet& sds = is.subDomainSet(*it);
      hvid[idx] = idx;
      sdv0[idx] = sds.contains(0) ? is.subIndexForSubDomain(0,*it) : -1;
      sdv1[idx] = sds.contains(1) ? is.subIndexForSubDomain(1,*it) : -1;
    }
    Dune::VTKWriter<GridView> vtkWriter(gv);
    vtkWriter.addCellData(sdc0,"cell_subdomain0");
    vtkWriter.addCellData(sdc1,"cell_subdomain1");
    vtkWriter.addCellData(hcid,"cell_hostIndex");
    vtkWriter.addVertexData(sdv0,"vertex_subdomain0");
    vtkWriter.addVertexData(sdv1,"vertex_subdomain1");
    vtkWriter.addVertexData(hvid,"vertex_hostIndex");
    vtkWriter.write("testindexset",Dune::VTKOptions::binary);
  } catch (Dune::Exception& e) {
    std::cout << e << std::endl;
  }
  return 0;
}
