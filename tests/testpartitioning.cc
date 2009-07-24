#include "config.h"
#include <iostream>
#include <dune/common/mpihelper.hh>
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/multidomaingrid.hh>

template<typename GridView, typename InterfaceIterator>
void vtkOut(GridView gv,std::string filename, InterfaceIterator iit, InterfaceIterator iend) {
    typedef typename GridView::template Codim<0>::Iterator Iterator;
    typedef typename GridView::template Codim<2>::Iterator VIterator;
    std::vector<int> hcid(gv.indexSet().size(0),0);
    std::vector<int> sdc0(gv.indexSet().size(0),0);

    std::vector<int> sdc1(gv.indexSet().size(0),0);
    for (Iterator it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
      typename GridView::IndexSet::IndexType idx = gv.indexSet().index(*it);
      const typename GridView::Grid::SubDomainSet& sds = gv.indexSet().subDomains(*it);
      hcid[idx] = idx;
      sdc0[idx] = sds.contains(0) ? gv.indexSet().index(0,*it) : -1;
      sdc1[idx] = sds.contains(1) ? gv.indexSet().index(1,*it) : -1;
    }
    std::vector<int> hvid(gv.indexSet().size(2),0);
    std::vector<int> sdv0(gv.indexSet().size(2),0);
    std::vector<int> sdv1(gv.indexSet().size(2),0);
    for (VIterator it = gv.template begin<2>(); it != gv.template end<2>(); ++it) {
      typename GridView::IndexSet::IndexType idx = gv.indexSet().index(*it);
      const typename GridView::Grid::SubDomainSet& sds = gv.indexSet().subDomains(*it);
      hvid[idx] = idx;
      sdv0[idx] = sds.contains(0) ? gv.indexSet().index(0,*it) : -1;
      sdv1[idx] = sds.contains(1) ? gv.indexSet().index(1,*it) : -1;
    }

    std::vector<int> borderCells(gv.indexSet().size(0),0);
    std::vector<int> borderVertices(gv.indexSet().size(2),0);

    for(; iit != iend; ++iit) {
      borderCells[gv.indexSet().index(*iit->firstCell())] = 1;
      borderCells[gv.indexSet().index(*iit->secondCell())] = 2;
    }

    Dune::VTKWriter<GridView> vtkWriter(gv);
    vtkWriter.addCellData(sdc0,"cell_subdomain0");
    vtkWriter.addCellData(sdc1,"cell_subdomain1");
    vtkWriter.addCellData(hcid,"cell_hostIndex");
    vtkWriter.addCellData(borderCells,"borderCells");
    vtkWriter.addVertexData(sdv0,"vertex_subdomain0");
    vtkWriter.addVertexData(sdv1,"vertex_subdomain1");
    vtkWriter.addVertexData(hvid,"vertex_hostIndex");
    vtkWriter.write(filename,Dune::VTKOptions::binary);
}

template<typename GridView>
void vtkOut2(GridView gv,std::string filename) {
    typedef typename GridView::template Codim<0>::Iterator Iterator;
    typedef typename GridView::template Codim<2>::Iterator VIterator;
    std::vector<int> cid(gv.indexSet().size(0),0);
    for (Iterator it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
      typename GridView::IndexSet::IndexType idx = gv.indexSet().index(*it);
      cid[idx] = idx;
    }
    std::vector<int> vid(gv.indexSet().size(2),0);
    for (VIterator it = gv.template begin<2>(); it != gv.template end<2>(); ++it) {
      typename GridView::IndexSet::IndexType idx = gv.indexSet().index(*it);
      vid[idx] = idx;
    }
    Dune::VTKWriter<GridView> vtkWriter(gv);
    vtkWriter.addCellData(cid,"cellIndex");
    vtkWriter.addVertexData(vid,"vertexIndex");
    vtkWriter.write(filename,Dune::VTKOptions::binary);
}

int main(int argc, char** argv) {
  try {
    Dune::MPIHelper::instance(argc,argv);
    //typedef Dune::YaspGrid<2> GridType;
    Dune::GridPtr<GridType> gridPtr("/Users/muethisn/Documents/dune/ws/dune-grid-howto/grids/unitcube2.dgf");
    GridType& wgrid = *gridPtr;
    Dune::MultiDomainGrid<GridType> grid(wgrid);
    grid.globalRefine(5);
    typedef Dune::MultiDomainGrid<GridType>::LeafGridView GridView;
    GridView gv = grid.leafView();
    typedef GridView::Codim<0>::Iterator Iterator;
    typedef GridView::Codim<2>::Iterator VIterator;
    typedef GridView::Codim<0>::Entity Entity;
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

    const Dune::MultiDomainGrid<GridType>::SubDomainGrid& sd0 = grid.subDomain(0);
    vtkOut(gv,"leafView",grid.leafSubDomainInterfaceBegin(0,1),grid.leafSubDomainInterfaceEnd(0,1));
    vtkOut(grid.levelView(0),"levelView0",grid.levelSubDomainInterfaceBegin(0,1,0),grid.levelSubDomainInterfaceEnd(0,1,0));
    vtkOut(grid.levelView(1),"levelView1",grid.levelSubDomainInterfaceBegin(0,1,1),grid.levelSubDomainInterfaceEnd(0,1,1));
    vtkOut(grid.levelView(2),"levelView2",grid.levelSubDomainInterfaceBegin(0,1,2),grid.levelSubDomainInterfaceEnd(0,1,2));
    vtkOut(grid.levelView(3),"levelView3",grid.levelSubDomainInterfaceBegin(0,1,3),grid.levelSubDomainInterfaceEnd(0,1,3));
    vtkOut(grid.levelView(4),"levelView4",grid.levelSubDomainInterfaceBegin(0,1,4),grid.levelSubDomainInterfaceEnd(0,1,4));
    vtkOut(grid.levelView(5),"levelView5",grid.levelSubDomainInterfaceBegin(0,1,5),grid.levelSubDomainInterfaceEnd(0,1,5));

    vtkOut2(grid.subDomain(0).leafView(),"subdomain0");
    vtkOut2(grid.subDomain(1).leafView(),"subdomain1");

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
