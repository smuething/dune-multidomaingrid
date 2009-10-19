#include "config.h"
#include <iostream>
#include <dune/common/mpihelper.hh>
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/multidomaingrid.hh>

template<int dim>
struct AllLayout
{
  bool contains(Dune::GeometryType gt) {
    return true;
  }
};

template<typename GridView, typename InterfaceIterator>
void vtkOut(GridView gv,std::string filename, InterfaceIterator iit, InterfaceIterator iend) {

    Dune::MultiDomainMCMGMapper<GridView,AllLayout> mapper(gv);

    typedef typename GridView::template Codim<0>::Iterator Iterator;
    typedef typename GridView::template Codim<2>::Iterator VIterator;
    std::vector<int> hcid(gv.indexSet().size(0),0);
    std::vector<int> sdc0(gv.indexSet().size(0),0);

    std::vector<int> sdc1(gv.indexSet().size(0),0);
    for (Iterator it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
      typename GridView::IndexSet::IndexType idx = gv.indexSet().index(*it);
      const typename GridView::Grid::MDGridTraits::template Codim<0>::SubDomainSet& sds = gv.indexSet().subDomains(*it);
      hcid[idx] = idx;
      sdc0[idx] = sds.contains(0) ? gv.indexSet().index(0,*it) : -1;
      sdc1[idx] = sds.contains(1) ? gv.indexSet().index(1,*it) : -1;
      typename GridView::IndexSet::IndexType tmp;
      assert(sdc0[idx] > -1 == mapper.contains(0,*it,tmp));
      assert(sdc1[idx] > -1 == mapper.contains(1,*it,tmp));
    }
    std::vector<int> hvid(gv.indexSet().size(2),0);
    std::vector<int> sdv0(gv.indexSet().size(2),0);
    std::vector<int> sdv1(gv.indexSet().size(2),0);
    for (VIterator it = gv.template begin<2>(); it != gv.template end<2>(); ++it) {
      typename GridView::IndexSet::IndexType idx = gv.indexSet().index(*it);
      const typename GridView::Grid::MDGridTraits::template Codim<2>::SubDomainSet& sds = gv.indexSet().subDomains(*it);
      hvid[idx] = idx;
      sdv0[idx] = sds.contains(0) ? gv.indexSet().index(0,*it) : -1;
      sdv1[idx] = sds.contains(1) ? gv.indexSet().index(1,*it) : -1;
      typename GridView::IndexSet::IndexType tmp;
      assert(sdv0[idx] > -1 == mapper.contains(0,*it,tmp));
      assert(sdv1[idx] > -1 == mapper.contains(1,*it,tmp));
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
    typedef Dune::ALUSimplexGrid<2> GridType;
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

    vtkOut(gv,"0_leafView",grid.leafSubDomainInterfaceBegin(0,1),grid.leafSubDomainInterfaceEnd(0,1));
    vtkOut(grid.levelView(0),"0_levelView0",grid.levelSubDomainInterfaceBegin(0,1,0),grid.levelSubDomainInterfaceEnd(0,1,0));
    vtkOut(grid.levelView(1),"0_levelView1",grid.levelSubDomainInterfaceBegin(0,1,1),grid.levelSubDomainInterfaceEnd(0,1,1));
    vtkOut(grid.levelView(2),"0_levelView2",grid.levelSubDomainInterfaceBegin(0,1,2),grid.levelSubDomainInterfaceEnd(0,1,2));
    vtkOut(grid.levelView(3),"0_levelView3",grid.levelSubDomainInterfaceBegin(0,1,3),grid.levelSubDomainInterfaceEnd(0,1,3));
    vtkOut(grid.levelView(4),"0_levelView4",grid.levelSubDomainInterfaceBegin(0,1,4),grid.levelSubDomainInterfaceEnd(0,1,4));
    vtkOut(grid.levelView(5),"0_levelView5",grid.levelSubDomainInterfaceBegin(0,1,5),grid.levelSubDomainInterfaceEnd(0,1,5));

    vtkOut2(grid.subDomain(0).leafView(),"0_subdomain0");
    vtkOut2(grid.subDomain(1).leafView(),"0_subdomain1");

    grid.globalrefine(1);

    vtkOut(gv,"1_leafView",grid.leafSubDomainInterfaceBegin(0,1),grid.leafSubDomainInterfaceEnd(0,1));
    vtkOut(grid.levelView(0),"1_levelView0",grid.levelSubDomainInterfaceBegin(0,1,0),grid.levelSubDomainInterfaceEnd(0,1,0));
    vtkOut(grid.levelView(1),"1_levelView1",grid.levelSubDomainInterfaceBegin(0,1,1),grid.levelSubDomainInterfaceEnd(0,1,1));
    vtkOut(grid.levelView(2),"1_levelView2",grid.levelSubDomainInterfaceBegin(0,1,2),grid.levelSubDomainInterfaceEnd(0,1,2));
    vtkOut(grid.levelView(3),"1_levelView3",grid.levelSubDomainInterfaceBegin(0,1,3),grid.levelSubDomainInterfaceEnd(0,1,3));
    vtkOut(grid.levelView(4),"1_levelView4",grid.levelSubDomainInterfaceBegin(0,1,4),grid.levelSubDomainInterfaceEnd(0,1,4));
    vtkOut(grid.levelView(5),"1_levelView5",grid.levelSubDomainInterfaceBegin(0,1,5),grid.levelSubDomainInterfaceEnd(0,1,5));
    vtkOut(grid.levelView(6),"1_levelView6",grid.levelSubDomainInterfaceBegin(0,1,5),grid.levelSubDomainInterfaceEnd(0,1,5));

    vtkOut2(grid.subDomain(0).leafView(),"1_subdomain0");
    vtkOut2(grid.subDomain(1).leafView(),"1_subdomain1");

    for (Iterator it = gv.begin<0>(); it != gv.end<0>(); ++it) {
      const Entity& e = *it;
      Dune::FieldVector<GridType::ctype,2> c = e.geometry().global(Dune::GenericReferenceElements<GridType::ctype,2>::general(e.type()).position(0,0));
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

    vtkOut(gv,"2_leafView",grid.leafSubDomainInterfaceBegin(0,1),grid.leafSubDomainInterfaceEnd(0,1));
    vtkOut(grid.levelView(0),"2_levelView0",grid.levelSubDomainInterfaceBegin(0,1,0),grid.levelSubDomainInterfaceEnd(0,1,0));
    vtkOut(grid.levelView(1),"2_levelView1",grid.levelSubDomainInterfaceBegin(0,1,1),grid.levelSubDomainInterfaceEnd(0,1,1));
    vtkOut(grid.levelView(2),"2_levelView2",grid.levelSubDomainInterfaceBegin(0,1,2),grid.levelSubDomainInterfaceEnd(0,1,2));
    vtkOut(grid.levelView(3),"2_levelView3",grid.levelSubDomainInterfaceBegin(0,1,3),grid.levelSubDomainInterfaceEnd(0,1,3));
    vtkOut(grid.levelView(4),"2_levelView4",grid.levelSubDomainInterfaceBegin(0,1,4),grid.levelSubDomainInterfaceEnd(0,1,4));
    vtkOut(grid.levelView(5),"2_levelView5",grid.levelSubDomainInterfaceBegin(0,1,5),grid.levelSubDomainInterfaceEnd(0,1,5));
    vtkOut(grid.levelView(6),"2_levelView6",grid.levelSubDomainInterfaceBegin(0,1,5),grid.levelSubDomainInterfaceEnd(0,1,5));
    vtkOut(grid.levelView(7),"2_levelView7",grid.levelSubDomainInterfaceBegin(0,1,5),grid.levelSubDomainInterfaceEnd(0,1,5));

    vtkOut2(grid.subDomain(0).leafView(),"2_subdomain0");
    vtkOut2(grid.subDomain(1).leafView(),"2_subdomain1");

    for (Iterator it = gv.begin<0>(); it != gv.end<0>(); ++it) {
      const Entity& e = *it;
      Dune::FieldVector<GridType::ctype,2> c = e.geometry().global(Dune::GenericReferenceElements<GridType::ctype,2>::general(e.type()).position(0,0));
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

    vtkOut(gv,"3_leafView",grid.leafSubDomainInterfaceBegin(0,1),grid.leafSubDomainInterfaceEnd(0,1));
    vtkOut(grid.levelView(0),"3_levelView0",grid.levelSubDomainInterfaceBegin(0,1,0),grid.levelSubDomainInterfaceEnd(0,1,0));
    vtkOut(grid.levelView(1),"3_levelView1",grid.levelSubDomainInterfaceBegin(0,1,1),grid.levelSubDomainInterfaceEnd(0,1,1));
    vtkOut(grid.levelView(2),"3_levelView2",grid.levelSubDomainInterfaceBegin(0,1,2),grid.levelSubDomainInterfaceEnd(0,1,2));
    vtkOut(grid.levelView(3),"3_levelView3",grid.levelSubDomainInterfaceBegin(0,1,3),grid.levelSubDomainInterfaceEnd(0,1,3));
    vtkOut(grid.levelView(4),"3_levelView4",grid.levelSubDomainInterfaceBegin(0,1,4),grid.levelSubDomainInterfaceEnd(0,1,4));
    vtkOut(grid.levelView(5),"3_levelView5",grid.levelSubDomainInterfaceBegin(0,1,5),grid.levelSubDomainInterfaceEnd(0,1,5));
    vtkOut(grid.levelView(6),"3_levelView6",grid.levelSubDomainInterfaceBegin(0,1,5),grid.levelSubDomainInterfaceEnd(0,1,5));

    vtkOut2(grid.subDomain(0).leafView(),"3_subdomain0");
    vtkOut2(grid.subDomain(1).leafView(),"3_subdomain1");

  } catch (Dune::Exception& e) {
    std::cout << e << std::endl;
  }
  return 0;
}
