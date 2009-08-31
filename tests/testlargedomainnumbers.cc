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
    std::vector<int> sdi(gv.indexSet().size(0),0);

    for (Iterator it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
      typename GridView::IndexSet::IndexType idx = gv.indexSet().index(*it);
      const typename GridView::Grid::MDGridTraits::template Codim<0>::SubDomainSet& sds = gv.indexSet().subDomains(*it);
      sdi[idx] = gv.indexSet().index(*sds.begin(),*it);
      typename GridView::IndexSet::IndexType tmp;
    }
    /*std::vector<int> hvid(gv.indexSet().size(2),0);
    std::vector<int> sdv0(gv.indexSet().size(2),0);
    */
    std::vector<int> sdc(gv.indexSet().size(2),0);
    for (VIterator it = gv.template begin<2>(); it != gv.template end<2>(); ++it) {
      typename GridView::IndexSet::IndexType idx = gv.indexSet().index(*it);
      const typename GridView::Grid::MDGridTraits::template Codim<2>::SubDomainSet& sds = gv.indexSet().subDomains(*it);
      sdc[idx] = sds.size();
    }

    std::vector<int> borderCells(gv.indexSet().size(0),0);
    std::vector<int> borderVertices(gv.indexSet().size(2),0);

    for(; iit != iend; ++iit) {
      borderCells[gv.indexSet().index(*iit->firstCell())] = 1;
      borderCells[gv.indexSet().index(*iit->secondCell())] = 2;
    }

    Dune::VTKWriter<GridView> vtkWriter(gv);
    vtkWriter.addCellData(sdi,"cell_subdomainIndex");
    //vtkWriter.addCellData(borderCells,"borderCells");
    vtkWriter.addVertexData(sdc,"vertex_subdomainCount");
    /*vtkWriter.addVertexData(sdv1,"vertex_subdomain1");
    vtkWriter.addVertexData(hvid,"vertex_hostIndex");*/
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
    typedef Dune::MultiDomainGrid<GridType,Dune::mdgrid::MDGridTraits<GridType::dimension,1> > Grid;
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
    vtkOut(gv,"leafView",grid.leafSubDomainInterfaceBegin(0,1),grid.leafSubDomainInterfaceEnd(0,1));

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
