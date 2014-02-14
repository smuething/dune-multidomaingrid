#ifndef DUNE_MULTIDOMAINGRID_TESTS_OUTPUT_HH
#define DUNE_MULTIDOMAINGRID_TESTS_OUTPUT_HH

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <iostream>
#include <iomanip>
#include <sstream>

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
      assert((sdc0[idx] > -1) == mapper.contains(0,*it,tmp));
      assert((sdc1[idx] > -1) == mapper.contains(1,*it,tmp));
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
      assert((sdv0[idx] > -1) == mapper.contains(0,*it,tmp));
      assert((sdv1[idx] > -1) == mapper.contains(1,*it,tmp));
    }

    std::vector<int> borderCells(gv.indexSet().size(0),0);
    std::vector<int> borderVertices(gv.indexSet().size(2),0);

    for(; iit != iend; ++iit) {
      gv.grid().subDomain(0).subDomainIntersectionIterator(iit);
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
    vtkWriter.write(filename,Dune::VTK::ascii);
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
    vtkWriter.write(filename,Dune::VTK::ascii);
}


template<typename stream>
stream& setup(stream& s, std::string prefix, int counter) {
  s.str("");
  s << prefix << "_";
  if (counter >= 0)
    s << counter << "_";
  return s;
}


template<typename Grid>
void printStatus(Grid& grid, std::string prefix, int counter = -1) {
  std::ostringstream s;

  setup(s,prefix,counter) << "leafGridView";
  vtkOut(grid.leafGridView(),s.str(),grid.leafSubDomainInterfaceBegin(0,1),grid.leafSubDomainInterfaceEnd(0,1));

  for (unsigned int level = 0; level <= grid.maxLevel(); ++level) {
    setup(s,prefix,counter) << "levelGridView_" << std::setw(2) << std::setfill('0') << level;
    vtkOut(grid.levelGridView(level),s.str(),grid.levelSubDomainInterfaceBegin(0,1,level),grid.levelSubDomainInterfaceEnd(0,1,level));
  }

  setup(s,prefix,counter) << "subdomain_0";
  vtkOut2(grid.subDomain(0).leafGridView(),s.str());
   setup(s,prefix,counter) << "subdomain_1";
  vtkOut2(grid.subDomain(1).leafGridView(),s.str());
}


#endif // DUNE_MULTIDOMAINGRID_TESTS_OUTPUT_HH
