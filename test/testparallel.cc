#include "config.h"

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/alugrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/multidomaingrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <iterator>

template<typename GV, typename DataVector>
class RankTransfer
  : public Dune::CommDataHandleIF<RankTransfer<GV,DataVector>,
                                  int
                                  >
{

public:

  bool contains(int dim, int codim) const
  {
    return codim == _codim;
  }

  bool fixedsize(int dim, int codim) const
  {
    return true;
  }

  template<typename Entity>
  std::size_t size(const Entity& e) const
  {
    return 1;
  }

  template<typename MessageBufferImp, typename Entity>
  void gather(MessageBufferImp& buf, const Entity& e) const
  {
    buf.write(Dune::MPIHelper::getCollectiveCommunication().rank());
  }

  template<typename MessageBufferImp, typename Entity>
  void scatter(MessageBufferImp& buf, const Entity& e, std::size_t n)
  {
    int i;
    buf.read(i);
    _data[_gv.indexSet().index(e)] |= 1 << i;
  }

  RankTransfer(const GV& gv, DataVector& data, int codim)
    : _gv(gv)
    , _data(data)
    , _codim(codim)
  {}

private:
  GV _gv;
  DataVector& _data;
  const int _codim;

};

template<typename HostGrid>
void testGrid(HostGrid& hostgrid, std::string prefix, Dune::MPIHelper& mpihelper)
{
  const int dim = HostGrid::dimension;
  typedef Dune::MultiDomainGrid<HostGrid,Dune::mdgrid::FewSubDomainsTraits<2,8> > MDGrid;
  //typedef Dune::MultiDomainGrid<HostGrid,Dune::mdgrid::ArrayBasedTraits<2,8,8> > MDGrid;

  hostgrid.leafIndexSet().index(*hostgrid.leafGridView().template begin<0>());

  MDGrid grid(hostgrid,true);
  typedef typename MDGrid::LeafGridView MDGV;
  typedef typename MDGrid::SubDomainIndex SubDomainIndex;
  MDGV mdgv = grid.leafGridView();

  grid.startSubDomainMarking();

  for (const auto& cell : elements(mdgv))
    {
      if (cell.partitionType() != Dune::InteriorEntity)
        continue;
      SubDomainIndex subdomain = 0;
      if (cell.geometry().center()[0] > 0.5)
        subdomain += 1;
      if (cell.geometry().center()[1] > 0.5)
        subdomain += 2;
      grid.addToSubDomain(subdomain,cell);
      if (cell.geometry().center()[1] > 0.5)
        grid.addToSubDomain(4,cell);
      if (cell.geometry().center()[0] > 0.5)
        grid.addToSubDomain(5,cell);
      if (cell.geometry().center()[0] > 0.5 && cell.geometry().center()[1] > 0.5)
        grid.addToSubDomain(6,cell);
      if (cell.geometry().center()[0] > 0.5 && cell.geometry().center()[1] < 0.5)
        grid.addToSubDomain(7,cell);
    }

  grid.preUpdateSubDomains();
  grid.updateSubDomains();
  grid.postUpdateSubDomains();

  for (const auto& cell : elements(mdgv))
    {
      std::cout << mpihelper.rank() << ": " << std::setw(2) << mdgv.indexSet().index(cell) << "  " << cell.geometry().center()
                << " subdomains:";
      for (auto subdomain : mdgv.indexSet().subDomains(cell))
        std::cout << " " << subdomain;
      std::cout << std::endl;
    }

  std::cout << std::endl;
  for (const auto& vertex : vertices(mdgv))
    {
      std::cout << mpihelper.rank() << ": " << std::setw(2) << mdgv.indexSet().index(vertex) << "  " << vertex.geometry().center()
                << " subdomains:";
      for (auto subdomain : mdgv.indexSet().subDomains(vertex))
        std::cout << " " << subdomain;
      std::cout << std::endl;
    }

  for (SubDomainIndex s = 0; s < 8; ++s)
    {
      typedef typename MDGrid::SubDomainGrid SDGrid;
      typedef typename SDGrid::LeafGridView SDGV;
      const SDGrid& sdgrid = grid.subDomain(s);
      SDGV sdgv = sdgrid.leafGridView();

      typedef std::vector<std::size_t> DataVector;
      DataVector celldata(sdgv.size(0));
      std::fill(celldata.begin(),celldata.end(),1 << mpihelper.rank());

      typedef RankTransfer<SDGV,DataVector> DataHandle;
      DataHandle celldatahandle(sdgv,celldata,0);

      sdgrid.communicate(celldatahandle,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);

      DataVector nodedata(sdgv.size(dim));
      std::fill(nodedata.begin(),nodedata.end(),1 << mpihelper.rank());

      DataHandle nodedatahandle(sdgv,nodedata,dim);

      sdgrid.communicate(nodedatahandle,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);

      bool dummy = false;
      for (const auto& cell : elements(sdgv))
        {
          for (const auto& intersection : intersections(sdgv,cell))
            dummy |= intersection.boundary() || intersection.neighbor();
        }

      std::stringstream sstr;
      sstr << prefix << (dummy ? "_subdomain_" : "_subdomain_") << (char)('a' + s);
      Dune::VTKWriter<SDGV> vtkwriter(sdgv);
      std::vector<int> cellIndices(sdgv.indexSet().size(0));
      std::iota(cellIndices.begin(),cellIndices.end(),int(0));
      std::cout << cellIndices.size() << " " << sdgv.indexSet().size(0) << std::endl;
      vtkwriter.addCellData(cellIndices,"cell index");
      vtkwriter.addCellData(celldata,"cell rank");
      vtkwriter.addVertexData(nodedata,"node rank");
      vtkwriter.write(sstr.str());
    }
}


int main(int argc, char** argv)
{
  try {
    const int dim = 2;
    Dune::MPIHelper& mpihelper = Dune::MPIHelper::instance(argc,argv);

    std::size_t N = 32;
    std::size_t overlap = 1;

    if (argc < 2)
      {
        std::cerr << "Usage: " << argv[0] <<
          " <number of elements per dim> <overlap>" << std::endl;
        std::cerr << std::endl
                  << "defaulting to " << N << " " << overlap << std::endl;
      }
    else
      {
        N = atoi(argv[1]);
        overlap = atoi(argv[2]);
      }

    /* helper code for attaching debugger
    if (mpihelper.rank() == 0)
      {
        volatile int i = 0;
        std::cout << getpid() << std::endl;
        while (i == 0)
          sleep(5);
      }
    */

    {
      const Dune::FieldVector<double,dim> h(1.0);
      Dune::array<int,dim> s;
      std::fill(s.begin(), s.end(), N);
      std::bitset<dim> p(false);
      typedef Dune::YaspGrid<dim> HostGrid;
      HostGrid hostgrid(h,s,p,overlap);

      testGrid(hostgrid,"YaspGrid_2",mpihelper);
    }
    /*
    {
      typedef Dune::ALUSimplexGrid<2,2> HostGrid;

      std::vector<int> boundaryid;
      std::vector<int> elementid;
      std::shared_ptr<HostGrid> gridptr
        (Dune::GmshReader<HostGrid>::read("square.msh",boundaryid, elementid,
                                          true, false));

      gridptr->loadBalance();

      testGrid(*gridptr,"AluSimplexGrid_2_2",mpihelper);
    }
    */
  }
  catch (std::exception & e) {
    std::cout << "STL ERROR: " << e.what() << std::endl;
    return 1;
  }
  catch (Dune::Exception & e) {
    std::cout << "DUNE ERROR: " << e.what() << std::endl;
    return 1;
  }
  catch (...) {
    std::cout << "Unknown ERROR" << std::endl;
    return 1;
  }
  return 0;
}
