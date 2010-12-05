#include "config.h"

#include <dune/common/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/multidomaingrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <iterator>
#include <boost/iterator/counting_iterator.hpp>

int main(int argc, char** argv)
{
  try {
    const int dim = 2;
    Dune::MPIHelper& mpihelper = Dune::MPIHelper::instance(argc,argv);
    const Dune::FieldVector<double,dim> h(1.0);
    const Dune::FieldVector<int,dim> s(atoi(argv[1]));
    const Dune::FieldVector<bool,dim> p(false);
    typedef Dune::YaspGrid<dim> HostGrid;
    HostGrid hostgrid(Dune::MPIHelper::getCommunicator(),h,s,p,atoi(argv[2]));
    //typedef HostGrid::LeafGridView GV;
    //GV gv = grid.leafView();

    //typedef Dune::MultiDomainGrid<HostGrid,Dune::mdgrid::FewSubDomainsTraits<2,16> > MDGrid;
    typedef Dune::MultiDomainGrid<HostGrid,Dune::mdgrid::ArrayBasedTraits<2,8,8> > MDGrid;

    MDGrid grid(hostgrid,true);
    typedef MDGrid::LeafGridView MDGV;
    typedef MDGrid::SubDomainIndexType SubDomainIndexType;
    MDGV mdgv = grid.leafView();

    grid.startSubDomainMarking();

    for (auto it = mdgv.begin<0>(); it != mdgv.end<0>(); ++it)
      {
        if (it->partitionType() != Dune::InteriorEntity)
          continue;
      SubDomainIndexType subdomain = 0;
      if (it->geometry().center()[0] > 0.5)
        subdomain += 1;
      if (it->geometry().center()[1] > 0.5)
        subdomain += 2;
      grid.addToSubDomain(subdomain,*it);
      if (it->geometry().center()[1] > 0.5)
        grid.addToSubDomain(4,*it);
      if (it->geometry().center()[0] > 0.5)
        grid.addToSubDomain(5,*it);
      if (it->geometry().center()[0] > 0.5 && it->geometry().center()[1] > 0.5)
        grid.addToSubDomain(6,*it);
      if (it->geometry().center()[0] > 0.5 && it->geometry().center()[1] < 0.5)
        grid.addToSubDomain(7,*it);
    }

    grid.preUpdateSubDomains();
    grid.updateSubDomains();
    grid.postUpdateSubDomains();

    for (auto it = mdgv.begin<0>(); it != mdgv.end<0>(); ++it)
      {
        std::cout << mpihelper.rank() << ": " << std::setw(2) << mdgv.indexSet().index(*it) << "  " << it->geometry().center()
                  << " subdomains:";
        auto end = mdgv.indexSet().subDomains(*it).end();
        for (auto sdit = mdgv.indexSet().subDomains(*it).begin(); sdit != end; ++sdit)
          std::cout << " " << *sdit;
        std::cout << std::endl;
      }
    std::cout << std::endl;
    for (auto it = mdgv.begin<dim>(); it != mdgv.end<dim>(); ++it)
      {
        std::cout << mpihelper.rank() << ": " << std::setw(2) << mdgv.indexSet().index(*it) << "  " << it->geometry().center()
                  << " subdomains:";
        auto end = mdgv.indexSet().subDomains(*it).end();
        for (auto sdit = mdgv.indexSet().subDomains(*it).begin(); sdit != end; ++sdit)
          std::cout << " " << *sdit;
        std::cout << std::endl;
      }
    for (SubDomainIndexType s = 0; s < 8; ++s)
      {
        typedef MDGrid::SubDomainGrid SDGrid;
        typedef SDGrid::LeafGridView SDGV;
        const SDGrid& sdgrid = grid.subDomain(s);
        SDGV sdgv = sdgrid.leafView();
        std::stringstream sstr;
        sstr << "subdomain_" << (char)('a' + s);
        Dune::VTKWriter<SDGV> vtkwriter(sdgv);
        std::vector<int> cellIndices;
        std::copy(boost::counting_iterator<int>(0),
                  boost::counting_iterator<int>(sdgv.indexSet().size(0)),
                  std::back_inserter(cellIndices));
        std::cout << cellIndices.size() << " " << sdgv.indexSet().size(0) << std::endl;
        vtkwriter.addCellData(cellIndices,"cell indices");
        vtkwriter.write(sstr.str());
      }
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
