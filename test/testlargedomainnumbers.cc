#include "config.h"

#include <cstdlib>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/grid/multidomaingrid.hh>
#include "output.hh"

template<typename MultiDomainGridTraits>
void run_test(MultiDomainGridTraits traits, int sdsize)
{
  typedef Dune::YaspGrid<2> GridType;
  Dune::FieldVector<double,2> L(1.0);
  std::array<int,2> n{{32,32}};
  GridType wgrid(L,n);
  typedef Dune::MultiDomainGrid<GridType,MultiDomainGridTraits> Grid;
  Grid grid(wgrid,traits,false);
  grid.globalRefine(2);
  typedef typename Grid::LeafGridView GridView;
  GridView gv = grid.leafGridView();
  grid.startSubDomainMarking();
  typename GridView::IndexSet::IndexType sd=0;
  int c = 0;
  for (const auto& cell : elements(gv)) {
    if (c++ == sdsize) {
      ++sd;
      c = 0;
    }
    grid.addToSubDomain(sd,cell);
  }
  std::cout << "Number of subdomains: " << sd + 1 << std::endl;
  grid.preUpdateSubDomains();
  grid.updateSubDomains();
  grid.postUpdateSubDomains();
  vtkOut(gv,"largedomainnumbers_leafGridView",grid.leafSubDomainInterfaceBegin(0,1),grid.leafSubDomainInterfaceEnd(0,1));
}


int main(int argc, char** argv) {
  try {
    Dune::MPIHelper::instance(argc,argv);

    const int sdsize = argc > 1 ? atoi(argv[1]) : 10;

    run_test(Dune::mdgrid::ArrayBasedTraits<2,1,65536>(),sdsize);
    run_test(Dune::mdgrid::DynamicSubDomainCountTraits<2,1>(65536),sdsize);

  } catch (Dune::Exception& e) {
    std::cout << e << std::endl;
  }
  return 0;
}
