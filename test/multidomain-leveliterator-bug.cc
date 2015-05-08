#include <config.h>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/multidomaingrid/multidomaingrid.hh>

// The grid dimension
const int dim = 2;

int main (int argc, char *argv[])
{
  try {
    Dune::MPIHelper::instance(argc,argv);

    // set up grid
    typedef Dune::YaspGrid<dim> GridType;
    Dune::FieldVector<double,dim> upperRight(1);
    Dune::array<int,2> elements = {{1, 1}};
    GridType grid(upperRight, elements);

    // set up MultiDomainGrid: a single subdomain
    typedef Dune::mdgrid::MultiDomainGrid<GridType, Dune::mdgrid::ArrayBasedTraits<dim,4,1000> > MultiDomainGridType;
    MultiDomainGridType multiDomainGrid(grid,true);

    typedef MultiDomainGridType::LeafGridView MultiDomainGridView;
    MultiDomainGridView gv = multiDomainGrid.leafGridView();

    multiDomainGrid.startSubDomainMarking();

    for (MultiDomainGridView::Codim<0>::Iterator it = gv.begin<0>(); it != gv.end<0>(); ++it)
        multiDomainGrid.addToSubDomain(0,*it);

    multiDomainGrid.preUpdateSubDomains();
    multiDomainGrid.updateSubDomains();
    multiDomainGrid.postUpdateSubDomains();

    // The actual test
    MultiDomainGridType::SubDomainGrid::Codim<0>::LevelIterator it    = multiDomainGrid.subDomain(0).lbegin<0>(0);
    MultiDomainGridType::SubDomainGrid::Codim<0>::LevelIterator endIt = multiDomainGrid.subDomain(0).lend<0>(0);

    for (; it!=endIt; ++it)
        std::cout << "foo" << std::endl;

    // should reset the iterator to the beginning
    it = multiDomainGrid.subDomain(0).lbegin<0>(0);

    for (; it!=endIt; ++it)
        std::cout << "bar" << std::endl;
    return 0;
  } catch (Dune::Exception& e) {
    std::cout << e << std::endl;
    return 1;
  } catch (...) {
    std::cout << "generic exception" << std::endl;
    return 2;
  }
}
