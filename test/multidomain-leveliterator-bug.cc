#include <config.h>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/multidomaingrid/multidomaingrid.hh>

// The grid dimension
const int dim = 2;

int main (int argc, char *argv[])
{
    // set up grid
    typedef Dune::YaspGrid<dim> GridType;
    Dune::FieldVector<double,dim> upperRight(1);
    Dune::FieldVector<int,2> elements(1);
    Dune::FieldVector<bool,2> periodic(false);
    GridType grid(upperRight, elements, periodic, 0);
    
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
}
