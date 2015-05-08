#include "config.h"

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/multidomaingrid.hh>
#include <iostream>
#include <cassert>

template<typename HostGrid>
void run_test(HostGrid& hostgrid)
{
  hostgrid.globalRefine(1);
  typedef Dune::MultiDomainGrid<HostGrid,Dune::mdgrid::ArrayBasedTraits<2,8,8> > MDGrid;
  MDGrid mdgrid(hostgrid,true);

  typedef typename MDGrid::LeafGridView MDGV;

  MDGV mdgv = mdgrid.leafGridView();

  mdgrid.startSubDomainMarking();

  for(const auto& cell : elements(mdgv))
      mdgrid.addToSubDomain(0,cell);

  mdgrid.preUpdateSubDomains();
  mdgrid.updateSubDomains();
  mdgrid.postUpdateSubDomains();

  typedef typename MDGrid::SubDomainGrid SDGrid;
  const SDGrid& sdgrid = mdgrid.subDomain(0);
  typename SDGrid::LeafGridView sdgv = sdgrid.leafGridView();

  for (const auto& cell : elements(sdgv))
    {
      cell.geometry();
      cell.geometryInFather();
      const typename MDGrid::LeafGridView::template Codim<0>::Entity& mde = sdgrid.multiDomainEntity(cell);
      mde.geometry();
      mde.geometryInFather();

      for (const auto& intersection : intersections(sdgv,cell))
        {
          intersection.geometry();
          intersection.geometryInInside();
          if (intersection.neighbor())
            intersection.geometryInOutside();
          const typename MDGrid::LeafGridView::Intersection& is1 = sdgrid.multiDomainIntersection(intersection);
          is1.geometry();
          is1.geometryInInside();
          if (intersection.neighbor())
            is1.geometryInOutside();
        }
    }
}


int main(int argc, char** argv)
{
  try {

    Dune::MPIHelper::instance(argc,argv);

    {
      Dune::FieldVector<double,2> L(1.0);
      Dune::array<int,2> s = {{2, 2}};

      typedef Dune::YaspGrid<2> HostGrid;
      HostGrid hostgrid(L,s);
      run_test(hostgrid);
    }

#if HAVE_UG
    {
      typedef Dune::UGGrid<2> HostGrid;
      Dune::FieldVector<double,2> lower_left(0.0);
      Dune::FieldVector<double,2> upper_right(1.0);
      std::array<unsigned int,2> elements = {{8,8}};
      Dune::shared_ptr<HostGrid> gridptr(
        Dune::StructuredGridFactory<HostGrid>::createSimplexGrid(
          lower_left,
          upper_right,
          elements
          )
        );
      run_test(*gridptr);
    }
#endif

    return 0;

  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (std::exception &e){
    std::cerr << "std reported error: " << e.what() << std::endl;
    return 2;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 3;
  }
}
