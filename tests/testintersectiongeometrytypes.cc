#include "config.h"

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/multidomaingrid.hh>
#include <iostream>
#include <cassert>

template<typename HostGrid>
void run_test(HostGrid& hostgrid)
{
  typedef Dune::MultiDomainGrid<HostGrid,Dune::mdgrid::ArrayBasedTraits<2,8,8> > MDGrid;
  MDGrid mdgrid(hostgrid,true);

  typedef typename MDGrid::LeafGridView MDGV;
  typedef typename MDGV::template Codim<0>::Iterator Iterator;
  typedef typename MDGrid::SubDomainIndexType SubDomainIndexType;

  MDGV mdgv = mdgrid.leafView();

  mdgrid.startSubDomainMarking();

  for(Iterator it = mdgv.template begin<0>(); it != mdgv.template end<0>(); ++it)
      mdgrid.addToSubDomain(0,*it);

  mdgrid.preUpdateSubDomains();
  mdgrid.updateSubDomains();
  mdgrid.postUpdateSubDomains();

  typedef typename MDGrid::SubDomainGrid SDGrid;
  const SDGrid& sdgrid = mdgrid.subDomain(0);
  typename SDGrid::LeafGridView sdgv = sdgrid.leafView();

  const typename SDGrid::LeafGridView::template Codim<0>::Iterator endit = sdgv.template end<0>();
  for (typename SDGrid::LeafGridView::template Codim<0>::Iterator it = sdgv.template begin<0>();
       it != endit;
       ++it)
    {
      const typename SDGrid::LeafGridView::IntersectionIterator endiit = sdgv.iend(*it);
      for (typename SDGrid::LeafGridView::IntersectionIterator iit = sdgv.ibegin(*it);
           iit != endiit;
           ++iit)
        {
          iit->geometry();
          iit->geometryInInside();
          if (iit->neighbor())
            iit->geometryInOutside();
          const typename MDGrid::LeafGridView::Intersection& is1 = sdgrid.multiDomainIntersection(*iit);
          is1.geometry();
          is1.geometryInInside();
          if (iit->neighbor())
            is1.geometryInOutside();
        }
    }
}


int main(int argc, char** argv)
{
  try {
    {
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> s(2);
      Dune::FieldVector<bool,2> p(false);
      int overlap = 0;

      typedef Dune::YaspGrid<2> HostGrid;
      HostGrid hostgrid(L,s,p,overlap);
      run_test(hostgrid);
    }

#if HAVE_UG
    {
      typedef Dune::UGGrid<2> HostGrid;
      Dune::shared_ptr<HostGrid> gridptr(Dune::GmshReader<HostGrid>::read("simple.msh",true));
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
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
