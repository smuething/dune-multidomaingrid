#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <dune/common/unused.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/multidomaingrid.hh>
#include <iostream>
#include <cassert>

int main(int argc, char** argv)
{

  try {

    Dune::MPIHelper::instance(argc,argv);

    Dune::FieldVector<double,2> L(1.0);
    std::array<int,2> s = {{2, 2}};

    typedef Dune::YaspGrid<2> HostGrid;
    HostGrid hostgrid(L,s);

    //typedef Dune::MultiDomainGrid<HostGrid,Dune::mdgrid::FewSubDomainsTraits<2,8> > MDGrid;
    typedef Dune::MultiDomainGrid<HostGrid,Dune::mdgrid::ArrayBasedTraits<2,8,8> > MDGrid;
    MDGrid mdgrid(hostgrid,true);

    typedef MDGrid::LeafGridView MDGV;
    typedef MDGV::Codim<0>::Iterator Iterator;

    MDGV mdgv = mdgrid.leafGridView();

    mdgrid.startSubDomainMarking();

    for(Iterator it = mdgv.begin<0>(); it != mdgv.end<0>(); ++it)
      mdgrid.addToSubDomain(0,*it);

    mdgrid.preUpdateSubDomains();
    mdgrid.updateSubDomains();
    mdgrid.postUpdateSubDomains();

    typedef MDGrid::SubDomainGrid SDGrid;
    const SDGrid& sdgrid = mdgrid.subDomain(0);
    SDGrid::LeafGridView sdgv = sdgrid.leafGridView();

    SDGrid::LeafGridView::Codim<0>::Iterator it = sdgv.begin<0>();
    SDGrid::LeafGridView::IntersectionIterator iit = sdgv.ibegin(*it);

    const MDGrid::LeafGridView::Intersection& is1 DUNE_UNUSED = sdgrid.multiDomainIntersection(*iit);
    const MDGrid::LeafGridView::Intersection& is2 DUNE_UNUSED = mdgrid.multiDomainIntersection(*iit);

    assert(is1.geometry().center() == is2.geometry().center());

    return 0;

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

}
