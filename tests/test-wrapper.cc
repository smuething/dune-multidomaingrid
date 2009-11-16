// $Id$

#include <config.h>

#include <iostream>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/multidomaingrid.hh>

#include <dune/grid/test/gridcheck.cc>
#include <dune/grid/test/checkcommunicate.cc>
#include <dune/grid/test/checkgeometryinfather.cc>
#include <dune/grid/test/checkintersectionit.cc>

int rank;

template <int dim>
void check_yasp(bool p0=false) {
  typedef Dune::FieldVector<int,dim> iTupel;
  typedef Dune::FieldVector<double,dim> fTupel;
  typedef Dune::FieldVector<bool,dim> bTupel;

  std::cout << std::endl << "YaspGrid<" << dim << ">";
  if (p0) std::cout << " periodic\n";
  std::cout << std::endl << std::endl;

  fTupel Len; Len = 1.0;
  iTupel s; s = 3;
  bTupel p; p = false;
  p[0] = p0;
  int overlap = 1;

  //#if HAVE_MPI
  //Dune::YaspGrid<dim> grid(MPI_COMM_WORLD,Len,s,p,overlap);
  //#else
  Dune::YaspGrid<dim> wgrid(Len,s,p,overlap);

  typedef Dune::MultiDomainGrid<Dune::YaspGrid<dim>,Dune::mdgrid::FewSubDomainsTraits<dim,4> > MDGrid;

  MDGrid grid(wgrid);
  //#endif

  grid.globalRefine(2);

  gridcheck(grid);

  typedef typename MDGrid::template Codim<0>::Entity Entity;
  typedef typename MDGrid::LeafGridView::template Codim<0>::Iterator Iterator;

  typename MDGrid::LeafGridView gv = grid.leafView();

  grid.startSubDomainMarking();
  for (Iterator it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    const Entity& e = *it;
    //IndexSet::SubDomainSet& sds = is.subDomainSet(e);
    Dune::FieldVector<typename MDGrid::ctype,dim> c = e.geometry().global(Dune::GenericReferenceElements<typename MDGrid::ctype,dim>::general(e.type()).position(0,0));
    double x = c[0];
    double y = dim > 1 ? c[1] : 0.5;
    if (x > 0.2) {
      if (y > 0.3 && y < 0.7) {
        if (x < 0.8)
          grid.addToSubDomain(1,e);
        if (x > 0.6)
          grid.addToSubDomain(0,e);
      } else {
        grid.addToSubDomain(0,e);
      }
    }
  }
  grid.preUpdateSubDomains();
  grid.updateSubDomains();
  grid.postUpdateSubDomains();

  // check communication interface
  checkCommunication(grid,-1,Dune::dvverb);
  for(int l=0; l<=grid.maxLevel(); ++l)
    checkCommunication(grid,l,Dune::dvverb);

  // check the method geometryInFather()
  checkGeometryInFather(grid);
  // check the intersection iterator and the geometries it returns
  checkIntersectionIterator(grid);

  gridcheck(grid.subDomain(0));
  checkGeometryInFather(grid.subDomain(0));
  checkIntersectionIterator(grid.subDomain(0));

  gridcheck(grid.subDomain(1));
  checkGeometryInFather(grid.subDomain(1));
  checkIntersectionIterator(grid.subDomain(1));

};

int main (int argc , char **argv) {
  try {
#if HAVE_MPI
    // initialize MPI
    MPI_Init(&argc,&argv);

    // get own rank
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif

    check_yasp<1>();
    //check_yasp<1>(true);
    check_yasp<2>();
    //check_yasp<2>(true);
    check_yasp<3>();
    //check_yasp<3>(true);
    //check_yasp<4>();

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

#if HAVE_MPI
  // Terminate MPI
  MPI_Finalize();
#endif

  return 0;
};
