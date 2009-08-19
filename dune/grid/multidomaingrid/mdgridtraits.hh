#ifndef DUNE_MULTIDOMAINGRID_MDGRIDTRAITS_HH
#define DUNE_MULTIDOMAINGRID_MDGRIDTRAITS_HH


#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/push_back.hpp>

namespace Dune {

namespace mdgrid {

namespace mpl = boost::mpl;

template<typename sequence, bool val, int codim>
struct makeBoolVectorHelper {

  typedef typename makeBoolVectorHelper<typename mpl::push_front<sequence,mpl::bool_<val> >::type,val,codim-1>::type type;

};

template<typename sequence, bool val>
struct makeBoolVectorHelper<sequence,val,0> {

  typedef typename mpl::push_front<sequence,mpl::bool_<val> >::type type;

};

template<int dim, bool val>
struct makeBoolVector {

  typedef typename makeBoolVectorHelper<mpl::vector_c<bool,val>,dim-1,val>::type type;

};

template<int dim, std::size_t subDomainsPerCell> //, typename supportedCodims = typename makeBoolVector<dim,true>::type >
struct MDGridTraits {

  typedef int SubDomainType;
  static const SubDomainType empty = -1;
  static const int dimension = dim;

  static const std::size_t maxSubDomainsPerCell = subDomainsPerCell;
  static const SubDomainType maxSubDomainIndex = maxSubDomainsPerCell - 1;

  //typedef supportedCodims supportedCodimensions;

  struct EmptyCodimBase {};

  template<int codim>
  struct Codim {
    static const bool supported = true;
    static const std::size_t maxSubDomainsPerEntity = maxSubDomainsPerCell;
    typedef Dune::mdgrid::IntegralTypeSubDomainSet<SubDomainType,maxSubDomainsPerEntity> SubDomainSet;
    typedef std::array<int,maxSubDomainsPerEntity> MultiIndexContainer;
    typedef std::array<int,maxSubDomainsPerEntity> SizeContainer;
  };
  /*
  template<int codim>
  struct Codim : public SelectType<mpl::at_c<supportedCodimensions,codim>::type::value,CodimBase<codim>,EmptyCodimBase>::Type {
    static const bool supported = mpl::at_c<supportedCodimensions,codim>::type::value;
  };
  */
};

} // namespace mdrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_MDGRIDTRAITS_HH
