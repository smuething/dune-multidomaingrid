#ifndef DUNE_MULTIDOMAINGRID_MDGRIDTRAITS_HH
#define DUNE_MULTIDOMAINGRID_MDGRIDTRAITS_HH

#include <vector>

#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/integer_traits.hpp>

#include <dune/common/deprecated.hh>

#include <dune/grid/multidomaingrid/subdomainset.hh>
#include <dune/grid/multidomaingrid/arraybasedset.hh>
#include <dune/grid/multidomaingrid/singlevalueset.hh>

namespace Dune {

namespace mdgrid {

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

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

  typedef typename fusion::result_of::as_vector<typename makeBoolVectorHelper<mpl::vector_c<bool,val>,dim-1,val>::type>::type type;

};

template<int dim, int codim>
struct AllCodims {
  static const bool supported = true;
};

template<int dim, int codim>
struct CellAndVertexCodims {
  static const bool supported = (codim == 0 || codim == dim);
};

template<int dim, std::size_t subDomainsPerCell, std::size_t subDomainCount, template<int dim_, int codim> class supportedCodims = AllCodims>
struct ArrayBasedTraits {

  typedef int SubDomainIndex;
  typedef SubDomainIndex SubDomainIndexType DUNE_DEPRECATED_MSG("Use SubDomainIndex instead.");
  typedef SubDomainIndex SubDomainType DUNE_DEPRECATED_MSG("Use SubDomainIndex instead.");
  static const SubDomainIndex empty = -1;
  static const int dimension = dim;

  static const std::size_t maxSubDomainsPerCell = subDomainsPerCell;

  static constexpr bool maxSubDomainIndexIsStatic()
  {
    return true;
  }

  static constexpr SubDomainIndex maxSubDomainIndex()
  {
    return subDomainCount;
  }

  struct EmptyCodimBase {
    typedef int SizeContainer;
    typedef int MultiIndexContainer;
    typedef int SubDomainSet;
  };

  template<int codim>
  struct CodimBase {
    static const std::size_t maxSubDomainsPerEntity = (2<<(codim)) * maxSubDomainsPerCell;
    typedef Dune::mdgrid::ArrayBasedSet<SubDomainIndex,maxSubDomainsPerEntity> SubDomainSet;
    typedef std::array<int,maxSubDomainsPerEntity> MultiIndexContainer; // TODO: really int??
    typedef std::array<int,subDomainCount> SizeContainer; // TODO: really int??
  };

  template<int codim>
  struct Codim : public conditional<supportedCodims<dim,codim>::supported,CodimBase<codim>,EmptyCodimBase>::type {
    static const bool supported = supportedCodims<dim,codim>::supported;
  };

  template<int codim, typename SizeContainer>
  void setupSizeContainer(SizeContainer&) const
  {}

};


template<int dim, std::size_t subDomainsPerCell, template<int dim_, int codim> class supportedCodims = AllCodims>
struct DynamicSubDomainCountTraits {

  typedef int SubDomainIndex;
  typedef SubDomainIndex SubDomainIndexType DUNE_DEPRECATED_MSG("Use SubDomainIndex instead.");
  typedef SubDomainIndex SubDomainType DUNE_DEPRECATED_MSG("Use SubDomainIndex instead.");
  static const SubDomainIndex empty = -1;
  static const int dimension = dim;

  static constexpr bool maxSubDomainIndexIsStatic()
  {
    return false;
  }

  static const std::size_t maxSubDomainsPerCell = subDomainsPerCell;

  SubDomainIndex maxSubDomainIndex() const
  {
    return _subDomainCount;
  }

  struct EmptyCodimBase {
    typedef int SizeContainer;
    typedef int MultiIndexContainer;
    typedef int SubDomainSet;

    template<typename SC>
    static void setupSizeContainer(const SC&, std::size_t)
    {}

  };

  template<int codim>
  struct CodimBase {
    static const std::size_t maxSubDomainsPerEntity = (2<<(codim)) * maxSubDomainsPerCell;
    typedef Dune::mdgrid::ArrayBasedSet<SubDomainIndex,maxSubDomainsPerEntity> SubDomainSet;
    typedef std::array<int,maxSubDomainsPerEntity> MultiIndexContainer; // TODO: really int??
    typedef std::vector<int> SizeContainer; // TODO: really int??

    static void setupSizeContainer(SizeContainer& container, std::size_t subDomainCount)
    {
      container.resize(subDomainCount);
    }

  };

  template<int codim>
  struct Codim : public conditional<supportedCodims<dim,codim>::supported,CodimBase<codim>,EmptyCodimBase>::type {
    static const bool supported = supportedCodims<dim,codim>::supported;
  };

  DynamicSubDomainCountTraits(std::size_t subDomainCount)
    : _subDomainCount(subDomainCount)
  {}

  template<int codim, typename SizeContainer>
  void setupSizeContainer(SizeContainer& container) const
  {
    Codim<codim>::setupSizeContainer(container,_subDomainCount);
  }

private:

  const std::size_t _subDomainCount;

};


template<int dim, std::size_t maxSubDomains, template<int dim_, int codim> class supportedCodims = AllCodims >
struct FewSubDomainsTraits {

  typedef unsigned int SubDomainIndex;
  typedef SubDomainIndex SubDomainIndexType DUNE_DEPRECATED_MSG("Use SubDomainIndex instead.");
  typedef SubDomainIndex SubDomainType DUNE_DEPRECATED_MSG("Use SubDomainIndex instead.");
  static const SubDomainIndex empty = ~SubDomainIndex(0); // this is not used, but has to be present to make the compiler happy
  static const int dimension = dim;

  static const std::size_t maxSubDomainsPerCell = maxSubDomains;

  static constexpr bool maxSubDomainIndexIsStatic()
  {
    return true;
  }

  static constexpr SubDomainIndex maxSubDomainIndex()
  {
    return maxSubDomains - 1;
  }

  struct EmptyCodimBase {
    typedef int SizeContainer;
    typedef int MultiIndexContainer;
    typedef int SubDomainSet;
  };

  template<int codim>
  struct CodimBase {
    static const std::size_t maxSubDomainsPerEntity = maxSubDomains;
    typedef Dune::mdgrid::IntegralTypeSubDomainSet<SubDomainIndex,maxSubDomainsPerEntity> SubDomainSet;
    typedef std::array<int,maxSubDomainsPerEntity> MultiIndexContainer;
    typedef std::array<int,maxSubDomains> SizeContainer;
  };

  template<int codim>
  struct Codim : public conditional<supportedCodims<dim,codim>::supported,CodimBase<codim>,EmptyCodimBase>::type {
    static const bool supported = supportedCodims<dim,codim>::supported;
  };

  template<int codim>
  void setupSizeContainer(typename Codim<codim>::SizeContainer&) const
  {}

};

} // namespace mdrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_MDGRIDTRAITS_HH
