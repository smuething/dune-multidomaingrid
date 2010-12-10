#ifndef DUNE_MULTIDOMAINGRID_UTILITY_HH
#define DUNE_MULTIDOMAINGRID_UTILITY_HH

#include <tuple>
#include <dune/common/geometrytype.hh>
#include <dune/common/iteratorfacades.hh>
#include <boost/swap.hpp>
#include <boost/fusion/support/is_sequence.hpp>
#include <boost/fusion/view/zip_view.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/fusion/sequence/intrinsic/front.hpp>
#include <boost/fusion/sequence/intrinsic/back.hpp>
#include <boost/mpl/and.hpp>

namespace Dune {

namespace mdgrid {

namespace util {

namespace mpl = boost::mpl;
namespace fusion = boost::fusion;

template<typename Iterator, typename Predicate>
inline bool all_of(Iterator begin, Iterator end, Predicate pred) {
  for ( ; begin != end; ++begin) {
    if (!pred(*begin)) {
      return false;
    }
  }
  return true;
}

struct GeometryTypeHash {

  std::size_t operator()(GeometryType gt) const {
    std::size_t hash = gt.dim() * 509;
    return gt.dim() < 2 ? hash : hash + static_cast<std::size_t>(gt.id());
  }

};

template<typename T, std::size_t I>
class TupleElementPickingIterator : public ForwardIteratorFacade<TupleElementPickingIterator<T,I>,
								 typename std::tuple_element<I,typename T::value_type>::type,
								 typename std::add_lvalue_reference<typename std::tuple_element<I,typename T::value_type>::type>::type,
								 typename T::difference_type
								 > {

public:

  typedef ForwardIteratorFacade<TupleElementPickingIterator<T,I>,
				typename std::tuple_element<I,typename T::value_type>::type,
				typename std::add_lvalue_reference<typename std::tuple_element<I,typename T::value_type>::type>::type,
				typename T::difference_type
				> BaseType;

  typedef typename BaseType::Reference Reference;
  typedef typename BaseType::Value Value;
  typedef TupleElementPickingIterator<T,I> ThisType;

  Reference dereference() const {
    return std::get<I>(*_it);
  }

  bool equals(const ThisType& rhs) const {
    return _it == rhs._it;
  }

  void increment() {
    ++_it;
  }

private:

  T _it;

public:

  TupleElementPickingIterator(T it) :
    _it(it)
  {}

};

template<std::size_t I, typename T>
TupleElementPickingIterator<T,I> pick_element(T it) {
  return TupleElementPickingIterator<T,I>(it);
}

template<typename T>
TupleElementPickingIterator<T,1> value_iterator(T it) {
  return TupleElementPickingIterator<T,1>(it);
}

template<typename T>
TupleElementPickingIterator<T,0> key_iterator(T it) {
  return TupleElementPickingIterator<T,0>(it);
}

template<typename T, typename binary_function>
struct collect_elementwise_struct {

  T& result;
  binary_function func;

  collect_elementwise_struct(T& r, binary_function f = binary_function()) :
    result(r),
    func(f)
  {}

  void operator()(T& val) {
    std::transform(val.begin(),val.end(),result.begin(),result.begin(),func);
  }
};

template<typename binary_function, typename T>
collect_elementwise_struct<T,binary_function> collect_elementwise(T& result, binary_function f = binary_function()) {
  return collect_elementwise_struct<T,binary_function>(result,f);
}

namespace detail {

  struct swap
  {
    template<typename Elem>
    struct result
    {
      typedef void type;
    };

    template<typename Elem>
    void operator()(Elem const& e) const
    {
      front(e).swap(back(e));
    }
  };
}

template<typename Seq1, typename Seq2>
typename boost::enable_if<mpl::and_<fusion::traits::is_sequence<Seq1>, fusion::traits::is_sequence<Seq2> >, void>::type
swap(Seq1& lhs, Seq2& rhs)
{
  typedef fusion::vector<Seq1&, Seq2&> references;
  for_each(fusion::zip_view<references>(references(lhs, rhs)), detail::swap());
}

} // namespace util

} // namespace mdgrid

} // namespace Dune


#endif // DUNE_MULTIDOMAINGRID_UTILITY_HH
