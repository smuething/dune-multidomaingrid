#ifndef DUNE_MULTIDOMAINGRID_UTILITY_HH
#define DUNE_MULTIDOMAINGRID_UTILITY_HH

#include <tuple>
#include <dune/geometry/type.hh>
#include <dune/common/iteratorfacades.hh>

namespace Dune {

namespace mdgrid {

namespace util {

struct GeometryTypeHash {

  std::size_t operator()(GeometryType gt) const {
    std::size_t hash = gt.dim() * 509;
    return gt.dim() < 2 ? hash : hash + static_cast<std::size_t>(gt.id());
  }

};

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

} // namespace util

} // namespace mdgrid

} // namespace Dune


#endif // DUNE_MULTIDOMAINGRID_UTILITY_HH
