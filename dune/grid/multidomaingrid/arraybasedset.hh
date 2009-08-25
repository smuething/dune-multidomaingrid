#ifndef DUNE_MULTIDOMAINGRID_ARRAYBASEDSET_HH
#define DUNE_MULTIDOMAINGRID_ARRAYBASEDSET_HH

#include <limits>
#include <cstddef>
#include <cstdint>
#include <dune/common/static_assert.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/exceptions.hh>
#include <boost/integer_traits.hpp>
#include <type_traits>
#include <cassert>
#include <strings.h>

namespace Dune {

namespace mdgrid {

template<typename SubDomainType, std::size_t capacity>
class ArrayBasedSet;


template<typename SubDomainType, std::size_t capacity>
bool setContains(const ArrayBasedSet<SubDomainType,capacity>& a,
                 const ArrayBasedSet<SubDomainType,capacity>& b);


template<typename SubDomainType, std::size_t capacity>
void setAdd(ArrayBasedSet<SubDomainType,capacity>& a,
            const ArrayBasedSet<SubDomainType,capacity>& b);


template<typename SubDomainType, std::size_t capacity>
class ArrayBasedSet {

  friend bool setContains<>(const ArrayBasedSet<SubDomainType,capacity>& a,
                        const ArrayBasedSet<SubDomainType,capacity>& b);

  friend void setAdd<>(ArrayBasedSet<SubDomainType,capacity>& a,
                   const ArrayBasedSet<SubDomainType,capacity>& b);

public:
  static const std::size_t maxSize = capacity;
  static const SubDomainType emptyTag = boost::integer_traits<SubDomainType>::const_max;
  typedef SubDomainType DomainType;
  typedef typename std::array<DomainType,maxSize>::iterator ArrayIterator;
  typedef typename std::array<DomainType,maxSize>::const_iterator Iterator;
  typedef ArrayBasedSet<SubDomainType,capacity> This;

  enum SetState {emptySet,simpleSet,multipleSet};

  Iterator begin() const {
    return _set.begin();
  }

  Iterator end() const {
    return _set.begin() + _size;
  }

  bool contains(DomainType domain) const {
    return std::binary_search(_set.begin(),_set.begin() + _size,domain);
  }

  template<typename Set>
  bool containsAll(const Set& set) const {
    return setContains(*this,set);
  }

  bool simple() const {
    return _size == 1;
  }

  bool empty() const {
    return _size == 0;
  }

  SetState state() const {
    switch (_size) {
    case 0:
      return emptySet;
    case 1:
      return simpleSet;
    default:
      return multipleSet;
    }
  }

  std::size_t size() const {
    return _size;
  }

  void clear() {
    _size = 0;
  }

  void add(DomainType domain) {
    assert(_size < maxSize);
    assert(!std::binary_search(_set.begin(),_set.begin()+_size,domain));
    _set[_size++] = domain;
    std::sort(_set.begin(),_set.begin()+_size);
  }

  void remove(DomainType domain) {
    ArrayIterator it = std::lower_bound(_set.begin(),_set.begin()+_size,domain);
    assert(*it == domain);
    *it = emptyTag;
    std::sort(_set.begin(),_set.end() + (_size--));
  }

  void set(DomainType domain) {
    _size = 1;
    _set[0] = domain;
  }

  template<typename Set>
  void addAll(const Set& set) {
    setAdd(*this,set);
  }

  int domainOffset(DomainType domain) const {
    Iterator it = std::lower_bound(_set.begin(),_set.begin()+_size,domain);
    assert(*it == domain);
    return it - _set.begin();
  }

  ArrayBasedSet() :
    _size(0)
  {}

private:
  std::size_t _size;
  std::array<DomainType,maxSize> _set;

};


template<typename SubDomainType, std::size_t capacity>
inline bool setContains(const ArrayBasedSet<SubDomainType,capacity>& a,
                        const ArrayBasedSet<SubDomainType,capacity>& b) {
  return std::includes(a._set.begin(),a._set.begin() + a._size,b._set.begin(),b._set.begin() + b._size);
}

template<typename SubDomainType, std::size_t capacity>
inline void setAdd(ArrayBasedSet<SubDomainType,capacity>& a,
                   const ArrayBasedSet<SubDomainType,capacity>& b) {
  std::array<SubDomainType,2*capacity> tmp;
  typename std::array<SubDomainType,2*capacity>::iterator it = std::set_union(a._set.begin(), a._set.begin() + a._size,
                                                                              b._set.begin(), b._set.begin() + b._size,
                                                                              tmp.begin());
  a._size = it - tmp.begin();
  assert(a._size <= capacity);
  std::copy(tmp.begin(),++it,a._set.begin());
}

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_ARRAYBASEDSET_HH
