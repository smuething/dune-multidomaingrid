#ifndef DUNE_MULTIDOMAINGRID_ARRAYBASEDSET_HH
#define DUNE_MULTIDOMAINGRID_ARRAYBASEDSET_HH

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <strings.h>

#include <dune/common/typetraits.hh>
#include <dune/common/exceptions.hh>

namespace Dune {

namespace mdgrid {

template<typename SI, std::size_t capacity>
class ArrayBasedSet;


template<typename SI, std::size_t capacity>
bool setContains(const ArrayBasedSet<SI,capacity>& a,
                 const ArrayBasedSet<SI,capacity>& b);


template<typename SI, std::size_t capacity>
void setAdd(ArrayBasedSet<SI,capacity>& a,
            const ArrayBasedSet<SI,capacity>& b);


template<typename SI, std::size_t capacity>
class ArrayBasedSet {

  friend bool setContains<>(const ArrayBasedSet<SI,capacity>& a,
                        const ArrayBasedSet<SI,capacity>& b);

  friend void setAdd<>(ArrayBasedSet<SI,capacity>& a,
                   const ArrayBasedSet<SI,capacity>& b);

public:
  typedef SI SubDomainIndex;
  typedef SubDomainIndex SubDomainIndexType DUNE_DEPRECATED_MSG("Use SubDomainIndex instead.");
  typedef SubDomainIndex DomainType DUNE_DEPRECATED_MSG("Use SubDomainIndex instead.");

  static const std::size_t maxSize = capacity;
  static const SubDomainIndex emptyTag = std::numeric_limits<SubDomainIndex>::max();
  typedef typename std::array<SubDomainIndex,maxSize>::iterator ArrayIterator;
  typedef typename std::array<SubDomainIndex,maxSize>::const_iterator Iterator;
  typedef ArrayBasedSet<SubDomainIndex,capacity> This;

  enum SetState {emptySet,simpleSet,multipleSet};

  struct DataHandle
  {
    typedef SubDomainIndex DataType;

    static bool fixedsize(int dim, int codim)
    {
      return false;
    }

    static std::size_t size(const ArrayBasedSet& sds)
    {
      return sds.size();
    }

    template<typename MessageBufferImp>
    static void gather(MessageBufferImp& buf, const ArrayBasedSet& sds)
    {
      for(Iterator it = sds.begin(); it != sds.end(); ++it)
        buf.write(*it);
    }

    template<typename MessageBufferImp>
    static void scatter(MessageBufferImp& buf, ArrayBasedSet& sds, std::size_t n)
    {
      ArrayBasedSet h;
      h._size = n;
      ArrayIterator end = h._set.begin() + n;
      for (ArrayIterator it = h._set.begin(); it != end; ++it)
        buf.read(*it);
      sds.addAll(h);
    }

  };

  Iterator begin() const {
    return _set.begin();
  }

  Iterator end() const {
    return _set.begin() + _size;
  }

  bool contains(SubDomainIndex domain) const {
    return std::binary_search(_set.begin(),_set.begin() + _size,domain);
  }

  template<typename Set>
  bool containsAll(const Set& set) const {
    return setContains(*this,set);
  }

  void difference(const ArrayBasedSet& minuend, const ArrayBasedSet& subtrahend)
  {
    Iterator res = std::set_difference(minuend.begin(),minuend.end(),
                                       subtrahend.begin(),subtrahend.end(),
                                       _set.begin());
    _size = res - _set.begin();
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

  void add(SubDomainIndex domain) {
    if (!std::binary_search(_set.begin(),_set.begin()+_size,domain)) {
      assert(_size < maxSize);
      _set[_size++] = domain;
      std::sort(_set.begin(),_set.begin()+_size);
    }
  }

  void remove(SubDomainIndex domain) {
    ArrayIterator it = std::lower_bound(_set.begin(),_set.begin()+_size,domain);
    assert(*it == domain);
    *it = emptyTag;
    std::sort(_set.begin(),_set.end() + (_size--));
  }

  void set(SubDomainIndex domain) {
    _size = 1;
    _set[0] = domain;
  }

  template<typename Set>
  void addAll(const Set& set) {
    setAdd(*this,set);
  }

  int domainOffset(SubDomainIndex domain) const {
    Iterator it = std::lower_bound(_set.begin(),_set.begin()+_size,domain);
    assert(*it == domain);
    return it - _set.begin();
  }

  ArrayBasedSet() :
    _size(0)
  {}

  bool operator==(const ArrayBasedSet& r) const {
    return _size == r._size && std::equal(_set.begin(),_set.begin()+_size,r._set.begin());
  }

  bool operator!=(const ArrayBasedSet& r) const {
    return !operator==(r);
  }

private:
  std::size_t _size;
  std::array<SubDomainIndex,maxSize> _set;

};


template<typename SubDomainIndex, std::size_t capacity>
inline bool setContains(const ArrayBasedSet<SubDomainIndex,capacity>& a,
                        const ArrayBasedSet<SubDomainIndex,capacity>& b) {
  return std::includes(a._set.begin(),a._set.begin() + a._size,b._set.begin(),b._set.begin() + b._size);
}

template<typename SubDomainIndex, std::size_t capacity>
inline void setAdd(ArrayBasedSet<SubDomainIndex,capacity>& a,
                   const ArrayBasedSet<SubDomainIndex,capacity>& b) {
  std::array<SubDomainIndex,2*capacity> tmp;
  typename std::array<SubDomainIndex,2*capacity>::iterator it = std::set_union(a._set.begin(), a._set.begin() + a._size,
                                                                                   b._set.begin(), b._set.begin() + b._size,
                                                                                   tmp.begin());
  a._size = it - tmp.begin();
  assert(a._size <= capacity);
  std::copy(tmp.begin(),++it,a._set.begin());
}

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_ARRAYBASEDSET_HH
