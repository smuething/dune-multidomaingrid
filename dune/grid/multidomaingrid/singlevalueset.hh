#ifndef DUNE_MULTIDOMAINGRID_SINGLEVALUESET_HH
#define DUNE_MULTIDOMAINGRID_SINGLEVALUESET_HH

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

template<typename SubDomainIndexType>
class SingleValueSet;


template<typename SubDomainIndexType>
bool setContains(const SingleValueSet<SubDomainIndexType>& a,
                 const SingleValueSet<SubDomainIndexType>& b);


template<typename SubDomainIndexType>
void setAdd(SingleValueSet<SubDomainIndexType>& a,
            const SingleValueSet<SubDomainIndexType>& b);


template<typename SubDomainIndexT>
class SingleValueSet {

  friend bool setContains<>(const SingleValueSet<SubDomainIndexT>& a,
                        const SingleValueSet<SubDomainIndexT>& b);

  friend void setAdd<>(SingleValueSet<SubDomainIndexT>& a,
                   const SingleValueSet<SubDomainIndexT>& b);

public:

  typedef SubDomainIndexT SubDomainIndexType;
  typedef SubDomainIndexType DomainType DUNE_DEPRECATED;

  static const std::size_t maxSize = 1;
  static const SubDomainIndexType emptyTag = boost::integer_traits<SubDomainIndexType>::const_max;
  typedef const SubDomainIndexType* Iterator;
  typedef SingleValueSet<SubDomainIndexType> This;

  enum SetState {emptySet,simpleSet,multipleSet};

  struct DataHandle
  {
    typedef SubDomainIndexType DataType;

    static bool fixedsize(int dim, int codim)
    {
      return true;
    }

    static std::size_t size(const SingleValueSet& sds)
    {
      return 1;
    }

    template<typename MessageBufferImp>
    static void gather(MessageBufferImp& buf, const SingleValueSet& sds)
    {
      buf.write(sds._set);
    }

    template<typename MessageBufferImp>
    static void scatter(MessageBufferImp& buf, SingleValueSet& sds, std::size_t n)
    {
      SubDomainIndexType i;
      buf.read(i);
      assert(sds._set == emptyTag || sds._set == i);
      sds._set = i;
    }

  };

  Iterator begin() const {
    return &_set;
  }

  Iterator end() const {
    return &_set + (_set != emptyTag ? 1 : 0);
  }

  bool contains(SubDomainIndexType domain) const {
    return _set == domain;
  }

  template<typename Set>
  bool containsAll(const Set& set) const {
    return setContains(*this,set);
  }

  bool simple() const {
    return _set != emptyTag;
  }

  bool empty() const {
    return _set == emptyTag;
  }

  SetState state() const {
    return _set == emptyTag ? emptySet : simpleSet;
  }

  std::size_t size() const {
    return _set == emptyTag ? 0 : 1;
  }

  void clear() {
    _set = emptyTag;
  }

  void add(SubDomainIndexType domain) {
    assert(_set == emptyTag || _set == domain);
    _set = domain;
  }

  void remove(SubDomainIndexType domain) {
    assert(_set == domain);
    _set = emptyTag;
  }

  void set(SubDomainIndexType domain) {
    _set = domain;
  }

  template<typename Set>
  void addAll(const Set& set) {
    setAdd(*this,set);
  }

  int domainOffset(SubDomainIndexType domain) const {
    assert(_set == domain);
    return 0;
  }

  SingleValueSet() :
    _set(emptyTag)
  {}

private:
  SubDomainIndexType _set;

};


template<typename SubDomainIndexType>
inline bool setContains(const SingleValueSet<SubDomainIndexType>& a,
                        const SingleValueSet<SubDomainIndexType>& b) {
  return a._set == b._set;
}

template<typename SubDomainIndexType>
inline void setAdd(SingleValueSet<SubDomainIndexType>& a,
                   const SingleValueSet<SubDomainIndexType>& b) {
  assert(a._set == SingleValueSet<SubDomainIndexType>::emptyTag || a._set == b._set);
  a._set = b._set;
}

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SINGLEVALUESET_HH
