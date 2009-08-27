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

template<typename SubDomainType>
class SingleValueSet;


template<typename SubDomainType>
bool setContains(const SingleValueSet<SubDomainType>& a,
                 const SingleValueSet<SubDomainType>& b);


template<typename SubDomainType>
void setAdd(SingleValueSet<SubDomainType>& a,
            const SingleValueSet<SubDomainType>& b);


template<typename SubDomainType>
class SingleValueSet {

  friend bool setContains<>(const SingleValueSet<SubDomainType>& a,
                        const SingleValueSet<SubDomainType>& b);

  friend void setAdd<>(SingleValueSet<SubDomainType>& a,
                   const SingleValueSet<SubDomainType>& b);

public:
  static const std::size_t maxSize = 1;
  static const SubDomainType emptyTag = boost::integer_traits<SubDomainType>::const_max;
  typedef SubDomainType DomainType;
  typedef const SubDomainType* Iterator;
  typedef SingleValueSet<SubDomainType> This;

  enum SetState {emptySet,simpleSet,multipleSet};

  Iterator begin() const {
    return &_set;
  }

  Iterator end() const {
    return &_set + (_set != emptyTag ? 1 : 0);
  }

  bool contains(DomainType domain) const {
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

  void add(DomainType domain) {
    assert(_set == emptyTag);
    _set = domain;
  }

  void remove(DomainType domain) {
    assert(_set == domain);
    _set = emptyTag;
  }

  void set(DomainType domain) {
    _set = domain;
  }

  template<typename Set>
  void addAll(const Set& set) {
    setAdd(*this,set);
  }

  int domainOffset(DomainType domain) const {
    assert(_set == domain);
    return 0;
  }

  SingleValueSet() :
    _set(emptyTag)
  {}

private:
  SubDomainType _set;

};


template<typename SubDomainType>
inline bool setContains(const SingleValueSet<SubDomainType>& a,
                        const SingleValueSet<SubDomainType>& b) {
  return a._set == b._set;
}

template<typename SubDomainType>
inline bool setAdd(SingleValueSet<SubDomainType>& a,
                   const SingleValueSet<SubDomainType>& b) {
  assert(a._set == SingleValueSet<SubDomainType>::emptyTag);
  a._set = b._set;
}

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SINGLEVALUESET_HH
