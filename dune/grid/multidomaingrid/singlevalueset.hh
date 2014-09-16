#ifndef DUNE_MULTIDOMAINGRID_SINGLEVALUESET_HH
#define DUNE_MULTIDOMAINGRID_SINGLEVALUESET_HH

#include <limits>
#include <cstddef>
#include <cstdint>
#include <dune/common/typetraits.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/deprecated.hh>
#include <type_traits>
#include <cassert>
#include <strings.h>

namespace Dune {

namespace mdgrid {

template<typename SI>
class SingleValueSet;


template<typename SI>
bool setContains(const SingleValueSet<SI>& a,
                 const SingleValueSet<SI>& b);


template<typename SI>
void setAdd(SingleValueSet<SI>& a,
            const SingleValueSet<SI>& b);


template<typename SI>
class SingleValueSet {

  friend bool setContains<>(const SingleValueSet<SI>& a,
                        const SingleValueSet<SI>& b);

  friend void setAdd<>(SingleValueSet<SI>& a,
                   const SingleValueSet<SI>& b);

public:

  typedef SI SubDomainIndex;
  typedef SI SubDomainIndexType DUNE_DEPRECATED_MSG("Use SubDomainIndex instead.");
  typedef SI DomainType DUNE_DEPRECATED_MSG("Use SubDomainIndex instead.");

  static const std::size_t maxSize = 1;
  static const SubDomainIndex emptyTag = ~SubDomainIndex(0);
  typedef const SubDomainIndex* Iterator;
  typedef SingleValueSet<SubDomainIndex> This;

  enum SetState {emptySet,simpleSet,multipleSet};

  struct DataHandle
  {
    typedef SubDomainIndex DataType;

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
      SubDomainIndex i;
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

  bool contains(SubDomainIndex domain) const {
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

  void add(SubDomainIndex domain) {
    assert(_set == emptyTag || _set == domain);
    _set = domain;
  }

  void remove(SubDomainIndex domain) {
    assert(_set == domain);
    _set = emptyTag;
  }

  void set(SubDomainIndex domain) {
    _set = domain;
  }

  template<typename Set>
  void addAll(const Set& set) {
    setAdd(*this,set);
  }

  int domainOffset(SubDomainIndex domain) const {
    assert(_set == domain);
    return 0;
  }

  SingleValueSet() :
    _set(emptyTag)
  {}

private:
  SubDomainIndex _set;

};


template<typename SI>
inline bool setContains(const SingleValueSet<SI>& a,
                        const SingleValueSet<SI>& b) {
  return a._set == b._set;
}

template<typename SI>
inline void setAdd(SingleValueSet<SI>& a,
                   const SingleValueSet<SI>& b) {
  assert(a._set == SingleValueSet<SI>::emptyTag || a._set == b._set);
  a._set = b._set;
}

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SINGLEVALUESET_HH
