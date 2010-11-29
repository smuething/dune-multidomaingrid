#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINSET_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINSET_HH

#include <limits>
#include <cstddef>
#include <cstdint>
#include <dune/common/static_assert.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/iteratorfacades.hh>
#include <dune/common/exceptions.hh>
#include <algorithm>
#include <boost/bind.hpp>
#include <cassert>
#include <strings.h>
#include <dune/grid/multidomaingrid/utility.hh>

namespace Dune {

namespace mdgrid {

// forward declarations
template<typename SubDomainIndexType, std::size_t capacity>
class IntegralTypeSubDomainSet;

template<typename SubDomainIndexType, std::size_t capacity>
bool setContains(const IntegralTypeSubDomainSet<SubDomainIndexType,capacity>& a,
                 const IntegralTypeSubDomainSet<SubDomainIndexType,capacity>& b);

template<typename SubDomainIndexType, std::size_t capacity>
void setAdd(IntegralTypeSubDomainSet<SubDomainIndexType,capacity>& a,
            const IntegralTypeSubDomainSet<SubDomainIndexType,capacity>& b);

//! @cond DEV_DOC

// \internal
namespace sds_detail {

  //! \internal
  template<typename T>
  struct Candidate;

  //! \internal
  template<>
  struct Candidate<uint8_t> {
    typedef uint8_t type;
    typedef Candidate<uint16_t> next_candidate;
  };

  //! \internal
  template<>
  struct Candidate<uint16_t> {
    typedef uint16_t type;
    typedef Candidate<uint32_t> next_candidate;
  };

  //! \internal
  template<>
  struct Candidate<uint32_t> {
    typedef uint32_t type;
    typedef Candidate<uint64_t> next_candidate;
  };

  //! \internal
  template<>
  struct Candidate<uint64_t> {
    typedef uint64_t type;
    typedef void next_candidate;
  };

  //! \internal
  template<std::size_t capacity, typename candidate>
  struct SetStorageTester {

    dune_static_assert(std::numeric_limits<typename candidate::type>::is_specialized,"numeric_limits<> lacks specialization");

    typedef typename SelectType<capacity <= static_cast<std::size_t>(std::numeric_limits<typename candidate::type>::digits),
					    typename candidate::type,
					    typename SetStorageTester<capacity,
								      typename candidate::next_candidate
								      >::type
					    >::Type type;

    dune_static_assert((!std::is_same<type,void>::value),"unsupported maximum number of subdomains");

  };

  /**
   * \internal
   */
  template<std::size_t capacity>
  struct SetStorageTester<capacity,void> {

    typedef void type;

  };

  //! \internal
  template<std::size_t capacity>
  struct SetStorageChooser {
    typedef typename SetStorageTester<capacity,Candidate<uint8_t> >::type type;
  };

  //! \internal
  template<typename T>
  class Log2 {

  private:
    inline static std::size_t doCalculation(T value);

  public:

    inline static std::size_t calculate(T value) {
      assert(value != 0 && (value & (value-1)) == 0); // make sure value is a power of 2
      return doCalculation(value);
    }
  };

  //! \internal asdf
  template<>
  std::size_t Log2<uint8_t>::doCalculation(uint8_t value) {
    return __builtin_ffs(value)-1;
  }

  //! \internal
  template<>
  std::size_t Log2<uint16_t>::doCalculation(uint16_t value) {
    return __builtin_ffs(value)-1;
  }

  //! \internal
  template<>
  std::size_t Log2<uint32_t>::doCalculation(uint32_t value) {
    return __builtin_ffsl(value)-1;
  }

  //! \internal
  template<>
  std::size_t Log2<uint64_t>::doCalculation(uint64_t value) {
    return __builtin_ffsll(value)-1;
  }

  //! \internal
  template<typename SubDomainIndexType, typename SetStorage>
  class Iterator : public ForwardIteratorFacade<Iterator<SubDomainIndexType,SetStorage>,
						SubDomainIndexType,
						SubDomainIndexType,
						std::ptrdiff_t> {

    template<std::size_t capacity>
    friend class ::Dune::mdgrid::IntegralTypeSubDomainSet;

  public:

    typedef Iterator<SubDomainIndexType,SetStorage> ThisType;
    static const SetStorage base = 1;

    SubDomainIndexType dereference() const {
      assert(_state != 0);
      return _value;
    }

    bool equals(const ThisType& rhs) const {
      return _state == rhs._state;
    }

    void increment() {
      _state &= ~(base << _value);
      if (_state != 0) {
	findNextValue();
      }
    }

  private:

    void findNextValue() {
      SetStorage lowestBit = _state & ((~_state) + 1);
      _value = Log2<SetStorage>::calculate(lowestBit);
    }

    explicit Iterator(SetStorage state) :
      _state(state),
      _value(0)
    {
      if (_state != 0)
	findNextValue();
    }

    // optimized constructor for end iterator
    explicit Iterator() :
      _state(0),
      _value(0)
    {}

    SetStorage _state;
    SubDomainIndexType _value;

  };

}

//! @endcond

template<typename SubDomainIndexT, std::size_t capacity>
class IntegralTypeSubDomainSet {

  friend bool setContains<>(const IntegralTypeSubDomainSet<SubDomainIndexT,capacity>& a,
                          const IntegralTypeSubDomainSet<SubDomainIndexT,capacity>& b);

  friend void setAdd<>(IntegralTypeSubDomainSet<SubDomainIndexT,capacity>& a,
                     const IntegralTypeSubDomainSet<SubDomainIndexT,capacity>& b);

  typedef typename sds_detail::SetStorageChooser<capacity>::type SetStorage;
  static const SetStorage base = 1;

public:
  static const std::size_t maxSize = capacity;
  typedef SubDomainIndexT SubDomainIndexType;
  typedef SubDomainIndexType DomainType DUNE_DEPRECATED;
  typedef sds_detail::Iterator<SubDomainIndexType,SetStorage> Iterator;
  typedef IntegralTypeSubDomainSet<SubDomainIndexType,capacity> This;

  enum SetState {emptySet,simpleSet,multipleSet};

  Iterator begin() const {
    return Iterator(_set);
  }

  Iterator end() const {
    return Iterator();
  }

  bool contains(SubDomainIndexType domain) const {
    assert(domain < maxSize);
    return (base << domain) & _set;
  }

  template<typename Set>
  bool containsAll(const Set& set) const {
    return setContains(*this,set);
  }

  bool simple() const {
    return (!empty()) && ((_set & (_set - 1)) == 0);
  }

  bool empty() const {
    return _set == 0;
  }

  SetState state() const {
    return (_set == 0 ? emptySet : (_set & (_set - 1)) == 0 ? simpleSet : multipleSet);
  }

  std::size_t size() const {
    std::size_t c = 0;
    for (SetStorage t = _set; t; ++c) {
      t &= t - 1;
    }
    return c;
  }

  void clear() {
    _set = 0;
  }

  void add(SubDomainIndexType domain) {
    assert(domain < maxSize);
    _set |= base << domain;
  }

  void remove(SubDomainIndexType domain) {
    assert(domain < maxSize);
    _set &= ~(base << domain);
  }

  void set(SubDomainIndexType domain) {
    assert(domain < maxSize);
    _set = base << domain;
  }

  template<typename Set>
  void addAll(const Set& rhs) {
    setAdd(*this,rhs);
  }

  int domainOffset(SubDomainIndexType domain) const {
    assert(domain >= 0 && domain < maxSize);
    return domain;
  }

  IntegralTypeSubDomainSet() :
    _set(0)
  {}

  bool operator==(const IntegralTypeSubDomainSet& r) const {
    return _set == r._set;
  }

  bool operator!=(const IntegralTypeSubDomainSet& r) const {
    return !operator==(r);
  }

private:
  SetStorage _set;

};


template<typename A, typename B>
inline bool setContains(const A& a, const B& b) {
  return util::all_of(b.begin(),b.end(),boost::bind(&A::contains,boost::ref(a),_1));
}

template<typename A, typename B>
inline void setAdd(A& a, const B& b) {
  std::for_each(b.begin(),b.end(),boost::bind(&A::add,boost::ref(a),_1));
}


template<typename SubDomainIndexType, std::size_t capacity>
inline bool setContains(const IntegralTypeSubDomainSet<SubDomainIndexType,capacity>& a,
                        const IntegralTypeSubDomainSet<SubDomainIndexType,capacity>& b) {
  return (a._set & b._set) == b._set;
}

template<typename SubDomainIndexType, std::size_t capacity>
inline void setAdd(IntegralTypeSubDomainSet<SubDomainIndexType,capacity>& a,
                   const IntegralTypeSubDomainSet<SubDomainIndexType,capacity>& b) {
  a._set |= b._set;
}

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINSET_HH
