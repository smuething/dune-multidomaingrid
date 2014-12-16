#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINSET_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINSET_HH

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <strings.h>
#include <dune/common/typetraits.hh>
#include <dune/common/iteratorfacades.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/multidomaingrid/utility.hh>

namespace Dune {

namespace mdgrid {

// forward declarations
template<typename SubDomainIndex, std::size_t capacity>
class IntegralTypeSubDomainSet;

template<typename SubDomainIndex, std::size_t capacity>
bool setContains(const IntegralTypeSubDomainSet<SubDomainIndex,capacity>& a,
                 const IntegralTypeSubDomainSet<SubDomainIndex,capacity>& b);

template<typename SubDomainIndex, std::size_t capacity>
void setAdd(IntegralTypeSubDomainSet<SubDomainIndex,capacity>& a,
            const IntegralTypeSubDomainSet<SubDomainIndex,capacity>& b);

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

    static_assert(std::numeric_limits<typename candidate::type>::is_specialized,"numeric_limits<> lacks specialization");

    typedef typename conditional<capacity <= static_cast<std::size_t>(std::numeric_limits<typename candidate::type>::digits),
                                 typename candidate::type,
                                 typename SetStorageTester<capacity,
                                                           typename candidate::next_candidate
                                                           >::type
                                 >::type type;

    static_assert((!std::is_same<type,void>::value),"unsupported maximum number of subdomains");

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
  inline std::size_t Log2<uint8_t>::doCalculation(uint8_t value) {
    return __builtin_ffs(value)-1;
  }

  //! \internal
  template<>
  inline std::size_t Log2<uint16_t>::doCalculation(uint16_t value) {
    return __builtin_ffs(value)-1;
  }

  //! \internal
  template<>
  inline std::size_t Log2<uint32_t>::doCalculation(uint32_t value) {
    return __builtin_ffsl(value)-1;
  }

  //! \internal
  template<>
  inline std::size_t Log2<uint64_t>::doCalculation(uint64_t value) {
    return __builtin_ffsll(value)-1;
  }

  //! \internal
  template<typename SubDomainIndex, typename SetStorage>
  class Iterator : public ForwardIteratorFacade<Iterator<SubDomainIndex,SetStorage>,
						SubDomainIndex,
						SubDomainIndex,
						std::ptrdiff_t> {

    template<typename,std::size_t>
    friend class ::Dune::mdgrid::IntegralTypeSubDomainSet;

  public:

    typedef Iterator<SubDomainIndex,SetStorage> ThisType;
    static const SetStorage base = 1;

    SubDomainIndex dereference() const {
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
    SubDomainIndex _value;

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
  typedef SubDomainIndexT SubDomainIndex;
  typedef sds_detail::Iterator<SubDomainIndex,SetStorage> Iterator;
  typedef IntegralTypeSubDomainSet<SubDomainIndex,capacity> This;

  struct DataHandle
  {
    typedef SetStorage DataType;

    static bool fixedsize(int dim, int codim)
    {
      return true;
    }

    static std::size_t size(const IntegralTypeSubDomainSet& sds)
    {
      return 1;
    }

    template<typename MessageBufferImp>
    static void gather(MessageBufferImp& buf, const IntegralTypeSubDomainSet& sds)
    {
      buf.write(sds._set);
    }

    template<typename MessageBufferImp>
    static void scatter(MessageBufferImp& buf, IntegralTypeSubDomainSet& sds, std::size_t n)
    {
      IntegralTypeSubDomainSet h;
      buf.read(h._set);
      sds.addAll(h);
    }

  };

  enum SetState {emptySet,simpleSet,multipleSet};

  Iterator begin() const {
    return Iterator(_set);
  }

  Iterator end() const {
    return Iterator();
  }

  bool contains(SubDomainIndex domain) const {
    assert(domain < maxSize);
    return (base << domain) & _set;
  }

  template<typename Set>
  bool containsAll(const Set& set) const {
    return setContains(*this,set);
  }

  void difference(const IntegralTypeSubDomainSet& minuend, const IntegralTypeSubDomainSet& subtrahend)
  {
    _set = minuend._set & ~subtrahend._set;
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

  void add(SubDomainIndex domain) {
    assert(domain < maxSize);
    _set |= base << domain;
  }

  void remove(SubDomainIndex domain) {
    assert(domain < maxSize);
    _set &= ~(base << domain);
  }

  void set(SubDomainIndex domain) {
    assert(domain < maxSize);
    _set = base << domain;
  }

  template<typename Set>
  void addAll(const Set& rhs) {
    setAdd(*this,rhs);
  }

  int domainOffset(SubDomainIndex domain) const {
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
  return std::all_of(b.begin(),b.end(),[&a](decltype(*(b.begin())) i) { return a.contains(i); });
}

template<typename A, typename B>
inline void setAdd(A& a, const B& b) {
  std::for_each(b.begin(),b.end(),[&a](decltype(*(b.begin())) i) { a.add(i); });
}


template<typename SubDomainIndex, std::size_t capacity>
inline bool setContains(const IntegralTypeSubDomainSet<SubDomainIndex,capacity>& a,
                        const IntegralTypeSubDomainSet<SubDomainIndex,capacity>& b) {
  return (a._set & b._set) == b._set;
}

template<typename SubDomainIndex, std::size_t capacity>
inline void setAdd(IntegralTypeSubDomainSet<SubDomainIndex,capacity>& a,
                   const IntegralTypeSubDomainSet<SubDomainIndex,capacity>& b) {
  a._set |= b._set;
}

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINSET_HH
