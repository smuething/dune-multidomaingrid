#ifndef DUNE_MULTIDOMAINGRID_INDEXSET_HH
#define DUNE_MULTIDOMAINGRID_INDEXSET_HH

#include <unordered_map>
#include <vector>
#include <array>
#include <algorithm>
#include <type_traits>
#include <tuple>
#include <boost/scoped_ptr.hpp>
#include <dune/common/iteratorfacades.hh>

namespace Dune {

namespace multidomaingrid {

struct GeometryTypeHash {

  std::size_t operator()(GeometryType gt) const {
    std::size_t hash = gt.dim() * 509;
    return gt.dim() < 2 ? hash : hash + static_cast<std::size_t>(gt.basicType());
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
struct collect_elementwise_struct {

  T& result;

  collect_elementwise_struct(T& r) :
    result(r)
  {}

  void operator()(T& val) {
    typedef typename T::iterator Iterator;
    Iterator end = result.end();
    Iterator rit = result.begin();
    Iterator vit = val.begin();
    for (; rit != end; ++rit, ++vit)
      *rit += *vit;
  }
};

template<typename T>
collect_elementwise_struct<T> collect_elementwise(T& result) {
  return collect_elementwise_struct<T>(result);
}

template<typename HostGridView, typename SubDomainSetType>
class IndexSet {

public:

  typedef SubDomainSetType SubDomainSet;

  typedef typename HostGridView::IndexSet HostIndexSet;
  typedef typename HostGridView::template Codim<0>::Iterator HostEntityIterator;
  typedef typename HostGridView::template Codim<0>::Entity HostEntity;
  typedef typename HostIndexSet::IndexType IndexType;
  typedef typename SubDomainSet::DomainType DomainType;
  typedef typename HostGridView::Grid::ctype ctype;
  static const std::size_t maxSubDomains = SubDomainSet::maxSize;
  static const std::size_t dimension = HostGridView::dimension;

  struct MapEntry {
    SubDomainSet domains;
    IndexType index;
  };

  typedef std::array<IndexType,maxSubDomains> SizeContainer;
  typedef std::unordered_map<GeometryType,std::vector<MapEntry>,GeometryTypeHash> IndexMap;
  typedef std::unordered_map<GeometryType,SizeContainer,GeometryTypeHash> SizeMap;
  typedef std::array<IndexType,maxSubDomains> MultiIndexContainer;

  void update(bool full) {
    const HostIndexSet& his = _gridView.indexSet();
    if (full) {
      for (int codim = 0; codim <= dimension; ++codim) {
	indexMap(codim).swap(IndexMap());
	indexMap(codim).rehash(his.geomTypes(codim).size());
	sizeMap(codim).swap(SizeMap());
	sizeMap(codim).rehash(his.geomTypes(codim).size());
      }
    }
    for (int codim = 0; codim <= dimension; ++codim) {
      for (std::vector<GeometryType>::const_iterator it = his.geomTypes(codim).begin(); it != his.geomTypes(codim).end(); ++it) {
	indexMap(codim)[*it].resize(his.size(*it));
	sizeMap(codim)[*it].assign(0);
      }
    }
    _multiIndexMap.clear();
    HostEntityIterator end = _gridView.template end<0>();
    IndexMap& im = indexMap(0);
    SizeMap& sm = sizeMap(0);
    for (HostEntityIterator it  = _gridView.template begin<0>(); it != end; ++it) {
      const HostEntity& he = *it;
      const GeometryType hgt = he.type();
      IndexType hostIndex = his.index(he);
      MapEntry& me = im[hgt][hostIndex];
      updateMapEntry(me,sm[hgt]);
      markSubIndices(he,me.domains,his,GenericReferenceElements<ctype,dimension>::general(hgt),1);
    }
    updateSubIndices(1);
  }

  void updateMapEntry(MapEntry& me, SizeContainer& sizes) {
    switch (me.domains.state()) {
    case SubDomainSet::simpleSet:
      me.index = sizes[*me.domains.begin()]++;
      break;
    case SubDomainSet::multipleSet:
      me.index = _multiIndexMap.size();
      _multiIndexMap.push_back(MultiIndexContainer());
      MultiIndexContainer& mic = _multiIndexMap.back();
      for (typename SubDomainSet::Iterator it = me.domains.begin(); it != me.domains.end(); ++it) {
	mic[*it] = sizes[*it]++;
      }
    }
  }

  void markSubIndices(const HostEntity& e, const SubDomainSet& domains, const HostIndexSet& his, const GenericReferenceElement<ctype,dimension>& refEl, int codim) {
    const int size = refEl.size(codim);
    IndexMap& im = indexMap(codim);
    for (int i = 0; i < size; ++i) {
      IndexType hostIndex = his.subIndex(e,i,codim);
      GeometryType gt = refEl.type(i,codim);
      im[gt][hostIndex].domains.addAll(domains);
    }
    if (codim < dimension) {
      markSubIndices(e,domains,his,refEl,codim+1);
    }
  }

  void updateSubIndices(int codim) {
    const typename IndexMap::iterator end = indexMap(codim).end();
    for (typename IndexMap::iterator it = indexMap(codim).begin(); it != end; ++it) {
      const GeometryType gt = it->first;
      std::vector<MapEntry>& indices = it->second;
      SizeContainer& sizes = sizeMap(codim)[gt];
      const typename std::vector<MapEntry>::iterator end = indices.end();
      for (typename std::vector<MapEntry>::iterator iit = indices.begin(); iit != end; ++iit) {
	updateMapEntry(*iit,sizes);
      }
    }
    if (codim < dimension) {
      updateSubIndices(codim+1);
    }
  }

  SubDomainSet& subDomainSet(const HostEntity& e) {
    return indexMap(0)[e.type()][_gridView.indexSet().index(e)].domains;
  }

  IndexType indexForSubDomain(DomainType subDomain, const HostEntity& e) {
    GeometryType gt = e.type();
    IndexType hostIndex = _gridView.indexSet().index(e);
    const MapEntry& me = indexMap(0)[gt][hostIndex];
    assert(me.domains.contains(subDomain));
    if (me.domains.simple()) {
      return me.index;
    } else {
      return _multiIndexMap[me.index][subDomain];
    }
  }

  explicit IndexSet(HostGridView gv) :
    _gridView(gv)
  {
    for (int codim = 0; codim <= dimension; ++codim) {
      _indexMap[codim].reset(new IndexMap());
      _sizeMap[codim].reset(new SizeMap());
    }
  }

private:
  std::array<boost::scoped_ptr<IndexMap>,dimension+1> _indexMap;
  std::array<boost::scoped_ptr<SizeMap>,dimension+1> _sizeMap;
  std::vector<MultiIndexContainer> _multiIndexMap;
  HostGridView _gridView;

  IndexMap& indexMap(std::size_t codim) {
    return *(_indexMap[codim]);
  }

  SizeMap& sizeMap(std::size_t codim) {
    return *(_sizeMap[codim]);
  }

};

} // namespace multidomaingrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_INDEXSET_HH
