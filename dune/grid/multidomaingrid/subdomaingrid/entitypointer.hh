#ifndef DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_ENTITYPOINTER_HH
#define DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_ENTITYPOINTER_HH

namespace Dune {

namespace mdgrid {

template<int,PartitionIteratorType,typename>
class LeafIteratorWrapper;

template<int,PartitionIteratorType,typename>
class LevelIteratorWrapper;

namespace subdomain {

template<int codim, typename GridImp>
class EntityPointerWrapper
{

  static const int dim = GridImp::dimension;

  struct Invalid {};

public:

  typedef EntityPointerWrapper EntityPointerImp;

  static const int codimension = codim;

  typedef typename GridImp::Traits::template Codim<codim>::Entity Entity;
  typedef EntityPointerWrapper<codim,GridImp> Base;

  typedef typename GridImp::MDGridType::Traits::template Codim<codim>::EntityPointer MultiDomainEntityPointer;
  typedef typename GridImp::MDGridType::Traits::HierarchicIterator MultiDomainHierarchicIterator;

  EntityPointerWrapper(const GridImp& grid, const MultiDomainEntityPointer& multiDomainEntityPointer) :
    _entityWrapper(grid,multiDomainEntityPointer)
  {}

  template<PartitionIteratorType pitype>
  EntityPointerWrapper(const GridImp& grid,
                       const EntityIterator<codim,const typename GridImp::MDGridType,Dune::mdgrid::LeafIteratorWrapper<codim,pitype,const typename GridImp::MDGridType> >& multiDomainEntityPointer) :
    _entityWrapper(grid,multiDomainEntityPointer)
  {}

  template<PartitionType pitype>
  EntityPointerWrapper(const GridImp& grid,
                       const EntityIterator<codim,const typename GridImp::MDGridType,Dune::mdgrid::LevelIteratorWrapper<codim,pitype,const typename GridImp::MDGridType> >& multiDomainEntityPointer) :
    _entityWrapper(grid,multiDomainEntityPointer)
  {}

  EntityPointerWrapper(const GridImp& grid,
                       typename SelectType<codim==0,const MultiDomainHierarchicIterator&,Invalid>::Type multiDomainEntityPointer) :
    _entityWrapper(grid,multiDomainEntityPointer)
  {}

  EntityPointerWrapper(const EntityWrapper<codim,dim,GridImp>& entity) :
    _entityWrapper(entity._grid,entity._multiDomainEntityPointer)
  {}

  bool equals(const EntityPointerWrapper<codim,GridImp>& rhs) const {
    return _entityWrapper.multiDomainEntityPointer() == rhs._entityWrapper.multiDomainEntityPointer();
  }

  Entity& dereference() const {
    return _entityWrapper;
  }

  void compactify() {
    _entityWrapper.compactify();
  }

  int level() const {
    return _entityWrapper.level();
  }

protected:

  mutable MakeableEntityWrapper<codim,dim,GridImp> _entityWrapper;

};

} // namespace subdomain

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_SUBDOMAINGRID_ENTITYPOINTER_HH
