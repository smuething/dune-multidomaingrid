#ifndef DUNE_MULTIDOMAINGRID_MULTDIDOMAINMCMGMAPPER_HH
#define DUNE_MULTIDOMAINGRID_MULTDIDOMAINMCMGMAPPER_HH

#include <iostream>
#include <map>
#include <dune/grid/common/mcmgmapper.hh>

/**
 * @file
 * @brief  Mapper for multiple codim and multiple geometry types
 * @author Peter Bastian
 */

namespace Dune {

namespace mdgrid {

/**
 * @addtogroup Mapper
 *
 * @{
 */

/** @brief Implementation class for a multiple codim and multiple geometry type mapper.
 *
 * In this implementation of a mapper the entity set used as domain for the map consists
 * of the entities of a subset of codimensions in the given index set. The index
 * set may contain entities of several geometry types. This
 * version is usually not used directly but is used to implement versions for leafwise and levelwise
 * entity sets.
 *
 * Template parameters are:
 *
 * \par GV
 *    A Dune GridView type.
 * \par Layout
 *  A helper class with a method contains(), that returns true for all geometry
 *  types that are in the domain of the map.  The class should be of the following
 *  shape
 \code
 template<int dim>
 struct LayoutClass {
 bool contains (Dune::GeometryType gt) const {
 // Return true if gt is in the domain of the map
 }
 };
 \endcode
 *
 * If you don't want to use the default constructor of the LayoutClass you can construct it yourself
 * and hand it to the respective constructor.
 */
template <typename GV, template<int> class Layout>
class MultiDomainMCMGMapper : public MultipleCodimMultipleGeomTypeMapper<GV,Layout > {

  typedef MultipleCodimMultipleGeomTypeMapper<GV,Layout > Base;

public:

  typedef typename GV::IndexSet::IndexType IndexType;
  typedef typename GV::Grid::SubDomainIndexType SubDomainIndexType;
  typedef SubDomainIndexType SubDomainType DUNE_DEPRECATED;

  MultiDomainMCMGMapper (const GV& gridView, const Layout<GV::dimension> layout) :
    Base(gridView,layout),
    _is(gridView.indexSet()),
    _layout(layout)
  {
    update();
  }

  /** @brief Construct mapper from grid and one of its index sets.

      \param grid A Dune grid object.
      \param indexset IndexSet object returned by grid.

  */
  MultiDomainMCMGMapper (const GV& gridView) :
    Base(gridView),
    _is(gridView.indexSet())
  {
    update();
  }

  /** @brief Map entity to array index.

      \param e Reference to codim cc entity, where cc is the template parameter of the function.
      \return An index in the range 0 ... Max number of entities in set - 1.
  */
  template<class EntityType>
  int map (SubDomainIndexType subDomain, const EntityType& e) const
  {
    return _is.index(subDomain,e) + _offset[subDomain].find(e.type())->second;
  }

  /** @brief Map subentity of codim 0 entity to array index.

      \param e Reference to codim 0 entity.
      \param i Number of subentity of e
      \param codim Codimension of the subendity
      \return An index in the range 0 ... Max number of entities in set - 1.
  */
  int map (SubDomainIndexType subDomain, const typename GV::template Codim<0>::Entity& e, int i, unsigned int codim) const
  {
    GeometryType gt=GenericReferenceElements<double,GV::dimension>::general(e.type()).type(i,codim);
    return _is.subIndex(subDomain,e,i,codim) + _offset[subDomain].find(gt)->second;
  }

  /** @brief Return total number of entities in the entity set managed by the mapper.

      This number can be used to allocate a vector of data elements associated with the
      entities of the set. In the parallel case this number is per process (i.e. it
      may be different in different processes).

      \return Size of the entity set.
  */
  int size (SubDomainIndexType subDomain) const
  {
    return _n[subDomain];
  }

  /** @brief Returns true if the entity is contained in the index set

      \param e Reference to entity
      \param result integer reference where corresponding index is  stored if true
      \return true if entity is in entity set of the mapper
  */
  template<class EntityType>
  bool contains (SubDomainIndexType subDomain, const EntityType& e, IndexType& result) const
  {
    if(!_is.contains(subDomain,e) || !_layout.contains(e.type()))
      {
        result = 0;
        return false;
      }
    result = map(subDomain,e);
    return true;
  }

  /** @brief Returns true if the entity is contained in the index set

      \param e Reference to codim 0 entity
      \param i subentity number
      \param result integer reference where corresponding index is  stored if true
      \return true if entity is in entity set of the mapper
  */
  template<int cc> // this is now the subentity's codim
  bool contains (SubDomainIndexType subDomain, const typename GV::template Codim<0>::Entity& e, int i, IndexType& result) const
  {
    result = this->template map<cc>(subDomain,e,i);
    return true;
  }


  /** @brief Recalculates map after mesh adaptation
   */
  void update()
  {
    static_cast<Base*>(this)->update();
    for (SubDomainIndexType subDomain = 0; subDomain < GV::Grid::maxSubDomainIndex; ++subDomain) {
      std::size_t& n = _n[subDomain];
      std::map<GeometryType,IndexType>& offset = _offset[subDomain];
      n=0; // zero data elements
      offset.clear();

      // Compute offsets for the different geometry types.
      // Note that mapper becomes invalid when the grid is modified.
      for (int c=0; c<=GV::dimension; c++)
        for (size_t i=0; i<_is.geomTypes(subDomain,c).size(); i++)
          if (_layout.contains(_is.geomTypes(subDomain,c)[i]))
            {
              offset[_is.geomTypes(subDomain,c)[i]] = n;
              n += _is.size(subDomain,_is.geomTypes(subDomain,c)[i]);
            }
    }
  }

private:
  std::array<std::size_t,GV::Grid::maxSubDomainIndex> _n;
  const typename GV::IndexSet& _is;
  std::array<std::map<GeometryType,IndexType>,GV::Grid::maxSubDomainIndex> _offset; // provide a map with all geometry types
  mutable Layout<GV::dimension> _layout; // get layout object
};

#if 0

/** @brief Multiple codim and multiple geometry type mapper for leaf entities.

    This mapper uses all leaf entities of a certain codimension as its entity set.

    Template parameters are:

    \par G
    A %Dune grid type.
    \par Layout
    A helper class with a method contains(), that returns true for all geometry
    types that are in the domain of the map.  The class should be of the following
    shape
    \code
    template<int dim>
    struct LayoutClass {
    bool contains (Dune::GeometryType gt) {
    // Return true if gt is in the domain of the map
    }
    };
    \endcode
*/

template <typename G, template<int> class Layout>
class LeafMultipleCodimMultipleGeomTypeMapper
  : public MultipleCodimMultipleGeomTypeMapper<typename G::LeafGridView,Layout>
{
public:
  /** @brief The constructor
      @param grid A reference to a grid.
  */
  LeafMultipleCodimMultipleGeomTypeMapper (const G& grid)
    : MultipleCodimMultipleGeomTypeMapper<typename G::LeafGridView,Layout>(grid.leafView())
  {}

  /** @brief The constructor
   *
   * Use this constructor to provide a custom layout object e.g. not
   * using the default constructor.
   *
   * @param grid A reference to a grid.
   * @param layout A layout object
   */
  LeafMultipleCodimMultipleGeomTypeMapper (const G& grid, const Layout<G::dimension> layout)
    : MultipleCodimMultipleGeomTypeMapper<typename G::Traits::LeafIndexSet,Layout>(grid,grid.leafIndexSet(),layout)
  {}

};

/** @brief Multiple codim and multiple geometry type mapper for entities of one level.


    This mapper uses all entities of a certain codimension on a given level as its entity set.

    Template parameters are:

    \par G
    A %Dune grid type.
    \par Layout
    A helper class with a method contains(), that returns true for all geometry
    types that are in the domain of the map.  The class should be of the following
    shape
    \code
    template<int dim>
    struct LayoutClass {
    bool contains (Dune::GeometryType gt) {
    // Return true if gt is in the domain of the map
    }
    };
    \endcode
*/
template <typename G, template<int> class Layout>
class LevelMultipleCodimMultipleGeomTypeMapper
  : public MultipleCodimMultipleGeomTypeMapper<typename G::LevelGridView,Layout> {
public:
  /** @brief The constructor
      @param grid A reference to a grid.
      @param level A valid level of the grid.
  */
  LevelMultipleCodimMultipleGeomTypeMapper (const G& grid, int level)
    : MultipleCodimMultipleGeomTypeMapper<typename G::LevelGridView,Layout>(grid.levelView(level))
  {}

  /** @brief The constructor
   *
   * Use this constructor to provide a custom layout object e.g. not
   * using the default constructor.
   *
   * @param grid A reference to a grid.
   * @param layout A layout object
   */
  LevelMultipleCodimMultipleGeomTypeMapper (const G& grid, int level, const Layout<G::dimension> layout)
    : MultipleCodimMultipleGeomTypeMapper<typename G::Traits::LevekIndexSet,Layout>(grid,grid.levelIndexSet(level),layout)
  {}

};

#endif

/** @} */

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_MULTDIDOMAINMCMGMAPPER_HH
