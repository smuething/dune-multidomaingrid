MultiDomainGrid
===============

[dune-multidomaingrid][1] is a meta grid built on top of the [DUNE][2]
grid interface. It can be used to carve out subdomains from an underlying
host grid which are then available as fully featured DUNE grids in their
own right.


Version
-------

This is version 3.0-dev of MultiDomainGrid. It is compatible with the 3.0-git
development series of the DUNE core modules. The versioning scheme of MultiDomainGrid
has recently been synchronized to that of the core modules to make it easier for users
to find a compatible release for their DUNE distribution.

An overview of changes to the library can be found in the file
[RELEASE_NOTES.md][13].


Features
--------

* Support for arbitrary subdomain topologies, including non-contiguous domains.

* Full support for MPI-parallel computations on the subdomains.

* Policy-driven choice of different storage backends tailored to the application
  scenario. The library currently ships with three different implementations:

  * A version optimized for a small number of domains that may overlap in arbitrary
    ways and that maintains O(1) storage and runtime complexity irrespective of the
    amount of subdomain overlap. This backend is best suited to standard multi-physics
    applications that combine a small number of regions with different physics.

  * A version optimized for a very large number of subdomains, with the total number
    known at compile time. This implementation places a compile-time limit on the
    number of subdomains that a grid entity may belong to; several of its algorithms
    are O(log N) in the number of overlapping subdomains N.

  * A variant of the former with a run-time configurable number of total subdomains.
    This version traits added flexibility with a slightly less optimal storage scheme.
    Apart from the switch from a compile-time to a run-time subdomain limit, it is mostly
    identical to the second implementation in terms of algorithmic and storage complexity.

* The subdomain layout can be modified during a simulation, with an API that aids in
  transferring solution data similar to the grid adaptivity interface of a DUNE grid.

* DUNE grid API extensions to query the subdomain membership of a given entity as well
  as conversion methods between `MultiDomainGrid` entities and `SubDomainGrid` entities.

* Iterators for automatic extraction of subdomain - subdomain interfaces, both for a
  given pair of subdomains and for all defined subdomains at once (i.e. in a single host
  grid traversal).

If you have downloaded a release tarball, you can find the autogenerated Doxygen
API documentation in doc/doxygen/html. Otherwise, you can build this documentation
yourself by calling "make doc". Note that you need Doxygen and GraphViz available at
configure time to be able to build the documentation.

If you need help, please ask via [GitHub][1]. If you find bugs, you can also submit
them to the [bugtracker][3]. Even better, if you have managed to fix a problem, open
a [pull request][4] to get your patch merged into the library.


High-level interface for multi-domain simulations
-------------------------------------------------

While dune-multidomaingrid is a very useful tool in its own right and is used in
a standalone fashion by a number of people, it was originally designed as a building
block for extending the high-level DUNE-based PDE solver toolbox [PDELab][5] with
support for multi-physics and multi-domain simulations. This support is contained in
the [dune-multidomain][6] library, which is an add-on module for PDELab that extends
the latter with concepts for subproblems and couplings between those subproblems and
uses MultiDomainGrid to provide the spatial layout of those subproblems and couplings.


Dependencies
------------

dune-multidomaingrid depends on the following software packages:

* [DUNE core libraries][1] (dune-common, dune-geometry, dune-grid, dune-istl,
  dune-localfunctions) version 3.0-git, and their respective dependencies.

* The [Boost][7] C++ libraries, in particular Boost.MPL and Boost.Fusion.

* MultiDomain's compiler requirements differ slightly from the underlying DUNE
  libraries: We require at least GCC 4.7 in C++11 mode. MultiDomainGrid should
  also work with very recent versions of ICC (icpc 14.0.3+) and current
  clang (3.3+).


License
-------

The MultiDomainGrid library, headers and test programs are free open-source software,
dual-licensed under version 3 or later of the GNU Lesser General Public License
and version 2 of the GNU General Public License with a [special run-time exception][8].

See the file [COPYING.md][9] for full copying permissions.


Installation
------------

Short installation instructions can be found in the file [README.GIT.md][10].
For a full explanation of the DUNE installation process please read
the [installation notes][11] or the [build system HOWTO][12].


Links
-----

[1]:  http://github.com/smuething/dune-multidomaingrid
[2]:  http://dune-project.org
[3]:  https://github.com/smuething/dune-multidomaingrid/issues
[4]:  https://github.com/smuething/dune-multidomaingrid/pulls
[5]:  http://dune-project.org/pdelab/
[6]:  http://github.com/smuething/dune-multidomain
[7]:  http://boost.org
[8]:  http://gcc.gnu.org/onlinedocs/libstdc++/faq.html#faq.license
[9]:  COPYING.md
[10]: README.GIT.md
[11]: http://dune-project.org/doc/installation-notes.html
[12]: http://dune-project.org/doc/buildsystem/buildsystem.pdf
[13]: RELEASES_NOTES.md