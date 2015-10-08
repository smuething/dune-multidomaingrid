MultiDomainGrid
===============

This is the 2.4.0 release of MultiDomainGrid, a meta grid build on top of the
[DUNE][1] grid interface. See the file README for further details including
license information.

If you need help, please ask via [GitHub][2]) If you find bugs, you can also
submit them to the [bugtracker][3]. Even better, if you have managed to fix a
problem, open a [pull request][4] to get your patch merged into the library.


Changes
-------


### MultiDomainGrid 2.4

* Bump compiler requirements to GCC 4.7.

* Compatiblity with DUNE 2.4.0.

* New implementation of per-codim data containers on top of std::tuple.

* Removed Boost dependency.

* New dependency on TypeTree for tuple iteration infrastructure.

* Copyable entities and intersections as required by dune-grid 2.4.0. Apart from
  the new API presented to the user, there are two important additional consequences:

  * The characteristics of the performance overhead of using MultiDomainGrid are different:
    Initial measurements show that using just a MultiDomainGrid is slightly slower than
    before (about 5-10%), but when stacking multiple meta grids is much less problematic.

  * The underlying host grid of a MultiDomainGrid must have been ported to the new copyable
    entity and intersection APIs. Unported grids will cause compilation errors. Apart from
    the old ALUGrid bindings (which are deprecated in DUNE 2.4), all important grids have
    been ported for the 2.4 release.

* Relesae history

  * `2.4.0-rc1` Initial release candidate.

### MultiDomainGrid 2.3

* First "official" release.

* New storage policy that supports specifying the maximum number of subdomains at runtime.

* Support for building with CMake.

* Compatibility with DUNE 2.3.x.

* Some minor bugfixes.

* Release history

  * `2.3.1` Final release, compatible with DUNE 2.3.1.

  * `2.3.0-rc2` Fix stupid build system problems
    * distribute documentation files
    * make sure there is a README file; otherwise `make dist` breaks
    * fix pkg-config info

  * `2.3.0-rc1` Initial release candidate


Caveats
-------

The following list is a non-exhaustive overview of possible problems you might
encounter with the current release.


### General

* Compile times can be really long for non-trivial problems. Some developers
  have had good success with using the clang compiler instead of GCC during
  development and bug-testing to reduce compile times.


Links
-----

[1]: http://dune-project.org
[2]: http://github.com/smuething/dune-multidomaingrid
[3]: https://github.com/smuething/dune-multidomaingrid/issues
[4]: https://github.com/smuething/dune-multidomaingrid/pulls
