MultiDomainGrid
===============

This is the 2.3 release of MultiDomainGrid, a meta grid build on top of the
[DUNE][1] grid interface. See the file README for further details including
license information.

MultiDomainGrid 2.3 is a minor release, mainly focusing on bug fixes and keeping
current with changes to the underlying DUNE modules.

If you need help, please ask via [GitHub][2]) If you find bugs, you can also
submit them to the [bugtracker][3]. Even better, if you have managed to fix a
problem, open a [pull request][4] to get your patch merged into the library.


Changes
-------

### MultiDomainGrid 2.3

* First "official" release.

* New storage policy that supports specifying the maximum number of subdomains at runtime.

* Support for building with CMake.

* Compatibility with DUNE 2.3.x.

* Some minor bugfixes.

* Release history

  * `2.3.0-rc1` Initial release candidate


Caveats
-------

The following list is a non-exhaustive overview of possible problems you might
encounter with the current release.


### General

* Compile times can be really long for non-trivial problems. Some developers
  have had good success with using the clang compiler instead of GCC during
  development and bug-testing to reduce compile times.

* After MultiDomainGrid 2.3, the minimum compiler requirement of MultiDomainGrid
  will be increased to GCC 4.7. Please be aware of this change in minimum requirements.


Links
-----

[1]: http://dune-project.org
[2]: http://github.com/smuething/dune-multidomaingrid
[3]: https://github.com/smuething/dune-multidomaingrid/issues
[4]: https://github.com/smuething/dune-multidomaingrid/pulls