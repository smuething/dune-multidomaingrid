Preparing the Sources
=====================

Additional to the software mentioned in README you'll need the
following programs installed on your system:

* automake >= 1.10

* autoconf >= 2.63

* libtool


Getting started
---------------

If these preliminaries are met, you should run

```
dunecontrol all
```

which will find all installed dune modules as well as all dune modules
(not installed) which sources reside in a subdirectory of the current
directory. Note that if dune is not installed properly you will either
have to add the directory where the dunecontrol script resides (probably
./dune-common/bin) to your path or specify the relative path of the script.

On your project and all uninstalled DUNE source modules found the script
will then calls the GNU autoconf/automake to create a ./configure-script
and the Makefiles. Afterwards that configure script will be called and the
modules will be build using make all

Most probably you'll have to provide additional information to dunecontrol
(e. g. compilers, configure options) and/or make options.

The most convenient way is to use options files in this case. The files
defining three variables:

AUTOGEN_FLAGS    flags passed to autogen
CONFIGURE_FLAGS  flags passed to configure
MAKE_FLAGS       flags passed to make

An example options file might look like this:

#use this options to autogen, configure and make if no other options are given
AUTOGEN_FLAGS="--ac=2.50 --ac=1.8" #Forces automake 2,50 and autoconf 1.8
CONFIGURE_FLAGS="CXX=g++-3.4 --prefix=/install/path" #force g++-3.4 as compiler
MAKE_FLAGS=install #Per default run make install instead of simply make

If you save this information into example.opts you can path the opts file to
dunecontrol via the --opts option, e. g.

```
dunecontrol --opts=example.opts all
```

To get a full list of available configure flags just run

```
dunecontrol configure --help
```

after running at least `dunecontrol autogen`.


CMake support
-------------

If you have an unexplicable aversion to autotools, you can also build the library
using [CMake](cmake.org). If you want to do so, add `--use-cmake` to any dunecontrol
calls. Further information can be found at http://users.dune-project.org/projects/cmakedune.


Contributing patches
--------------------

If you have fixed a bug and would like to contribute it back to the project, we're happy
to take your patches! Please be aware that you will have to accept our license arrangements
(see file COPYING). Moreover, we have a small set of code style guidelines:

* No trailing whitespace

* Indentation with 2 spaces - no tabs please!

If you call `dunecontrol all` on your module, the script will automatically install a Git
commit hook that checks any commits you create locally for those guidelines. Alternatively,
you can explicitly request installation of the hook by calling `dunecontrol vcsetup`.

The easiest way to contribute changes is by opening a [pull request](https://github.com/smuething/dune-multidomaingrid/pulls) on GitHub, but if you don't like GitHub, just create the commits
locally and generate patch files with `git format-patch` (see https://www.kernel.org/pub/software/scm/git/docs/git-format-patch.html for details). Please make sure to set a valid Git identity; that way, your name will be recorded in the Git history for stuff like copyright notices.


More info
---------

See

```
dunecontrol --help
```

for further options.


The full build-system is described in dune-common/doc/buildsystem (Git version) or under share/doc/dune-common/buildsystem if you installed DUNE!
