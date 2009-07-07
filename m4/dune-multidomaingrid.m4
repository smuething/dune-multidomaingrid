# Additional checks needed to build the module
AC_DEFUN([DUNE_MULTIDOMAINGRID_CHECKS],[
  AX_BOOST_BASE([1.38.0])
])
# Additional checks needed to find the module
AC_DEFUN([DUNE_MULTIDOMAINGRID_CHECK_MODULE],[
  DUNE_CHECK_MODULES([dune-multidomaingrid], [grid/multidomaingrid.hh])
])
