# Additional checks needed to build the module
AC_DEFUN([DUNE_MULTIDOMAINGRID_CHECKS],[
  DUNE_BOOST_BASE(, [ DUNE_BOOST_FUSION ] , [] )
])
# Additional checks needed to find the module
AC_DEFUN([DUNE_MULTIDOMAINGRID_CHECK_MODULE],[
  DUNE_CHECK_MODULES([dune-multidomaingrid], [grid/multidomaingrid/singlevalueset.hh])
])
