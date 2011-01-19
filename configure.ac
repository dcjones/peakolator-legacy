
AC_INIT( [peakolator], [0.1], [dcjones@cs.washington.edu] )
AM_INIT_AUTOMAKE( [foreign -Wall -Werror] )
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

AC_CONFIG_MACRO_DIR( [m4] )

AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CXX

# optimization flags
opt_CFLAGS="-g -O3 -ffast-math"
dbg_CFLAGS="-g -O0"

AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
	                      [enable debugging info (default is no)])],
              [], [enable_debug=no])

AS_IF([test "x$enable_debug" = xyes],
      [CFLAGS="$dbg_CFLAGS"],
      [CFLAGS="$opt_CFLAGS"])
CXXFLAGS=$CFLAGS


AC_DISABLE_STATIC
AC_PROG_LIBTOOL


AC_DEFINE(_GNU_SOURCE) # for asprintf, in particular
AC_DEFINE(_FILE_OFFSET_BITS, 64)
AC_DEFINE(GSL_RANGE_CHECK_OFF)


# zlib
AX_CHECK_ZLIB
               
# python
AX_PYTHON_DEVEL( [>= '2.6'] ) 
AX_PYTHON_MODULE( numpy, 1 )
AX_PYTHON_MODULE( argparse, 1 )


# cython
AM_CHECK_CYTHON

# gsl
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
AX_PATH_GSL

# nlopt
AC_CHECK_LIB(nlopt, nlopt_create, , [AC_MSG_ERROR(
[NLopt not found, see http://ab-initio.mit.edu/wiki/index.php/NLopt])])


AC_CONFIG_FILES( [Makefile
                  src/Makefile
                  src/samtools/Makefile
                  src/yaml-cpp/Makefile
                  ])
AC_OUTPUT



