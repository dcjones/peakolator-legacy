
AC_INIT( [peakolator], [0.1], [dcjones@cs.washington.edu] )
AM_INIT_AUTOMAKE( [foreign -Wall -Werror] )
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

AC_CONFIG_MACRO_DIR( [m4] )

AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CXX

AC_DISABLE_STATIC
AC_PROG_LIBTOOL


AC_DEFINE( _FILE_OFFSET_BITS, 64 )
AC_DEFINE( GSL_RANGE_CHECK_OFF )


# zlib
AX_CHECK_ZLIB
               
# python
AX_PYTHON_DEVEL( [>= '2.6'] ) 
AX_PYTHON_MODULE( numpy, 1 )
AX_PYTHON_MODULE( argparse, 1 )


# cython
AM_CHECK_CYTHON

# boost
AX_BOOST_BASE

# gsl
AX_PATH_GSL

# gmp
AC_CHECK_LIB(gmp, __gmpz_init, , [AC_MSG_ERROR(
[GNU MP not found, see http://swox.com/gmp])])

# mpfr
AC_CHECK_LIB(mpfr, mpfr_add, , [AC_MSG_ERROR(
[Need MPFR either from GNU MP 4 or separate MPFR package.
See http://www.mpfr.org or http://swox.com/gmp])])

# nlopt
AC_CHECK_LIB(nlopt, nlopt_create, , [AC_MSG_ERROR(
[NLopt not found, see http://ab-initio.mit.edu/wiki/index.php/NLopt])])


# TODO: make sure large file support is turned on

AC_CONFIG_FILES( [Makefile
                  src/Makefile
                  src/gmpfrxx/Makefile
                  src/samtools/Makefile
                  src/yaml-cpp/Makefile
                  ])
AC_OUTPUT



