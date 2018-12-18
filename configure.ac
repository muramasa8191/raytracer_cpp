#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([raytracer], [1.0], [muramasa8191@gmail.com])
AC_CONFIG_MACRO_DIRS([m4])

AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_SRCDIR([lib/stb/stb_image_write.h])
AC_CONFIG_HEADERS([config.h])

m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

LT_INIT

# Checks for programs.
AC_PROG_CXX([clang++ g++])
AC_PROG_CC([clang gcc])
AC_PROG_INSTALL

AC_LANG([C++], [C])
# Checks for libraries.

# Check for C++11 supports
#   ext: use extended mode (e.g. -std=gnu++11)
#   noext: use non-extended mode (e.g. -std=c++11)
AX_CXX_COMPILE_STDCXX_11([noext])

# Checks for header files.
AC_CHECK_HEADERS([omp.h])
AC_CHECK_HEADERS([OpenCL/cl.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_CHECK_FUNCS([memset pow sqrt strtol])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
