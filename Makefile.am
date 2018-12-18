ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS= -Wall
AM_CXXFLAGS= -Wall

AM_CPPFLAGS = \
  -I$(top_srcdir)/lib \
  -I$(top_srcdir)/src

lib_LTLIBRARIES = libraytracer.la
libraytracer_la_SOURCES = \
  src/camera.hpp \
  src/image.hpp \
  src/image.cpp \
  src/vector.hpp \
  src/ray.hpp \
  src/utils.hpp \
  src/material/material.hpp \
  src/shape/sphere.hpp \
  src/shape/sphere.cpp \
  src/shape/hitable.hpp \
  src/shape/shape_list.hpp \
  src/shape/shape_list.cpp \
  src/material/lambertian.hpp \
  src/material/dielectric.hpp \
  src/material/metal.hpp

bin_PROGRAMS = raytracer
raytracer_SOURCES = \
  stb/stb_image.h \
  stb/stb_image_write.h \
  src/raytracer.cpp

raytracer_CPPFLAGS = \
  -I$(top_srcdir)/src \
  -I$(top_srcdir)/lib \
  -fopenmp

raytracer_LDFLAGS = /usr/local/opt/llvm/lib/libgomp.dylib

raytracer_LDADD = libraytracer.la
