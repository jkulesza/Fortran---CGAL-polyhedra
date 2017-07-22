#!/bin/sh
g++ cgal_polyhedra.cc -std=c++11 -o cc_cgal_polyhedra.o -c -I/sw/include

gfortran-mp-5 cc_cgal_polyhedra.o cgal_polyhedra.f90 test.f90 -lc++ -L/sw/lib -lCGAL -o test

