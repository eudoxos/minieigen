# this makefile is only useful for testing
DCSRC=double-conversion/bignum.cc double-conversion/bignum-dtoa.cc double-conversion/cached-powers.cc double-conversion/diy-fp.cc double-conversion/double-conversion.cc double-conversion/fast-dtoa.cc double-conversion/fixed-dtoa.cc double-conversion/strtod.cc
MESRC=minieigen.cpp expose-boxes.cpp expose-complex.cpp expose-converters.cpp expose-matrices.cpp expose-quaternion.cpp expose-vectors.cpp 

default:
	clang++ -DEIGEN_DONT_ALIGN -ansi $(MESRC) $(DCSRC) -o minieigen.so -shared -fPIC `pkg-config python --cflags` -ftemplate-backtrace-limit=0 -lboost_python -I/usr/include/eigen3
	python basic-test.py
	# PYTHONPATH=. python -c 'from minieigen import *; print Vector6.Zero, Matrix6(Vector6(1.232222,3.111,2.333,1.1,1.3,1.44)); print AlignedBox3((0,1,2),(3.12,4.34,5.22));'
	# for testing complex numbers
	# a=Vector6c(); a[0]=2+1j; a[1]=2j; a[2]=3; a[3]=0; print a
debug:
	g++ -DEIGEN_DONT_ALIGN -g -O0 -ansi $(MESRC) $(DCSRC) -o minieigen.so -shared -fPIC `pkg-config python --cflags` -lboost_python -I/usr/include/eigen3
	PYTHONPATH=. gdb --args python -c 'from minieigen import *; print Vector6.Zero, Matrix6(Vector6(1.232222,3.111,2.333,1.1,1.3,1.44)); print AlignedBox3((0,1,2),(3.12,4.34,5.22));'

