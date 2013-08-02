# this makefile is only useful for testing
DCSRC=double-conversion/bignum.cc double-conversion/bignum-dtoa.cc double-conversion/cached-powers.cc double-conversion/diy-fp.cc double-conversion/double-conversion.cc double-conversion/fast-dtoa.cc double-conversion/fixed-dtoa.cc double-conversion/strtod.cc

default:
	clang++ -DEIGEN_DONT_ALIGN -ansi minieigen.cpp $(DCSRC) -o minieigen.so -shared -fPIC `pkg-config python --cflags` -lboost_python -I/usr/include/eigen3
	PYTHONPATH=. python -c 'from minieigen import *; print Vector6.Zero, Matrix6(Vector6(1.232222,3.111,2.333,1.1,1.3,1.44)); print AlignedBox3((0,1,2),(3.12,4.34,5.22));'
debug:
	g++ -DEIGEN_DONT_ALIGN -g -O0 -ansi minieigen.cpp $(DCSRC) -o minieigen.so -shared -fPIC `pkg-config python --cflags` -lboost_python -I/usr/include/eigen3
	PYTHONPATH=. gdb --args python -c 'from minieigen import *; print Vector6.Zero, Matrix6(Vector6(1.232222,3.111,2.333,1.1,1.3,1.44)); print AlignedBox3((0,1,2),(3.12,4.34,5.22));'

