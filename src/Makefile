# this makefile is only useful for testing
default:
	g++ -DEIGEN_DONT_ALIGN -ansi miniEigen.cpp -o miniEigen.so -shared -fPIC `pkg-config python --cflags` -lboost_python -I/usr/include/eigen3
	PYTHONPATH=. python -c 'from miniEigen import *; print Vector6.Zero()'

