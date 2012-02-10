default:
	g++ -ansi miniEigen.cpp -o miniEigen.so -shared -fPIC `pkg-config python --cflags` -lboost_python -I/usr/include/eigen3
