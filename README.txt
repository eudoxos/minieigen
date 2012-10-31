This is miniEigen, small wrapper for the http://eigen.tuxfamily.org library.
Home of this project is at http://www.launchpad.net/minieigen/. The code is
licensed under the LGPLv3. The author is Václav Šmilauer <eu@doxos.eu>.

Run `python setup.py install` to build & install this module. You will need the
boost_python library to be installed, and Eigen (in /usr/include/eigen3, or
somewhere where the compiler finds it).

Windows is experimentally supported, but paths must be tweaked in setup.py
to match your system. Suggestions for a better way are welcome.

To compile by hand under Linux, run something like

	g++ -ansi src/miniEigen.cpp src/double-conversion/*.cc -o minieigen.so -shared -fPIC `pkg-config python --cflags` -lboost_python -I/usr/include/eigen3

Enjoy.
