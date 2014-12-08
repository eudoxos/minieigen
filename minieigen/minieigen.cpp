// 2009-2012 © Václav Šmilauer <eu@doxos.eu>
// licensed under the Lesser General Public License version 3 (LGPLv3)

/* TODO:
	* Figure out if aligned types can be wrapped	(we failed at this previously; unaligned types force the c++ part not to align, making code perhaps less efficient numerically)
	* Add converters from 1-column MatrixX to VectorX so that matrix eqs work as expected
	* Figure out if integer types are ints or longs
*/

/*
The code is split to live in several files to reduce the amount of RAM necessary for compilation -- see http://www.boost.org/doc/libs/1_52_0/libs/python/doc/v2/faq.html#slow_compilation for the suggestion of this technique.
*/

#include"common.hpp"
#include"expose.hpp"

BOOST_PYTHON_MODULE(minieigen){
	py::scope().attr("__doc__")="miniEigen is wrapper for a small part of the `Eigen <http://eigen.tuxfamily.org>`_ library. Refer to its documentation for details. All classes in this module support pickling.";

	py::docstring_options docopt;
	docopt.enable_all();
	docopt.disable_cpp_signatures();



	expose_converters(); // in expose-converters.cpp

	expose_vectors();
	expose_matrices(); // must come after vectors
	expose_complex();
	expose_quaternion();
	expose_boxes();

	py::def("float2str",&doubleToShortest,(py::arg("f"),py::arg("pad")=0),"Return the shortest string representation of *f* which will is equal to *f* when converted back to float. This function is only useful in Python prior to 3.0; starting from that version, standard string conversion does just that.");

	#ifdef EIGEN_DONT_ALIGN
		py::scope().attr("vectorize")=false;
	#else
		py::scope().attr("vectorize")=true;
	#endif
};








