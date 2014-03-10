#include"visitors.hpp"
void expose_matrices(){
	py::class_<Matrix3r>("Matrix3","3x3 float matrix.\n\nSupported operations (``m`` is a Matrix3, ``f`` if a float/int, ``v`` is a Vector3): ``-m``, ``m+m``, ``m+=m``, ``m-m``, ``m-=m``, ``m*f``, ``f*m``, ``m*=f``, ``m/f``, ``m/=f``, ``m*m``, ``m*=m``, ``m*v``, ``v*m``, ``m==m``, ``m!=m``.\n\nStatic attributes: ``Zero``, ``Ones``, ``Identity``.",py::init<>())
		.def(py::init<Quaternionr const &>((py::arg("q"))))
		.def(MatrixVisitor<Matrix3r>())
	;
	py::class_<Matrix6r>("Matrix6","6x6 float matrix. Constructed from 4 3x3 sub-matrices, from 6xVector6 (rows).\n\nSupported operations (``m`` is a Matrix6, ``f`` if a float/int, ``v`` is a Vector6): ``-m``, ``m+m``, ``m+=m``, ``m-m``, ``m-=m``, ``m*f``, ``f*m``, ``m*=f``, ``m/f``, ``m/=f``, ``m*m``, ``m*=m``, ``m*v``, ``v*m``, ``m==m``, ``m!=m``.\n\nStatic attributes: ``Zero``, ``Ones``, ``Identity``.",py::init<>())
		.def(MatrixVisitor<Matrix6r>())
	;
	py::class_<MatrixXr>("MatrixX","XxX (dynamic-sized) float matrix. Constructed from list of rows (as VectorX).\n\nSupported operations (``m`` is a MatrixX, ``f`` if a float/int, ``v`` is a VectorX): ``-m``, ``m+m``, ``m+=m``, ``m-m``, ``m-=m``, ``m*f``, ``f*m``, ``m*=f``, ``m/f``, ``m/=f``, ``m*m``, ``m*=m``, ``m*v``, ``v*m``, ``m==m``, ``m!=m``.",py::init<>())
		.def(MatrixVisitor<MatrixXr>())
	;
}
