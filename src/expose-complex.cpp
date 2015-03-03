#include"visitors.hpp"
void expose_complex(){
	#ifdef _COMPLEX_SUPPORT
		py::class_<Vector2cr>("Vector2c","/*TODO*/",py::init<>()).def(VectorVisitor<Vector2cr>());
		py::class_<Vector3cr>("Vector3c","/*TODO*/",py::init<>()).def(VectorVisitor<Vector3cr>());
		py::class_<Vector6cr>("Vector6c","/*TODO*/",py::init<>()).def(VectorVisitor<Vector6cr>());
		py::class_<VectorXcr>("VectorXc","/*TODO*/",py::init<>()).def(VectorVisitor<VectorXcr>());

		py::class_<Matrix3cr>("Matrix3c","/*TODO*/",py::init<>()).def(MatrixVisitor<Matrix3cr>());
		py::class_<Matrix6cr>("Matrix6c","/*TODO*/",py::init<>()).def(MatrixVisitor<Matrix6cr>());
		py::class_<MatrixXcr>("MatrixXc","/*TODO*/",py::init<>()).def(MatrixVisitor<MatrixXcr>());
		;
	#endif
}
