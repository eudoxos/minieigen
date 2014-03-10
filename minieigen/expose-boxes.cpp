#include"visitors.hpp"
void expose_boxes(){
	py::class_<AlignedBox3r>("AlignedBox3","Axis-aligned box object, defined by its minimum and maximum corners",py::init<>())
		.def(AabbVisitor<AlignedBox3r>())
	;

	py::class_<AlignedBox2r>("AlignedBox2","Axis-aligned box object in 2d, defined by its minimum and maximum corners",py::init<>())
		.def(AabbVisitor<AlignedBox2r>())
	;
}
