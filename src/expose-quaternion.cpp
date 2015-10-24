#include"visitors.hpp"
void expose_quaternion(){
	py::class_<Quaternionr>("Quaternion","Quaternion representing rotation.\n\nSupported operations (``q`` is a Quaternion, ``v`` is a Vector3): ``q*q`` (rotation composition), ``q*=q``, ``q*v`` (rotating ``v`` by ``q``), ``q==q``, ``q!=q``.\n\nStatic attributes: ``Identity``.\n\n.. note:: Quaternion is represented as axis-angle when printed (e.g. ``Identity`` is ``Quaternion((1,0,0),0)``, and can also be constructed from the axis-angle representation. This is however different from the data stored inside, which can be accessed by indices ``[0]`` (:math:`x`), ``[1]`` (:math:`y`), ``[2]`` (:math:`z`), ``[3]`` (:math:`w`). To obtain axis-angle programatically, use :obj:`Quaternion.toAxisAngle` which returns the tuple.",py::init<>())
		.def(QuaternionVisitor<Quaternionr>())
	;
}
