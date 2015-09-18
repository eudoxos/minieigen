#pragma once
#include"common.hpp"
/*** automatic conversions from non-eigen types (sequences) ***/

/* template to define custom converter from sequence/list or approriate length and type, to eigen's Vector
   - length is stored in VT::RowsAtCompileTime
	- type is VT::Scalar
*/

template<typename T>
bool pySeqItemCheck(PyObject* o, int i){ return py::extract<T>(py::object(py::handle<>(PySequence_GetItem(o,i)))).check(); }

template<typename T>
T pySeqItemExtract(PyObject* o, int i){ return py::extract<T>(py::object(py::handle<>(PySequence_GetItem(o,i))))(); }

template<class VT>
struct custom_VectorAnyAny_from_sequence{
	custom_VectorAnyAny_from_sequence(){ py::converter::registry::push_back(&convertible,&construct,py::type_id<VT>()); }
	static void* convertible(PyObject* obj_ptr){ if(!PySequence_Check(obj_ptr) || (VT::RowsAtCompileTime!=Eigen::Dynamic && (PySequence_Size(obj_ptr)!=VT::RowsAtCompileTime))) return 0;
		// check that sequence items are convertible to scalars (should be done in other converters as well?!); otherwise Matrix3 is convertible to Vector3, but then we fail in *construct* very unclearly (TypeError: No registered converter was able to produce a C++ rvalue of type double from this Python object of type Vector3)
		size_t len=PySequence_Size(obj_ptr);
		for(size_t i=0; i<len; i++) if(!pySeqItemCheck<typename VT::Scalar>(obj_ptr,i)) return 0;
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, py::converter::rvalue_from_python_stage1_data* data){
		void* storage=((py::converter::rvalue_from_python_storage<VT>*)(data))->storage.bytes;
		new (storage) VT; size_t len;
		if(VT::RowsAtCompileTime!=Eigen::Dynamic){ len=VT::RowsAtCompileTime; }
		else{ len=PySequence_Size(obj_ptr); ((VT*)storage)->resize(len); }
		for(size_t i=0; i<len; i++) (*((VT*)storage))[i]=pySeqItemExtract<typename VT::Scalar>(obj_ptr,i);
		data->convertible=storage;
	}
};

template<class MT>
struct custom_MatrixAnyAny_from_sequence{
	custom_MatrixAnyAny_from_sequence(){ py::converter::registry::push_back(&convertible,&construct,py::type_id<MT>()); }
	static void* convertible(PyObject* obj_ptr){
		if(!PySequence_Check(obj_ptr)) return 0;
		bool isFlat=!PySequence_Check(py::handle<>(PySequence_GetItem(obj_ptr,0)).get());
		// mixed static/dynamic not handled (also not needed)
		BOOST_STATIC_ASSERT(
			(MT::RowsAtCompileTime!=Eigen::Dynamic && MT::ColsAtCompileTime!=Eigen::Dynamic)
			||
			(MT::RowsAtCompileTime==Eigen::Dynamic && MT::ColsAtCompileTime==Eigen::Dynamic)
		);
		int sz=PySequence_Size(obj_ptr);
		if(MT::RowsAtCompileTime!=Eigen::Dynamic){
			if(isFlat){
				// flat sequence (first item not sub-sequence), must contain exactly all items
				if(sz!=MT::RowsAtCompileTime*MT::ColsAtCompileTime) return 0;
			} else {
				// contains nested sequences, one per row
				if(sz!=MT::RowsAtCompileTime) return 0;
			}
		};
		return obj_ptr;
		// other checks done in the construct function
		// FIXME: it may be too late to do it there, as overloads are chosen based on *convertible*
		// (at least a clear message should be given when py::extract fails there)
	}
	static void construct(PyObject* obj_ptr, py::converter::rvalue_from_python_stage1_data* data){
		void* storage=((py::converter::rvalue_from_python_storage<MT>*)(data))->storage.bytes;
		new (storage) MT;
		MT &mx=*(MT*)storage;
		int sz=PySequence_Size(obj_ptr);
		bool isFlat=!PySequence_Check(py::handle<>(PySequence_GetItem(obj_ptr,0)).get());
		if(MT::RowsAtCompileTime!=Eigen::Dynamic){
			// do nothing
		} else {
			// find the right size
			if(isFlat) mx.resize(sz,1); // row vector, if flat
			else{ // find maximum size of items
				int rows=sz; int cols=0;
				for(int i=0; i<rows; i++){
					if(!PySequence_Check(py::handle<>(PySequence_GetItem(obj_ptr,i)).get())) throw std::runtime_error("Some elements of the array given are not sequences");
					int cols2=PySequence_Size(py::handle<>(PySequence_GetItem(obj_ptr,i)).get());
					if(cols==0) cols=cols2;
					if(cols!=cols2) throw std::runtime_error("Not all sub-sequences have the same length when assigning dynamic-sized matrix.");
				}
				mx.resize(rows,cols);
			}
		}
		if(isFlat){
			if(sz!=mx.rows()*mx.cols()) throw std::runtime_error("Assigning matrix "+lexical_cast<string>(mx.rows())+"x"+lexical_cast<string>(mx.cols())+" from flat vector of size "+lexical_cast<string>(sz));
			for(int i=0; i<sz; i++){
				mx(i/mx.rows(),i%mx.cols())=pySeqItemExtract<typename MT::Scalar>(obj_ptr,i);
			}
		} else {
			for(Index row=0; row<mx.rows(); row++){
				if(row>=PySequence_Size(obj_ptr)) throw std::runtime_error("Sequence rows of size "+lexical_cast<string>(sz)+" too short for assigning matrix with "+lexical_cast<string>(mx.rows())+" rows.");
				py::handle<> rowSeq(PySequence_GetItem(obj_ptr,row));
				if(!PySequence_Check(rowSeq.get())) throw std::runtime_error("Element of row sequence not a sequence.");
				if(mx.cols()!=PySequence_Size(rowSeq.get())) throw std::runtime_error("Row "+lexical_cast<string>(row)+": should specify exactly "+lexical_cast<string>(mx.cols())+" numbers, has "+lexical_cast<string>(PySequence_Size(rowSeq.get())));
				for(Index col=0; col<mx.cols(); col++){
					mx(row,col)=pySeqItemExtract<typename MT::Scalar>(rowSeq.get(),col);
				}
			}
		}
		data->convertible=storage;
	}
};

// create AlignedBoxNr from tuple of 2 Vector3r's
template<int dim>
struct custom_alignedBoxNr_from_seq{
	typedef Eigen::AlignedBox<Real,dim> AlignedBoxNr;
	typedef Eigen::Matrix<Real,dim,1> VectorNr;
	custom_alignedBoxNr_from_seq(){
		py::converter::registry::push_back(&convertible,&construct,py::type_id<AlignedBoxNr>());
	}
	static void* convertible(PyObject* obj_ptr){
		 if(!PySequence_Check(obj_ptr)) return 0;
		 if(PySequence_Size(obj_ptr)!=2) return 0;
		 if(!pySeqItemCheck<VectorNr>(obj_ptr,0) || !pySeqItemCheck<VectorNr>(obj_ptr,1)) return 0;
		 return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, py::converter::rvalue_from_python_stage1_data* data){
		void* storage=((py::converter::rvalue_from_python_storage<AlignedBoxNr>*)(data))->storage.bytes;
		new (storage) AlignedBoxNr(pySeqItemExtract<VectorNr>(obj_ptr,0),pySeqItemExtract<VectorNr>(obj_ptr,1));
		data->convertible=storage;
	}
};

struct custom_Quaternionr_from_axisAngle_or_angleAxis{
	custom_Quaternionr_from_axisAngle_or_angleAxis(){
		py::converter::registry::push_back(&convertible,&construct,py::type_id<Quaternionr>());
	}
	static void* convertible(PyObject* obj_ptr){
		if(!PySequence_Check(obj_ptr)) return 0;
		if(PySequence_Size(obj_ptr)!=2) return 0;
		py::object a(py::handle<>(PySequence_GetItem(obj_ptr,0))), b(py::handle<>(PySequence_GetItem(obj_ptr,1)));
		// axis-angle or angle-axis
		if((py::extract<Vector3r>(a).check() && py::extract<Real>(b).check()) || (py::extract<Real>(a).check() && py::extract<Vector3r>(b).check())) return obj_ptr;
		return 0;
	}
	static void construct(PyObject* obj_ptr, py::converter::rvalue_from_python_stage1_data* data){
		void* storage=((py::converter::rvalue_from_python_storage<Quaternionr>*)(data))->storage.bytes;
		py::object a(py::handle<>(PySequence_GetItem(obj_ptr,0))), b(py::handle<>(PySequence_GetItem(obj_ptr,1)));
		if(py::extract<Vector3r>(py::object(a)).check()) new (storage) Quaternionr(AngleAxisr(py::extract<Real>(b)(),py::extract<Vector3r>(a)().normalized()));
		else new (storage) Quaternionr(AngleAxisr(py::extract<Real>(a)(),py::extract<Vector3r>(b)().normalized()));
		data->convertible=storage;
	}
};



