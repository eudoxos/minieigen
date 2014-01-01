// 2009-2012 © Václav Šmilauer <eu@doxos.eu>
// licensed under the Lesser General Public License version 3 (LGPLv3)

/* TODO:
	* Figure out if aligned types can be wrapped	(we failed at this previously; unaligned types force the c++ part not to align, making code perhaps less efficient numerically)
	* Add converters from 1-column MatrixX to VectorX so that matrix eqs work as expected
	* Figure out if integer types are ints or longs
	* make indices properly use Eigen's Index type (where is it exposed in headers?) instead of ints -- to avoid narrowing warnings
*/


// if #defined, use visitors, which saves a lot of copied code and is necessary for complex numbers support
// once this is tested throughly, it will be the default and lots of old code will be removed :)
#define _NEW_VISITORS

// if #defined, complex support is enabled (requires _NEW_VISITORS)
#define _COMPLEX_SUPPORT

/* change to float for single-precision */
typedef double Real;

// BEGIN workaround for
// * http://eigen.tuxfamily.org/bz/show_bug.cgi?id=528
// * https://sourceforge.net/tracker/index.php?func=detail&aid=3584127&group_id=202880&atid=983354
// (only needed with gcc <= 4.7)
#include<stdlib.h>
#include<sys/stat.h>
// END workaround

#ifdef _COMPLEX_SUPPORT
	#ifndef _NEW_VISITORS
		#error _COMPLEX_SUPPORT requires _NEW_VISITORS
	#endif
#endif


#include<Eigen/Core>
#include<Eigen/Geometry>
#include<Eigen/Eigenvalues>
#include<Eigen/SVD>

/* exposed types */
typedef Eigen::Matrix<int ,2,1> Vector2i;
typedef Eigen::Matrix<Real,2,1> Vector2r;
typedef Eigen::Matrix<int ,3,1> Vector3i;
typedef Eigen::Matrix<Real,3,1> Vector3r;
typedef Eigen::Matrix<int ,6,1> Vector6i;
typedef Eigen::Matrix<Real,6,1> Vector6r;
typedef Eigen::Matrix<Real,3,3> Matrix3r;
typedef Eigen::Matrix<Real,6,6> Matrix6r;

typedef Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> MatrixXr;
typedef Eigen::Matrix<Real,Eigen::Dynamic,1> VectorXr;

typedef Eigen::Quaternion<Real> Quaternionr;
typedef Eigen::AngleAxis<Real> AngleAxisr;
typedef Eigen::AlignedBox<Real,3> AlignedBox3r;
typedef Eigen::AlignedBox<Real,2> AlignedBox2r;

#ifdef _COMPLEX_SUPPORT
#include<complex>
	using std::complex;
	typedef Eigen::Matrix<complex<Real>,2,1> Vector2cr;
	typedef Eigen::Matrix<complex<Real>,3,1> Vector3cr;
	typedef Eigen::Matrix<complex<Real>,6,1> Vector6cr;
	typedef Eigen::Matrix<complex<Real>,Eigen::Dynamic,1> VectorXcr;
	typedef Eigen::Matrix<complex<Real>,3,3> Matrix3cr;
	typedef Eigen::Matrix<complex<Real>,6,6> Matrix6cr;
	typedef Eigen::Matrix<complex<Real>,Eigen::Dynamic,Eigen::Dynamic> MatrixXcr;
#endif


#include<string>
using std::string;
#include<stdexcept>
#include<sstream>
#include<iomanip>
#include<vector>

#include<boost/python.hpp>
namespace py=boost::python;
#include<boost/lexical_cast.hpp>
using boost::lexical_cast;
#include<boost/static_assert.hpp>

/**** double-conversion helpers *****/
#include"double-conversion/double-conversion.h"

double_conversion::DoubleToStringConverter doubleToString(
	double_conversion::DoubleToStringConverter::NO_FLAGS,
	"inf", /* infinity symbol */
	"nan", /* NaN symbol */
	'e', /*exponent symbol*/
	-5, /* decimal_in_shortest_low: 0.0001, but 0.00001->1e-5 */
	7, /* decimal_in_shortest_high */
	/* the following are irrelevant for the shortest representation */
	6, /* max_leading_padding_zeroes_in_precision_mode */
	6 /* max_trailing_padding_zeroes_in_precision_mode */
);

/* optionally pad from the left */
string doubleToShortest(double d, int pad=0){
	/* 32 is perhaps wasteful */
	/* it would be better to write to the string's buffer itself, not sure how to do that */
	char buf[32];
	double_conversion::StringBuilder sb(buf,32);
	doubleToString.ToShortest(d,&sb);
	string ret(sb.Finalize());
	if(pad==0 || (int)ret.size()>=pad) return ret;
	return string(pad-ret.size(),' ')+ret; // left-padded if shorter
} 


/* generic function to print numbers, via lexical_cast plus padding -- used for ints */
template<typename T>
string num_to_string(const T& num, int pad=0){
	string ret(lexical_cast<string>(num));
	if(pad==0 || (int)ret.size()>=pad) return ret;
	return string(pad-ret.size(),' ')+ret; // left-pad with spaces
}

// for doubles, use the shortest representation
string num_to_string(const double& num, int pad=0){ return doubleToShortest(num,pad); }

#ifdef _COMPLEX_SUPPORT
	// for complex numbers (with any scalar type, though only doubles are really used)
	template<typename T>
	string num_to_string(const complex<T>& num, int pad=0){
		string ret;
		// both components non-zero
		if(num.real()!=0 && num.imag()!=0){
			// don't add "+" in the middle if imag is negative and will start with "-"
			string ret=num_to_string(num.real(),/*pad*/0)+(num.imag()>0?"+":"")+num_to_string(num.imag(),/*pad*/0)+"j";
			if(pad==0 || (int)ret.size()>=pad) return ret;
			return string(pad-ret.size(),' ')+ret; // left-pad with spaces
		}
		// only imaginary is non-zero: skip the real part, and decrease padding to accomoadate the trailing "j"
		if(num.imag()!=0){
			return num_to_string(num.imag(),/*pad*/pad>0?pad-1:0)+"j";
		}
		// non-complex (zero or not)
		return num_to_string(num.real(),pad);
	}
#endif


/*** getters and setters with bound guards ***/
void IDX_CHECK(int i,int MAX){ if(i<0 || i>=MAX) { PyErr_SetString(PyExc_IndexError,("Index "+lexical_cast<string>(i)+" out of range 0.." + lexical_cast<string>(MAX-1)).c_str()); py::throw_error_already_set(); } }
void IDX2_CHECKED_TUPLE_INTS(py::tuple tuple,const int max2[2], int arr2[2]) {int l=py::len(tuple); if(l!=2) { PyErr_SetString(PyExc_IndexError,"Index must be integer or a 2-tuple"); py::throw_error_already_set(); } for(int _i=0; _i<2; _i++) { py::extract<int> val(tuple[_i]); if(!val.check()){ PyErr_SetString(PyExc_ValueError,("Unable to convert "+lexical_cast<string>(_i)+"-th index to int.").c_str()); py::throw_error_already_set(); } int v=val(); IDX_CHECK(v,max2[_i]); arr2[_i]=v; }  }


/*** automatic conversions from non-eigen types (sequences) ***/

/* template to define custom converter from sequence/list or approriate length and type, to eigen's Vector
   - length is stored in VT::RowsAtCompileTime
	- type is VT::Scalar
*/
template<class VT>
struct custom_VectorAnyAny_from_sequence{
	custom_VectorAnyAny_from_sequence(){ py::converter::registry::push_back(&convertible,&construct,py::type_id<VT>()); }
	static void* convertible(PyObject* obj_ptr){ if(!PySequence_Check(obj_ptr) || (VT::RowsAtCompileTime!=Eigen::Dynamic && (PySequence_Size(obj_ptr)!=VT::RowsAtCompileTime))) return 0; return obj_ptr; }
	static void construct(PyObject* obj_ptr, py::converter::rvalue_from_python_stage1_data* data){
		void* storage=((py::converter::rvalue_from_python_storage<VT>*)(data))->storage.bytes;
		new (storage) VT; size_t len;
		if(VT::RowsAtCompileTime!=Eigen::Dynamic){ len=VT::RowsAtCompileTime; }
		else{ len=PySequence_Size(obj_ptr); ((VT*)storage)->resize(len); }
		for(size_t i=0; i<len; i++) (*((VT*)storage))[i]=py::extract<typename VT::Scalar>(PySequence_GetItem(obj_ptr,i));
		data->convertible=storage;
	}
};

template<class MT>
struct custom_MatrixAnyAny_from_sequence{
	custom_MatrixAnyAny_from_sequence(){ py::converter::registry::push_back(&convertible,&construct,py::type_id<MT>()); }
	static void* convertible(PyObject* obj_ptr){
		if(!PySequence_Check(obj_ptr)) return 0;
		PySequence_GetItem(obj_ptr,0);
		bool isFlat=!PySequence_Check(PySequence_GetItem(obj_ptr,0));
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
	}
	static void construct(PyObject* obj_ptr, py::converter::rvalue_from_python_stage1_data* data){
		void* storage=((py::converter::rvalue_from_python_storage<MT>*)(data))->storage.bytes;
		new (storage) MT;
		MT &mx=*(MT*)storage;
		int sz=PySequence_Size(obj_ptr);
		bool isFlat=!PySequence_Check(PySequence_GetItem(obj_ptr,0));
		if(MT::RowsAtCompileTime!=Eigen::Dynamic){
			// do nothing
		} else {
			// find the right size
			if(isFlat) mx.resize(sz,1); // row vector, if flat
			else{ // find maximum size of items
				int rows=sz; int cols=0;
				for(int i=0; i<rows; i++){
					if(!PySequence_Check(PySequence_GetItem(obj_ptr,i))) throw std::runtime_error("Some elements of the array given are not sequences");
					int cols2=PySequence_Size(PySequence_GetItem(obj_ptr,i));
					if(cols==0) cols=cols2;
					if(cols!=cols2) throw std::runtime_error("Not all sub-sequences have the same length when assigning dynamic-sized matrix.");
				}
				mx.resize(rows,cols);
			}
		}
		if(isFlat){
			if(sz!=mx.rows()*mx.cols()) throw std::runtime_error("Assigning matrix "+lexical_cast<string>(mx.rows())+"x"+lexical_cast<string>(mx.cols())+" from flat vector of size "+lexical_cast<string>(sz));
			for(int i=0; i<sz; i++){
				mx(i/mx.rows(),i%mx.cols())=py::extract<typename MT::Scalar>(PySequence_GetItem(obj_ptr,i));
			}
		} else {
			for(int row=0; row<mx.rows(); row++){
				if(row>=PySequence_Size(obj_ptr)) throw std::runtime_error("Sequence rows of size "+lexical_cast<string>(sz)+" too short for assigning matrix with "+lexical_cast<string>(mx.rows())+" rows.");
				PyObject* rowSeq=PySequence_GetItem(obj_ptr,row);
				if(!PySequence_Check(rowSeq)) throw std::runtime_error("Element of row sequence not a sequence.");
				if(mx.cols()!=PySequence_Size(rowSeq)) throw std::runtime_error("Row "+lexical_cast<string>(row)+": should specify exactly "+lexical_cast<string>(mx.cols())+" numbers, has "+lexical_cast<string>(PySequence_Size(rowSeq)));
				for(int col=0; col<mx.cols(); col++){
					mx(row,col)=py::extract<typename MT::Scalar>(PySequence_GetItem(rowSeq,col));
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
		 if(!py::extract<VectorNr>(PySequence_GetItem(obj_ptr,0)).check() || !py::extract<VectorNr>(PySequence_GetItem(obj_ptr,1)).check()) return 0;
		 return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, py::converter::rvalue_from_python_stage1_data* data){
		void* storage=((py::converter::rvalue_from_python_storage<AlignedBoxNr>*)(data))->storage.bytes;
		new (storage) AlignedBoxNr(py::extract<VectorNr>(PySequence_GetItem(obj_ptr,0))(),py::extract<VectorNr>(PySequence_GetItem(obj_ptr,1))());
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
		PyObject *a(PySequence_GetItem(obj_ptr,0)), *b(PySequence_GetItem(obj_ptr,1));
		// axis-angle or angle-axis
		if((py::extract<Vector3r>(a).check() && py::extract<Real>(b).check()) || (py::extract<Real>(a).check() && py::extract<Vector3r>(b).check())) return obj_ptr;
		return 0;
	}
	static void construct(PyObject* obj_ptr, py::converter::rvalue_from_python_stage1_data* data){
		void* storage=((py::converter::rvalue_from_python_storage<Quaternionr>*)(data))->storage.bytes;
		PyObject *a(PySequence_GetItem(obj_ptr,0)), *b(PySequence_GetItem(obj_ptr,1));
		if(py::extract<Vector3r>(a).check()) new (storage) Quaternionr(AngleAxisr(py::extract<Real>(b)(),py::extract<Vector3r>(a)().normalized()));
		else new (storage) Quaternionr(AngleAxisr(py::extract<Real>(a)(),py::extract<Vector3r>(b)().normalized()));
		data->convertible=storage;
	}
};


#ifdef _NEW_VISITORS

static string object_class_name(const py::object& obj){ return py::extract<string>(obj.attr("__class__").attr("__name__"))(); }

// methods common for vectors and matrices
template<typename MatrixBaseT>
class MatrixBaseVisitor: public py::def_visitor<MatrixBaseVisitor<MatrixBaseT> >{
	typedef typename MatrixBaseT::Scalar Scalar;
	public:
	template<class PyClass>
	void visit(PyClass& cl) const {
		cl
		.def(py::init<MatrixBaseT>(py::arg("other")))
		.def("__neg__",&MatrixBaseVisitor::__neg__)
		.def("__add__",&MatrixBaseVisitor::__add__).def("__iadd__",&MatrixBaseVisitor::__iadd__)
		.def("__sub__",&MatrixBaseVisitor::__sub__).def("__isub__",&MatrixBaseVisitor::__isub__)
		.def("__eq__",&MatrixBaseVisitor::__eq__).def("__ne__",&MatrixBaseVisitor::__ne__)
		.def("__mul__",&MatrixBaseVisitor::__mul__scalar<long>)
		.def("__imul__",&MatrixBaseVisitor::__imul__scalar<long>)
		.def("__rmul__",&MatrixBaseVisitor::__rmul__scalar<long>)
		.def("rows",&MatrixBaseT::rows,"Number of rows.")
		.def("cols",&MatrixBaseT::cols,"Number of columns.")
		;
		visit_if_float<Scalar,PyClass>(cl);
		visit_fixed_or_dynamic<MatrixBaseT,PyClass>(cl);

		// reductions
		cl
		.def("sum",&MatrixBaseT::sum,"Sum of all elements.")
		.def("maxAbsCoeff",&MatrixBaseVisitor::maxAbsCoeff,"Maximum absolute value over all elements.")
		;
	};
	private:
	template<class PyClass> static string name(PyClass& cl){ return py::extract<string>(cl.attr("__name__"))(); }
	// for dynamic matrices/vectors
	template<typename MatrixBaseT2, class PyClass> static void visit_fixed_or_dynamic(PyClass& cl, typename boost::enable_if_c<MatrixBaseT2::RowsAtCompileTime==Eigen::Dynamic>::type* dummy = 0){
		//std::cerr<<"MatrixBaseVisitor: dynamic MatrixBaseT for "<<name(cl)<<std::endl;
		;
	}
	// for static matrices/vectors
	template<typename MatrixBaseT2, class PyClass> static void visit_fixed_or_dynamic(PyClass& cl, typename boost::disable_if_c<MatrixBaseT2::RowsAtCompileTime==Eigen::Dynamic>::type* dummy = 0){
		//std::cerr<<"MatrixBaseVisitor: fixed MatrixBaseT for "<<name(cl)<<std::endl;
		cl
		.add_static_property("Ones",&MatrixBaseVisitor::Ones)
		.add_static_property("Zero",&MatrixBaseVisitor::Zero)
		.def("Random",&MatrixBaseVisitor::Random,"Return an object where all elements are randomly set to values between 0 and 1.").staticmethod("Random")
		.add_static_property("Identity",&MatrixBaseVisitor::Identity)
		;
	}
	template<typename Scalar, class PyClass> static	void visit_if_float(PyClass& cl, typename boost::enable_if<boost::is_integral<Scalar> >::type* dummy = 0){ /* do nothing */ }
	template<typename Scalar, class PyClass> static void visit_if_float(PyClass& cl, typename boost::disable_if<boost::is_integral<Scalar> >::type* dummy = 0){
		// operations with other scalars (Scalar is the floating type, long is the python integer type)
		cl.def("__mul__",&MatrixBaseVisitor::__mul__scalar<Scalar>)
		.def("__rmul__",&MatrixBaseVisitor::__rmul__scalar<Scalar>)
		.def("__imul__",&MatrixBaseVisitor::__imul__scalar<Scalar>)
		.def("__div__",&MatrixBaseVisitor::__div__scalar<long>)
		.def("__idiv__",&MatrixBaseVisitor::__idiv__scalar<long>)
		.def("__div__",&MatrixBaseVisitor::__div__scalar<Scalar>)
		.def("__idiv__",&MatrixBaseVisitor::__idiv__scalar<Scalar>)
		//
		.def("norm",&MatrixBaseT::norm,"Euclidean norm.")
		.def("__abs__",&MatrixBaseT::norm)
		.def("squaredNorm",&MatrixBaseT::squaredNorm,"Square of the Euclidean norm.")
		.def("normalize",&MatrixBaseT::normalize,"Normalize this object in-place.")
		.def("normalized",&MatrixBaseT::normalized,"Return normalized copy of this object")
		.def("pruned",&MatrixBaseVisitor::pruned,py::arg("absTol")=1e-6,"Zero all elements which are greater than *absTol*. Negative zeros are not pruned.")
		;
	}
	// for fixed-size matrices/vectors only
	static Scalar maxAbsCoeff(const MatrixBaseT& m){ return m.array().abs().maxCoeff(); }
	static MatrixBaseT Ones(){ return MatrixBaseT::Ones(); }
	static MatrixBaseT Zero(){ return MatrixBaseT::Zero(); }
	static MatrixBaseT Random(){ return MatrixBaseT::Random(); }
	static MatrixBaseT Identity(){ return MatrixBaseT::Identity(); }

	static bool __eq__(const MatrixBaseT& a, const MatrixBaseT& b){
		if(a.rows()!=b.rows() || a.cols()!=b.cols()) return false;
		return a.cwiseEqual(b).all();
	}
	static bool __ne__(const MatrixBaseT& a, const MatrixBaseT& b){ return !__eq__(a,b); }
	static MatrixBaseT __neg__(const MatrixBaseT& a){ return -a; };
	static MatrixBaseT __add__(const MatrixBaseT& a, const MatrixBaseT& b){ return a+b; }
	static MatrixBaseT __sub__(const MatrixBaseT& a, const MatrixBaseT& b){ return a-b; }
	static MatrixBaseT __iadd__(MatrixBaseT& a, const MatrixBaseT& b){ a+=b; return a; };
	static MatrixBaseT __isub__(MatrixBaseT& a, const MatrixBaseT& b){ a-=b; return a; };

	template<typename Scalar2> static MatrixBaseT __mul__scalar(const MatrixBaseT& a, const Scalar2& scalar){ return a*scalar; }
	template<typename Scalar2> static void __imul__scalar(MatrixBaseT& a, const Scalar2& scalar){ a*=scalar; }
	template<typename Scalar2> static MatrixBaseT __rmul__scalar(const MatrixBaseT& a, const Scalar2& scalar){ return a*scalar; }
	template<typename Scalar2> static MatrixBaseT __div__scalar(const MatrixBaseT& a, const Scalar2& scalar){ return a/scalar; }
	template<typename Scalar2> static void __idiv__scalar(MatrixBaseT& a, const Scalar2& scalar){ a/=scalar; }

	// we want to keep -0 (rather than replacing it by 0), but that does not work for complex numbers
	// hence two versions
	template<typename Scalar> static bool prune_element(const Scalar& num, double absTol, typename boost::disable_if<boost::is_complex<Scalar> >::type* dummy=0){ return std::abs(num)<=absTol || num!=-0; }
	template<typename Scalar> static bool prune_element(const Scalar& num, double absTol, typename boost::enable_if<boost::is_complex<Scalar> >::type* dummy=0){ return std::abs(num)<=absTol; }
	
	static MatrixBaseT pruned(const MatrixBaseT& a, double absTol=1e-6){ // typename MatrixBaseT::Scalar absTol=1e-6){
		MatrixBaseT ret(MatrixBaseT::Zero(a.rows(),a.cols()));
		for(int c=0;c<a.cols();c++){ for(int r=0;r<a.rows();r++){ if(!prune_element(a(c,r),absTol)) ret(c,r)=a(c,r); } }
		return ret;
	};
};

template<typename VectorT>
class VectorVisitor: public py::def_visitor<VectorVisitor<VectorT> >{
	friend class def_visitor_access;
	typedef typename VectorT::Scalar Scalar;
	typedef Eigen::Matrix<Scalar,VectorT::RowsAtCompileTime,VectorT::RowsAtCompileTime> CompatMatrixT;
	typedef Eigen::Matrix<Scalar,2,1> CompatVec2;
	typedef Eigen::Matrix<Scalar,3,1> CompatVec3;
	typedef Eigen::Matrix<Scalar,6,1> CompatVec6;
	typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1> CompatVecX;
	enum{Dim=VectorT::RowsAtCompileTime};
	public:
	template<class PyClass>
	void visit(PyClass& cl) const {
		MatrixBaseVisitor<VectorT>().visit(cl);
		cl
		.def_pickle(VectorPickle())
		.def("__setitem__",&VectorVisitor::set_item)
		.def("__getitem__",&VectorVisitor::get_item)
		.def("__str__",&VectorVisitor::__str__).def("__repr__",&VectorVisitor::__str__)
		.def("dot",&VectorVisitor::dot,py::arg("other"),"Dot product with *other*.")
		.def("outer",&VectorVisitor::outer,py::arg("other"),"Outer product with *other*.")
		.def("asDiagonal",&VectorVisitor::asDiagonal,"Return diagonal matrix with this vector on the diagonal.")
		;

		visit_fixed_or_dynamic<VectorT,PyClass>(cl);

		visit_special_sizes<VectorT,PyClass>(cl);
	};
	private:
	// for dynamic vectors
	template<typename VectorT2, class PyClass> static void visit_fixed_or_dynamic(PyClass& cl, typename boost::enable_if_c<VectorT2::RowsAtCompileTime==Eigen::Dynamic>::type* dummy = 0){
		//std::cerr<<"VectorVisitor: dynamic vector for "<<name()<<std::endl;
		cl
		.def("__len__",&VectorVisitor::dyn__len__)
		.def("resize",&VectorVisitor::resize)
		.def("Unit",&VectorVisitor::dyn_Unit).staticmethod("Unit")
		.def("Ones",&VectorVisitor::dyn_Ones).staticmethod("Ones")
		.def("Zero",&VectorVisitor::dyn_Zero).staticmethod("Zero")
		.def("Random",&VectorVisitor::dyn_Random,py::arg("len"),"Return vector of given length with all elements set to values between 0 and 1 randomly.").staticmethod("Random")
		;
	}
	// for fixed-size vectors
	template<typename VectorT2, class PyClass> static void visit_fixed_or_dynamic(PyClass& cl, typename boost::disable_if_c<VectorT2::RowsAtCompileTime==Eigen::Dynamic>::type* dummy = 0){
		//std::cerr<<"VectorVisitor: fixed vector for "<<name()<<std::endl;
		cl.def("__len__",&VectorVisitor::__len__).staticmethod("__len__")
		.def("Unit",&VectorVisitor::Unit).staticmethod("Unit")
		;
	}

	// handle specific sizes of vectors separately

	// 2-vector
	template<typename VectorT2, class PyClass> static void visit_special_sizes(PyClass& cl, typename boost::enable_if_c<VectorT2::RowsAtCompileTime==2>::type* dummy=0){
		cl
		.def(py::init<typename VectorT2::Scalar,typename VectorT2::Scalar>((py::arg("x"),py::arg("y"))))
		.add_static_property("UnitX",&VectorVisitor::Vec2_UnitX)
		.add_static_property("UnitY",&VectorVisitor::Vec2_UnitY)
		;
	}
	static CompatVec2 Vec2_UnitX(){ return CompatVec2::UnitX(); }
	static CompatVec2 Vec2_UnitY(){ return CompatVec2::UnitY(); }

	// 3-vector
	template<typename VectorT2, class PyClass> static void visit_special_sizes(PyClass& cl, typename boost::enable_if_c<VectorT2::RowsAtCompileTime==3>::type* dummy=0){
		cl
		.def(py::init<typename VectorT2::Scalar,typename VectorT2::Scalar,typename VectorT2::Scalar>((py::arg("x"),py::arg("y"),py::arg("z"))))
		.def("cross",&VectorVisitor::cross) // cross-product only meaningful for 3-sized vectors
		.add_static_property("UnitX",&VectorVisitor::Vec3_UnitX)
		.add_static_property("UnitY",&VectorVisitor::Vec3_UnitY)
		.add_static_property("UnitZ",&VectorVisitor::Vec3_UnitZ)
		// swizzles
		.def("xy",&VectorVisitor::Vec3_xy).def("yx",&VectorVisitor::Vec3_yx).def("xz",&VectorVisitor::Vec3_xz).def("zx",&VectorVisitor::Vec3_zx).def("yz",&VectorVisitor::Vec3_yz).def("zy",&VectorVisitor::Vec3_zy)
		;
	}
	static CompatVec3 cross(const CompatVec3& self, const CompatVec3& other){ return self.cross(other); }
	static CompatVec3 Vec3_UnitX(){ return CompatVec3::UnitX(); }
	static CompatVec3 Vec3_UnitY(){ return CompatVec3::UnitY(); }
	static CompatVec3 Vec3_UnitZ(){ return CompatVec3::UnitZ(); }

	static CompatVec2 Vec3_xy(const CompatVec3& v){ return CompatVec2(v[0],v[1]); }
	static CompatVec2 Vec3_yx(const CompatVec3& v){ return CompatVec2(v[1],v[0]); }
	static CompatVec2 Vec3_xz(const CompatVec3& v){ return CompatVec2(v[0],v[2]); }
	static CompatVec2 Vec3_zx(const CompatVec3& v){ return CompatVec2(v[2],v[0]); }
	static CompatVec2 Vec3_yz(const CompatVec3& v){ return CompatVec2(v[1],v[2]); }
	static CompatVec2 Vec3_zy(const CompatVec3& v){ return CompatVec2(v[2],v[1]); }
	
	// 6-vector
	template<typename VectorT2, class PyClass> static void visit_special_sizes(PyClass& cl, typename boost::enable_if_c<VectorT2::RowsAtCompileTime==6>::type* dummy=0){
		cl
		.def("__init__",py::make_constructor(&VectorVisitor::Vec6_fromElements,py::default_call_policies(),(py::arg("v0"),py::arg("v1"),py::arg("v2"),py::arg("v3"),py::arg("v4"),py::arg("v5"))))
		.def("__init__",py::make_constructor(&VectorVisitor::Vec6_fromHeadTail,py::default_call_policies(),(py::arg("head"),py::arg("tail"))))
		.def("head",&VectorVisitor::Vec6_head).def("tail",&VectorVisitor::Vec6_tail)
		;
	}
	// only used with dim==6
	static CompatVec6* Vec6_fromElements(const Scalar& v0, const Scalar& v1, const Scalar& v2, const Scalar& v3, const Scalar& v4, const Scalar& v5){ CompatVec6* v(new CompatVec6); (*v)<<v0,v1,v2,v3,v4,v5; return v; }
	static CompatVec6* Vec6_fromHeadTail(const CompatVec3& head, const CompatVec3& tail){ CompatVec6* v(new CompatVec6); (*v)<<head,tail; return v; }
	static CompatVec3 Vec6_head(const CompatVec6& v){ return v.template head<3>(); }
	static CompatVec3 Vec6_tail(const CompatVec6& v){ return v.template tail<3>(); }

	// ctor for dynamic vectors
	template<typename VectorT2, class PyClass> static void visit_special_sizes(PyClass& cl, typename boost::enable_if_c<VectorT2::RowsAtCompileTime==Eigen::Dynamic>::type* dummy=0){
		cl
		.def("__init__",py::make_constructor(&VecX_fromList,py::default_call_policies(),(py::arg("vv"))))
		;
	}
	static CompatVecX* VecX_fromList(const std::vector<Scalar>& ii){ CompatVecX* v(new CompatVecX(ii.size())); for(size_t i=0; i<ii.size(); i++) (*v)[i]=ii[i]; return v; }


	static VectorT dyn_Ones(int size){ return VectorT::Ones(size); }
	static VectorT dyn_Zero(int size){ return VectorT::Zero(size); }
	static VectorT dyn_Random(int size){ return VectorT::Random(size); }
	static VectorT Unit(int ix){ IDX_CHECK(ix,(int)Dim); return VectorT::Unit(ix); }
	static VectorT dyn_Unit(int size,int ix){ IDX_CHECK(ix,size); return VectorT::Unit(size,ix); }
	static bool dyn(){ return Dim==Eigen::Dynamic; }
	static int __len__(){ assert(!dyn()); return Dim; }
	static int dyn__len__(const VectorT& self){ return self.size(); }
	static void resize(VectorT& self, int size){ self.resize(size); }
	static Scalar dot(const VectorT& self, const VectorT& other){ return self.dot(other); }
	static CompatMatrixT outer(const VectorT& self, const VectorT& other){ return self*other.transpose(); }
	static CompatMatrixT asDiagonal(const VectorT& self){ return self.asDiagonal(); }
	static Scalar get_item(const VectorT& self, int ix){ IDX_CHECK(ix,dyn()?(int)self.size():(int)Dim); return self[ix]; }
	static void set_item(VectorT& self, int ix, Scalar value){ IDX_CHECK(ix,dyn()?(int)self.size():(int)Dim); self[ix]=value; }
	struct VectorPickle: py::pickle_suite{
		static py::tuple getinitargs(const VectorT& x){
			// if this fails, add supported size to the switch below
			BOOST_STATIC_ASSERT(Dim==2 || Dim==3 || Dim==6 || Dim==Eigen::Dynamic);
			switch((int)Dim){
				case 2: return py::make_tuple(x[0],x[1]);
				case 3: return py::make_tuple(x[0],x[1],x[2]);
				case 6: return py::make_tuple(x[0],x[1],x[2],x[3],x[4],x[5]);
				default: return py::make_tuple(py::list(x));
			}
		};
	};
	public:
	static string __str__(const py::object& obj){
		std::ostringstream oss;
		const VectorT& self=py::extract<VectorT>(obj)();
		bool list=(Dim==Eigen::Dynamic && self.size()>0);
		oss<<object_class_name(obj)<<(list?"([":"(");
		Vector_data_stream(self,oss);
		oss<<(list?"])":")");
		return oss.str();
	};

	// not sure why this must be templated now?!
	template<typename VectorType>
	static void Vector_data_stream(const VectorType& self, std::ostringstream& oss, int pad=0){
		for(int i=0; i<self.size(); i++) oss<<(i==0?"":(((i%3)!=0 || pad>0)?",":", "))<<num_to_string(self.row(i/self.cols())[i%self.cols()],/*pad*/pad);
	}
};

template<typename MatrixT>
class MatrixVisitor: public py::def_visitor<MatrixVisitor<MatrixT> >{
	friend class def_visitor_access;
	typedef typename MatrixT::Scalar Scalar;
	typedef typename Eigen::Matrix<Scalar,MatrixT::RowsAtCompileTime,1> CompatVectorT;
	typedef Eigen::Matrix<Scalar,3,3> CompatMat3;
	typedef Eigen::Matrix<Scalar,3,1> CompatVec3;
	typedef Eigen::Matrix<Scalar,6,6> CompatMat6;
	typedef Eigen::Matrix<Scalar,6,1> CompatVec6;
	typedef Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> CompatMatX;
	typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1> CompatVecX;
	enum{Dim=MatrixT::RowsAtCompileTime};
	public:
	template<class PyClass>
	void visit(PyClass& cl) const {
		MatrixBaseVisitor<MatrixT>().visit(cl);
		cl
		.def_pickle(MatrixPickle())
		.def("__init__",py::make_constructor(&MatrixVisitor::fromDiagonal,py::default_call_policies(),(py::arg("diag"))))

		.def("determinant",&MatrixT::determinant,"Return matrix determinant.")
		.def("trace",&MatrixT::trace,"Return sum of diagonal elements.")
		.def("transpose",&MatrixVisitor::transpose,"Return transposed matrix.")
		.def("diagonal",&MatrixVisitor::diagonal,"Return diagonal as vector.")
		.def("row",&MatrixVisitor::row,py::arg("row"),"Return row as vector.")
		.def("col",&MatrixVisitor::col,py::arg("col"),"Return column as vector.")
		// matrix-matrix product
		.def("__mul__",&MatrixVisitor::__mul__).def("__imul__",&MatrixVisitor::__imul__)
		.def("__mul__",&MatrixVisitor::__mul__vec).def("__rmul__",&MatrixVisitor::__mul__vec)
		.def("__setitem__",&MatrixVisitor::set_row).def("__getitem__",&MatrixVisitor::get_row)
		.def("__setitem__",&MatrixVisitor::set_item).def("__getitem__",&MatrixVisitor::get_item)
		.def("__str__",&MatrixVisitor::__str__).def("__repr__",&MatrixVisitor::__str__)
		;
		visit_if_float<Scalar,PyClass>(cl);
		visit_fixed_or_dynamic<MatrixT,PyClass>(cl);
		visit_special_sizes<MatrixT,PyClass>(cl);

		//std::cerr<<"MatrixVisitor: "<<name()<<std::endl;

	}
	private:
	// for dynamic matrices
	template<typename MatrixT2, class PyClass> static void visit_fixed_or_dynamic(PyClass& cl, typename boost::enable_if_c<MatrixT2::RowsAtCompileTime==Eigen::Dynamic>::type* dummy = 0){
		cl
		.def("__len__",&MatrixVisitor::__len__)
		.def("resize",&MatrixVisitor::resize,"Change size of the matrix, keep values of elements which exist in the new matrix",(py::arg("rows"),py::arg("cols")))
		.def("Ones",&MatrixVisitor::dyn_Ones,(py::arg("rows"),py::arg("cols")),"Create matrix of given dimensions where all elements are set to 1.").staticmethod("Ones")
		.def("Zero",&MatrixVisitor::dyn_Zero,(py::arg("rows"),py::arg("cols")),"Create zero matrix of given dimensions").staticmethod("Zero")
		.def("Random",&MatrixVisitor::dyn_Random,(py::arg("rows"),py::arg("cols")),"Create matrix with given dimensions where all elements are set to number between 0 and 1 (uniformly-distributed).").staticmethod("Random")
		.def("Identity",&MatrixVisitor::dyn_Identity,(py::arg("rank")),"Create identity matrix with given rank (square).").staticmethod("Identity")
		;
	}
	// for fixed-size matrices
	template<typename MatrixT2, class PyClass> static void visit_fixed_or_dynamic(PyClass& cl, typename boost::disable_if_c<MatrixT2::RowsAtCompileTime==Eigen::Dynamic>::type* dummy = 0){
		cl.def("__len__",&MatrixVisitor::__len__).staticmethod("__len__")
		;
	}

	template<typename Scalar, class PyClass> static	void visit_if_float(PyClass& cl, typename boost::enable_if<boost::is_integral<Scalar> >::type* dummy = 0){ /* do nothing */ }
	template<typename Scalar, class PyClass> static void visit_if_float(PyClass& cl, typename boost::disable_if<boost::is_integral<Scalar> >::type* dummy = 0){
		cl
		// matrix-matrix division?!
		//.def("__div__",&MatrixBaseVisitor::__div__).def("__idiv__",&MatrixBaseVisitor::__idiv__)
		.def("inverse",&MatrixT::inverse,"Return inverted matrix.");
		// decompositions are only meaningful on non-complex numbers
		visit_if_decompositions_meaningful<Scalar,PyClass>(cl);
	}
	// for complex numbers, do nothing
	template<typename Scalar, class PyClass> static	void visit_if_decompositions_meaningful(PyClass& cl, typename boost::enable_if<boost::is_complex<Scalar> >::type* dummy = 0){ /* do nothing */ }
	// for non-complex numbers, define decompositions
	template<typename Scalar, class PyClass> static	void visit_if_decompositions_meaningful(PyClass& cl, typename boost::disable_if<boost::is_complex<Scalar> >::type* dummy = 0){
		cl
		.def("jacobiSVD",&MatrixVisitor::jacobiSVD,"Compute SVD decomposition of square matrix, retuns (U,S,V) such that self=U*S*V.transpose()")
		.def("svd",&MatrixVisitor::jacobiSVD,"Alias for :obj:`jacobiSVD`.")
		.def("computeUnitaryPositive",&MatrixVisitor::computeUnitaryPositive,"Compute polar decomposition (unitary matrix U and positive semi-definite symmetric matrix P such that self=U*P).")
		.def("polarDecomposition",&MatrixVisitor::computeUnitaryPositive,"Alias for :obj:`computeUnitaryPositive`.")
		.def("selfAdjointEigenDecomposition",&MatrixVisitor::selfAdjointEigenDecomposition,"Compute eigen (spectral) decomposition of symmetric matrix, returns (eigVecs,eigVals). eigVecs is orthogonal Matrix3 with columns ar normalized eigenvectors, eigVals is Vector3 with corresponding eigenvalues. self=eigVecs*diag(eigVals)*eigVecs.transpose().")
		.def("spectralDecomposition",&MatrixVisitor::selfAdjointEigenDecomposition,"Alias for :obj:`selfAdjointEigenDecomposition`.")
		;
	}

	// handle specific matrix sizes
	// 3x3
	template<typename MatT2, class PyClass> static void visit_special_sizes(PyClass& cl, typename boost::enable_if_c<MatT2::RowsAtCompileTime==3>::type* dummy=0){
		cl
		.def("__init__",py::make_constructor(&MatrixVisitor::Mat3_fromElements,py::default_call_policies(),(py::arg("m00"),py::arg("m01"),py::arg("m02"),py::arg("m10"),py::arg("m11"),py::arg("m12"),py::arg("m20"),py::arg("m21"),py::arg("m22"))))
		.def("__init__",py::make_constructor(&MatrixVisitor::Mat3_fromRows,py::default_call_policies(),(py::arg("r0"),py::arg("r1"),py::arg("r2"),py::arg("cols")=false)))
		;
	}
	static CompatMat3* Mat3_fromElements(const Scalar& m00, const Scalar& m01, const Scalar& m02, const Scalar& m10, const Scalar& m11, const Scalar& m12, const Scalar& m20, const Scalar& m21, const Scalar& m22){ CompatMat3* m(new CompatMat3); (*m)<<m00,m01,m02,m10,m11,m12,m20,m21,m22; return m; }
	static CompatMat3* Mat3_fromRows(const CompatVec3& l0, const CompatVec3& l1, const CompatVec3& l2, bool cols=false){ CompatMat3* m(new CompatMat3); if(cols){m->col(0)=l0; m->col(1)=l1; m->col(2)=l2; } else {m->row(0)=l0; m->row(1)=l1; m->row(2)=l2;} return m; }

	// 6x6
	template<typename MatT2, class PyClass> static void visit_special_sizes(PyClass& cl, typename boost::enable_if_c<MatT2::RowsAtCompileTime==6>::type* dummy=0){
		cl
		.def("__init__",py::make_constructor(&MatrixVisitor::Mat6_fromBlocks,py::default_call_policies(),(py::arg("ul"),py::arg("ur"),py::arg("ll"),py::arg("lr"))))
		.def("__init__",py::make_constructor(&MatrixVisitor::Mat6_fromRows,py::default_call_policies(),(py::arg("l0"),py::arg("l1"),py::arg("l2"),py::arg("l3"),py::arg("l4"),py::arg("l5"),py::arg("cols")=false)))
		/* 3x3 blocks */
			.def("ul",&MatrixVisitor::Mat6_ul,"Return upper-left 3x3 block")
			.def("ur",&MatrixVisitor::Mat6_ur,"Return upper-right 3x3 block")
			.def("ll",&MatrixVisitor::Mat6_ll,"Return lower-left 3x3 block")
			.def("lr",&MatrixVisitor::Mat6_lr,"Return lower-right 3x3 block")
		;
	}
	static CompatMat6* Mat6_fromBlocks(const CompatMat3& ul, const CompatMat3& ur, const CompatMat3& ll, const CompatMat3& lr){ CompatMat6* m(new CompatMat6); (*m)<<ul,ur,ll,lr; return m; }
	static CompatMat6* Mat6_fromRows(const CompatVec6& l0, const CompatVec6& l1, const CompatVec6& l2, const CompatVec6& l3, const CompatVec6& l4, const CompatVec6& l5, bool cols=false){ CompatMat6* m(new CompatMat6); if(cols){ m->col(0)=l0; m->col(1)=l1; m->col(2)=l2; m->col(3)=l3; m->col(4)=l4; m->col(5)=l5; } else { m->row(0)=l0; m->row(1)=l1; m->row(2)=l2; m->row(3)=l3; m->row(4)=l4; m->row(5)=l5; } return m; }
		// see http://stackoverflow.com/questions/20870648/error-in-seemingly-correct-template-code-whats-wrong#20870661
		// on why the template keyword is needed after the m.
		static CompatMat3 Mat6_ul(const CompatMat6& m){ return m.template topLeftCorner<3,3>(); }
		static CompatMat3 Mat6_ur(const CompatMat6& m){ return m.template topRightCorner<3,3>(); }
		static CompatMat3 Mat6_ll(const CompatMat6& m){ return m.template bottomLeftCorner<3,3>(); }
		static CompatMat3 Mat6_lr(const CompatMat6& m){ return m.template bottomRightCorner<3,3>(); }

	// XxX
	template<typename MatT2, class PyClass> static void visit_special_sizes(PyClass& cl, typename boost::enable_if_c<MatT2::RowsAtCompileTime==Eigen::Dynamic>::type* dummy=0){
		cl
		.def("__init__",py::make_constructor(&MatrixVisitor::MatX_fromRows,py::default_call_policies(),(py::arg("r0")=CompatVecX(),py::arg("r1")=CompatVecX(),py::arg("r2")=CompatVecX(),py::arg("r3")=CompatVecX(),py::arg("r4")=CompatVecX(),py::arg("r5")=CompatVecX(),py::arg("r6")=CompatVecX(),py::arg("r7")=CompatVecX(),py::arg("r8")=CompatVecX(),py::arg("r9")=CompatVecX(),py::arg("cols")=false)))
		.def("__init__",py::make_constructor(&MatrixVisitor::MatX_fromRowSeq,py::default_call_policies(),(py::arg("rows"),py::arg("cols")=false)))
		;
	}

	static CompatMatX* MatX_fromRows(const CompatVecX& r0, const CompatVecX& r1, const CompatVecX& r2, const CompatVecX& r3, const CompatVecX& r4, const CompatVecX& r5, const CompatVecX& r6, const CompatVecX& r7, const CompatVecX& r8, const CompatVecX& r9, bool setCols){
		/* check vector dimensions */ CompatVecX rr[]={r0,r1,r2,r3,r4,r5,r6,r7,r8,r9};
		int cols=-1, rows=-1;
		for(int i=0; i<10; i++){
			if(rows<0 && rr[i].size()==0) rows=i;
			if(rows>=0 && rr[i].size()>0) throw std::invalid_argument("Matrix6r: non-empty rows not allowed after first empty row, which marks end of the matrix.");
		}
		cols=(rows>0?rr[0].size():0);
		for(int i=1; i<rows; i++) if(rr[i].size()!=cols) throw std::invalid_argument(("Matrix6: all non-empty rows must have the same length (0th row has "+lexical_cast<string>(rr[0].size())+" items, "+lexical_cast<string>(i)+"th row has "+lexical_cast<string>(rr[i].size())+" items)").c_str());
		CompatMatX* m;
		m=setCols?new CompatMatX(cols,rows):new CompatMatX(rows,cols);
		for(int i=0; i<rows; i++){ if(setCols) m->col(i)=rr[i]; else m->row(i)=rr[i]; }
		return m;
	}
	static CompatMatX* MatX_fromRowSeq(const std::vector<CompatVecX>& rr, bool setCols){
		int rows=rr.size(),cols=rr.size()>0?rr[0].size():0;
		for(int i=1; i<rows; i++) if(rr[i].size()!=cols) throw std::invalid_argument(("MatrixX: all rows must have the same length."));
		CompatMatX* m;
		m=setCols?new CompatMatX(cols,rows):new CompatMatX(rows,cols);
		for(int i=0; i<rows; i++){ if(setCols) m->col(i)=rr[i]; else m->row(i)=rr[i]; }
		return m;
	};


	static MatrixT dyn_Ones(int rows, int cols){ return MatrixT::Ones(rows,cols); }
	static MatrixT dyn_Zero(int rows, int cols){ return MatrixT::Zero(rows,cols); }
	static MatrixT dyn_Random(int rows, int cols){ return MatrixT::Random(rows,cols); }
	static MatrixT dyn_Identity(int rows, int cols){ return MatrixT::Identity(rows,cols); }
	static MatrixT Identity(){ return MatrixT::Identity(); }
	static MatrixT transpose(const MatrixT& m){ return m.transpose(); }
	static CompatVectorT diagonal(const MatrixT& m){ return m.diagonal(); }
	static MatrixT* fromDiagonal(const CompatVectorT& d){ MatrixT* m(new MatrixT); *m=d.asDiagonal(); return m; }
	static typename MatrixT::Index __len__(MatrixT& a){ return a.rows(); }
	static void resize(MatrixT& self, int rows, int cols){ self.resize(rows,cols); }
	static CompatVectorT get_row(const MatrixT& a, int ix){ return a.row(ix); }
	static void set_row(MatrixT& a, int ix, const CompatVectorT& r){ a.row(ix)=r; }
	static Scalar get_item(const MatrixT& a, py::tuple _idx){ int idx[2]; int mx[2]={a.rows(),a.cols()}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); return a(idx[0],idx[1]); }
	static void set_item(MatrixT& a, py::tuple _idx, const Scalar& value){ int idx[2]; int mx[2]={a.rows(),a.cols()}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); a(idx[0],idx[1])=value; }

	static void __imul__(MatrixT& a, const MatrixT& b){ a*=b; };
	static MatrixT __mul__(const MatrixT& a, const MatrixT& b){ return a*b; }
	static CompatVectorT __mul__vec(const MatrixT& m, const CompatVectorT& v){ return m*v; }
	// float matrices only
	static MatrixT inverse(const MatrixT& m){ return m.inverse(); }
	static MatrixT __div__(const MatrixT& a, const MatrixT& b){ return a/b; }
	static void __idiv__(MatrixT& a, const MatrixT& b){ a/=b; };
	static CompatVectorT row(const MatrixT& m, int ix){ IDX_CHECK(ix,m.rows()); return m.row(ix); }
	static CompatVectorT col(const MatrixT& m, int ix){ IDX_CHECK(ix,m.cols()); return m.col(ix); }

	static void ensureSquare(const MatrixT& m){ if(m.rows()!=m.cols()) throw std::runtime_error("Matrix is not square."); }
	static py::tuple jacobiSVD(const MatrixT& in) {
		ensureSquare(in);
		Eigen::JacobiSVD<MatrixT> svd(in, Eigen::ComputeThinU | Eigen::ComputeThinV);
		return py::make_tuple(svd.matrixU(),svd.matrixV(),MatrixT(svd.singularValues().asDiagonal()));
	};
	// polar decomposition
	static py::tuple computeUnitaryPositive(const MatrixT& in) {
		ensureSquare(in);
		Eigen::JacobiSVD<MatrixT> svd(in, Eigen::ComputeThinU | Eigen::ComputeThinV);
		const MatrixT& u=svd.matrixU(); const MatrixT& v=svd.matrixV();
		MatrixT s=svd.singularValues().asDiagonal();
		return py::make_tuple(u*v.transpose(),v*s*v.transpose());
	}
	// eigen decomposition
	static py::tuple selfAdjointEigenDecomposition(const MatrixT& in) {
		ensureSquare(in);
		Eigen::SelfAdjointEigenSolver<MatrixT> a(in);
		return py::make_tuple(a.eigenvectors(),a.eigenvalues());
	}
	static bool dyn(){ return Dim==Eigen::Dynamic; }
	static string __str__(const py::object& obj){
		std::ostringstream oss;
		const MatrixT& m=py::extract<MatrixT>(obj)();
		oss<<object_class_name(obj)<<"(";
		bool wrap=((dyn() && m.rows()>1) || (!dyn() && m.rows()>3));
		// non-wrapping fixed-size: flat list of numbers, not rows as tuples (Matrix3)
		if(!dyn() && !wrap){
			VectorVisitor<CompatVectorT>::template Vector_data_stream<MatrixT>(m,oss,/*pad=*/0);
		} else {
			if(wrap) oss<<"\n";
			for(int r=0; r<m.rows(); r++){
				oss<<(wrap?"\t":"")<<"(";
				VectorVisitor<CompatVectorT>::template Vector_data_stream<CompatVectorT>(m.row(r),oss,/*pad=*/(wrap?7:0));
				oss<<")"<<(r<m.rows()-1?",":"")<<(wrap?"\n":"");
			}
		}
		oss<<")";
		return oss.str();
	}
	struct MatrixPickle: py::pickle_suite{
		static py::tuple getinitargs(const MatrixT& x){
			// if this fails, add supported size to the switch below
			BOOST_STATIC_ASSERT(Dim==2 || Dim==3 || Dim==6 || Dim==Eigen::Dynamic);
			switch((int)Dim){
				case 2: return py::make_tuple(x(0,0),x(0,1),x(1,0),x(1,1));
				case 3: return py::make_tuple(x(0,0),x(0,1),x(0,2),x(1,0),x(1,1),x(1,2),x(2,0),x(2,1),x(2,2));
				case 6: return py::make_tuple(x.row(0),x.row(1),x.row(2),x.row(3),x.row(4),x.row(5));
				// should return list of rows, which are VectorX
				default: return py::make_tuple(py::list(x));
			}
		};
	};
};


template<typename Box>
class AabbVisitor: public py::def_visitor<AabbVisitor<Box> >{
	friend class def_visitor_access;
	typedef typename Box::VectorType VectorType;
	typedef typename Box::Scalar Scalar;
	public:
	template <class PyClass>
	void visit(PyClass& cl) const {
		cl
		.def(py::init<Box>(py::arg("other")))
		.def(py::init<VectorType,VectorType>((py::arg("min"),py::arg("max"))))
		.def_pickle(BoxPickle())
		.def("volume",&Box::volume)
		.def("empty",&Box::isEmpty)
		.def("center",&AabbVisitor::center)
		.def("sizes",&AabbVisitor::sizes)
		.def("contains",&AabbVisitor::containsPt)
		.def("contains",&AabbVisitor::containsBox)
		// for the "in" operator
		.def("__contains__",&AabbVisitor::containsPt) 
		.def("__contains__",&AabbVisitor::containsBox)
		.def("extend",&AabbVisitor::extendPt)
		.def("extend",&AabbVisitor::extendBox)
		.def("clamp",&AabbVisitor::clamp)
		// return new objects
		.def("intersection",&Box::intersection)
		.def("merged",&Box::merged)
		// those return internal references, which is what we want (FIXME: this is not true, they return copies!!)
		.add_property("min",&AabbVisitor::min) 
		.add_property("max",&AabbVisitor::max)
		.def("__len__",&AabbVisitor::len).staticmethod("__len__")
		.def("__setitem__",&AabbVisitor::set_item).def("__getitem__",&AabbVisitor::get_item)
		.def("__setitem__",&AabbVisitor::set_minmax).def("__getitem__",&AabbVisitor::get_minmax)
		.def("__str__",&AabbVisitor::__str__).def("__repr__",&AabbVisitor::__str__)
		;
	};
	private:
	static bool containsPt(const Box& self, const VectorType& pt){ return self.contains(pt); }
	static bool containsBox(const Box& self, const Box& other){ return self.contains(other); }
	static void extendPt(Box& self, const VectorType& pt){ self.extend(pt); }
	static void extendBox(Box& self, const Box& other){ self.extend(other); }
	static void clamp(Box& self, const Box& other){ self.clamp(other); }
	static VectorType min(const Box& self){ return self.min(); }
	static VectorType max(const Box& self){ return self.max(); }
	static VectorType center(const Box& self){ return self.center(); }
	static VectorType sizes(const Box& self){ return self.sizes(); }
	struct BoxPickle: py::pickle_suite{
		static py::tuple getinitargs(const Box& x){ return py::make_tuple(x.min(),x.max()); }
	};
	static int len(){ return Box::AmbientDimAtCompileTime; }
	// getters and setters 
	static Scalar get_item(const Box& self, py::tuple _idx){ int idx[2]; int mx[2]={2,Box::AmbientDimAtCompileTime}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); if(idx[0]==0) return self.min()[idx[1]]; return self.max()[idx[1]]; }
	static void set_item(Box& self, py::tuple _idx, Scalar value){ int idx[2]; int mx[2]={2,Box::AmbientDimAtCompileTime}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); if(idx[0]==0) self.min()[idx[1]]=value; else self.max()[idx[1]]=value; }
	static VectorType get_minmax(const Box& self, int idx){ IDX_CHECK(idx,2); if(idx==0) return self.min(); return self.max(); }
	static void set_minmax(Box& self, int idx, const VectorType& value){ IDX_CHECK(idx,2); if(idx==0) self.min()=value; else self.max()=value; }
	static string __str__(const py::object& obj){
		const Box& self=py::extract<Box>(obj)();
		std::ostringstream oss; oss<<object_class_name(obj)<<"((";
		VectorVisitor<VectorType>::template Vector_data_stream<VectorType>(self.min(),oss);
		oss<<"), (";
		VectorVisitor<VectorType>::template Vector_data_stream<VectorType>(self.max(),oss);
		oss<<"))";
		return oss.str();
	}
};

template<typename QuaternionT>
class QuaternionVisitor:  public py::def_visitor<QuaternionVisitor<QuaternionT> >{
	typedef typename QuaternionT::Scalar Scalar;
	typedef Eigen::Matrix<Scalar,3,1> CompatVec3;
	typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1> CompatVecX;
	typedef Eigen::Matrix<Scalar,3,3> CompatMat3;
	typedef Eigen::AngleAxis<Scalar> AngleAxisT;
	public:
	template<class PyClass>
	void visit(PyClass& cl) const {
		cl
		.def("__init__",py::make_constructor(&QuaternionVisitor::fromAxisAngle,py::default_call_policies(),(py::arg("axis"),py::arg("angle"))))
		.def("__init__",py::make_constructor(&QuaternionVisitor::fromAngleAxis,py::default_call_policies(),(py::arg("angle"),py::arg("axis"))))
		.def("__init__",py::make_constructor(&QuaternionVisitor::fromTwoVectors,py::default_call_policies(),(py::arg("u"),py::arg("v"))))
		.def(py::init<Scalar,Scalar,Scalar,Scalar>((py::arg("w"),py::arg("x"),py::arg("y"),py::arg("z")),"Initialize from coefficients.\n\n.. note:: The order of coefficients is *w*, *x*, *y*, *z*. The [] operator numbers them differently, 0...4 for *x* *y* *z* *w*!"))
		.def(py::init<CompatMat3>((py::arg("rotMatrix")))) //,"Initialize from given rotation matrix.")
		.def(py::init<QuaternionT>((py::arg("other"))))
		.def_pickle(QuaternionPickle())
		// properties
		.add_static_property("Identity",&QuaternionVisitor::Identity)
		// methods
		.def("setFromTwoVectors",&QuaternionVisitor::setFromTwoVectors,((py::arg("u"),py::arg("v"))))
		.def("conjugate",&QuaternionT::conjugate)
		.def("toAxisAngle",&QuaternionVisitor::toAxisAngle)
		.def("toAngleAxis",&QuaternionVisitor::toAngleAxis)
		.def("toRotationMatrix",&QuaternionT::toRotationMatrix)
		.def("toRotationVector",&QuaternionVisitor::toRotationVector)
		.def("Rotate",&QuaternionVisitor::Rotate,((py::arg("v"))))
		.def("inverse",&QuaternionT::inverse)
		.def("norm",&QuaternionT::norm)
		.def("normalize",&QuaternionT::normalize)
		.def("normalized",&QuaternionT::normalized)
		// .def("random",&QuaternionVisitor::random,"Assign random orientation to the quaternion.")
		// operators
		.def(py::self * py::self)
		.def(py::self *= py::self)
		.def(py::self * py::other<CompatVec3>())
		.def("__eq__",&QuaternionVisitor::__eq__).def("__ne__",&QuaternionVisitor::__ne__)
		.def("__sub__",&QuaternionVisitor::__sub__) 
		// specials
		.def("__abs__",&QuaternionT::norm)
		.def("__len__",&QuaternionVisitor::__len__).staticmethod("__len__")
		.def("__setitem__",&QuaternionVisitor::__setitem__).def("__getitem__",&QuaternionVisitor::__getitem__)
		.def("__str__",&QuaternionVisitor::__str__).def("__repr__",&QuaternionVisitor::__str__)
		;
	}
	private:
	static QuaternionT* fromAxisAngle(const CompatVec3& axis, const Scalar& angle){ return new QuaternionT(AngleAxisT(angle,axis)); }
	static QuaternionT* fromAngleAxis(const Scalar& angle, const CompatVec3& axis){ return new QuaternionT(AngleAxisT(angle,axis)); }
	static QuaternionT* fromTwoVectors(const CompatVec3& u, const CompatVec3& v){ QuaternionT* q(new QuaternionT); q->setFromTwoVectors(u,v); return q; }

	struct QuaternionPickle: py::pickle_suite{static py::tuple getinitargs(const QuaternionT& x){ return py::make_tuple(x.w(),x.x(),x.y(),x.z());} };
	static QuaternionT Identity(){ return QuaternionT::Identity(); }
	static Vector3r Rotate(const QuaternionT& self, const Vector3r& u){ return self*u; }
	static py::tuple toAxisAngle(const QuaternionT& self){ AngleAxisT aa(self); return py::make_tuple(aa.axis(),aa.angle());}
	static py::tuple toAngleAxis(const QuaternionT& self){ AngleAxisT aa(self); return py::make_tuple(aa.angle(),aa.axis());}
	static CompatVec3 toRotationVector(const QuaternionT& self){ AngleAxisT aa(self); return aa.angle()*aa.axis();}
	static void setFromTwoVectors(QuaternionT& self, const Vector3r& u, const Vector3r& v){ self.setFromTwoVectors(u,v); }

	static bool __eq__(const QuaternionT& u, const QuaternionT& v){ return u.x()==v.x() && u.y()==v.y() && u.z()==v.z() && u.w()==v.w(); }
	static bool __ne__(const QuaternionT& u, const QuaternionT& v){ return !__eq__(u,v); }
	static CompatVecX __sub__(const QuaternionT& a, const QuaternionT& b){ CompatVecX r(4); r<<a.w()-b.w(),a.x()-b.x(),a.y()-b.y(),a.z()-b.z(); return r; }

	static Scalar __getitem__(const QuaternionT & self, int idx){ IDX_CHECK(idx,4); if(idx==0) return self.x(); if(idx==1) return self.y(); if(idx==2) return self.z(); return self.w(); }
	static void __setitem__(QuaternionT& self, int idx, Real value){ IDX_CHECK(idx,4); if(idx==0) self.x()=value; else if(idx==1) self.y()=value; else if(idx==2) self.z()=value; else if(idx==3) self.w()=value; }
	static string __str__(const py::object& obj){
		const QuaternionT& self=py::extract<QuaternionT>(obj)();
		AngleAxisT aa(self);
		return string(object_class_name(obj)+"((")+num_to_string(aa.axis()[0])+","+num_to_string(aa.axis()[1])+","+num_to_string(aa.axis()[2])+"),"+num_to_string(aa.angle())+")";
	}
	static int __len__(){return 4;}
};

#endif /* _NEW_VISITORS */


// this all goes away with visitors
#ifndef _NEW_VISITORS

	Real Quaternionr_get_item(const Quaternionr & self, int idx){ IDX_CHECK(idx,4); if(idx==0) return self.x(); if(idx==1) return self.y(); if(idx==2) return self.z(); return self.w(); }
	void Quaternionr_set_item(Quaternionr & self, int idx, Real value){ IDX_CHECK(idx,4);  if(idx==0) self.x()=value; else if(idx==1) self.y()=value; else if(idx==2) self.z()=value; else if(idx==3) self.w()=value; }
	string Quaternionr_str(const Quaternionr & self){ AngleAxisr aa(self); return string("Quaternion((")+doubleToShortest(aa.axis()[0])+","+doubleToShortest(aa.axis()[1])+","+doubleToShortest(aa.axis()[2])+"),"+doubleToShortest(aa.angle())+")";}
	int Quaternionr_len(){return 4;}
	struct Quaternionr_pickle: py::pickle_suite{static py::tuple getinitargs(const Quaternionr& x){ return py::make_tuple(x.w(),x.x(),x.y(),x.z());} };

	// quaternion subtraction: useful for relative error computation using (q1-q2).norm()
	static VectorXr Quaternionr__sub__Quaternionr(const Quaternionr& a, const Quaternionr& b){  VectorXr r(4); r<<a.w()-b.w(),a.x()-b.x(),a.y()-b.y(),a.z()-b.z(); return r; }

	template<typename Scalar> bool operator==(const Eigen::Quaternion<Scalar>& u, const Eigen::Quaternion<Scalar>& v){ return u.x()==v.x() && u.y()==v.y() && u.z()==v.z() && u.w()==v.w(); }
	template<typename Scalar> bool operator!=(const Eigen::Quaternion<Scalar>& u, const Eigen::Quaternion<Scalar>& v){ return !(u==v); }
	/*** wrap methods which return expression templates - those must be converted to the primitive type ***/
	static Vector3r Quaternionr_Rotate(Quaternionr& q, const Vector3r& u){ return q*u; }
	// supposed to return raw pointer (or auto_ptr), boost::python takes care of the lifetime management
	static Quaternionr* Quaternionr_fromAxisAngle(const Vector3r& axis, const Real angle){ return new Quaternionr(AngleAxisr(angle,axis)); }
	static Quaternionr* Quaternionr_fromAngleAxis(const Real angle, const Vector3r& axis){ return new Quaternionr(AngleAxisr(angle,axis)); }
	static py::tuple Quaternionr_toAxisAngle(const Quaternionr& self){ AngleAxisr aa(self); return py::make_tuple(aa.axis(),aa.angle());}
	static py::tuple Quaternionr_toAngleAxis(const Quaternionr& self){ AngleAxisr aa(self); return py::make_tuple(aa.angle(),aa.axis());}
	static Vector3r Quaternionr_toRotationVector(const Quaternionr& self){ AngleAxisr aa(self); return aa.angle()*aa.axis();}

	static Quaternionr Quaternionr_setFromTwoVectors(Quaternionr& q, const Vector3r& u, const Vector3r& v){ return q.setFromTwoVectors(u,v); }
	static bool Quaternionr__eq__(const Quaternionr& q1, const Quaternionr& q2){ return q1==q2; }
	static bool Quaternionr__ne__(const Quaternionr& q1, const Quaternionr& q2){ return q1!=q2; }
	static Quaternionr Quaternionr_Identity(){ return Quaternionr::Identity(); }

	// vector getters
	Real VectorXr_get_item(const VectorXr & self, int idx){ IDX_CHECK(idx,self.size()); return self[idx]; }
	Real Vector6r_get_item(const Vector6r & self, int idx){ IDX_CHECK(idx,6); return self[idx]; }
	int  Vector6i_get_item(const Vector6i & self, int idx){ IDX_CHECK(idx,6); return self[idx]; }
	Real Vector3r_get_item(const Vector3r & self, int idx){ IDX_CHECK(idx,3); return self[idx]; }
	int  Vector3i_get_item(const Vector3i & self, int idx){ IDX_CHECK(idx,3); return self[idx]; }
	Real Vector2r_get_item(const Vector2r & self, int idx){ IDX_CHECK(idx,2); return self[idx]; }
	int  Vector2i_get_item(const Vector2i & self, int idx){ IDX_CHECK(idx,2); return self[idx]; }

	// matrix getters
	Real     Matrix3r_get_item(Matrix3r & self, py::tuple _idx){ int idx[2]; int mx[2]={3,3}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); return self(idx[0],idx[1]); }
	Real     Matrix6r_get_item(Matrix6r & self, py::tuple _idx){ int idx[2]; int mx[2]={6,6}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); return self(idx[0],idx[1]); }
	Real     MatrixXr_get_item(MatrixXr & self, py::tuple _idx){ int idx[2]; int mx[2]={self.rows(),self.cols()}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); return self(idx[0],idx[1]); }
	Vector3r Matrix3r_get_row (Matrix3r & self, int idx){ IDX_CHECK(idx,3); return self.row(idx); }
	Vector6r Matrix6r_get_row (Matrix6r & self, int idx){ IDX_CHECK(idx,6); return self.row(idx); }
	VectorXr MatrixXr_get_row (MatrixXr & self, int idx){ IDX_CHECK(idx,self.rows()); return self.row(idx); }

	// vector setters
	void VectorXr_set_item(VectorXr & self, int idx, Real value){ IDX_CHECK(idx,self.size()); self[idx]=value; }
	void Vector6r_set_item(Vector6r & self, int idx, Real value){ IDX_CHECK(idx,6); self[idx]=value; }
	void Vector6i_set_item(Vector6i & self, int idx, int  value){ IDX_CHECK(idx,6); self[idx]=value; }
	void Vector3r_set_item(Vector3r & self, int idx, Real value){ IDX_CHECK(idx,3); self[idx]=value; }
	void Vector3i_set_item(Vector3i & self, int idx, int  value){ IDX_CHECK(idx,3); self[idx]=value; }
	void Vector2r_set_item(Vector2r & self, int idx, Real value){ IDX_CHECK(idx,2); self[idx]=value; }
	void Vector2i_set_item(Vector2i & self, int idx, int  value){ IDX_CHECK(idx,2); self[idx]=value; }

	// matrix setters
	void Matrix3r_set_item(Matrix3r & self, py::tuple _idx, Real value){ int idx[2]; int mx[2]={3,3}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); self(idx[0],idx[1])=value; }
	void Matrix6r_set_item(Matrix6r & self, py::tuple _idx, Real value){ int idx[2]; int mx[2]={6,6}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); self(idx[0],idx[1])=value; }
	void MatrixXr_set_item(MatrixXr & self, py::tuple _idx, Real value){ int idx[2]; int mx[2]={self.rows(),self.cols()}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); self(idx[0],idx[1])=value; }
	void Matrix3r_set_row (Matrix3r & self, int idx, const Vector3r& row){ IDX_CHECK(idx,3); self.row(idx)=row; }
	void MatrixXr_set_row (MatrixXr & self, int idx, const VectorXr& row){ IDX_CHECK(idx,self.rows()); self.row(idx)=row; }
	void Matrix6r_set_row (Matrix6r & self, int idx, const Vector6r& row){ IDX_CHECK(idx,6); self.row(idx)=row; }

	// resizers for dynamic matrices/vectors
	static void MatrixXr_resize(MatrixXr& m, int rows, int cols){ m.resize(rows,cols); }
	static void VectorXr_resize(VectorXr& v, int n){ v.resize(n); }

	// aligned boxes
	Real AlignedBox3r_get_item(AlignedBox3r & self, py::tuple _idx){ int idx[2]; int mx[2]={2,3}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); if(idx[0]==0) return self.min()[idx[1]]; return self.max()[idx[1]]; }
	void AlignedBox3r_set_item(AlignedBox3r & self, py::tuple _idx, Real value){ int idx[2]; int mx[2]={2,3}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); if(idx[0]==0) self.min()[idx[1]]=value; else self.max()[idx[1]]=value; }
	Vector3r AlignedBox3r_get_minmax(const AlignedBox3r & self, int idx){ IDX_CHECK(idx,2); if(idx==0) return self.min(); return self.max(); }
	void AlignedBox3r_set_minmax(AlignedBox3r & self, int idx, const Vector3r& value){ IDX_CHECK(idx,2); if(idx==0) self.min()=value; else self.max()=value; }
	Real AlignedBox2r_get_item(AlignedBox2r & self, py::tuple _idx){ int idx[2]; int mx[2]={2,2}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); if(idx[0]==0) return self.min()[idx[1]]; return self.max()[idx[1]]; }
	void AlignedBox2r_set_item(AlignedBox2r & self, py::tuple _idx, Real value){ int idx[2]; int mx[2]={2,2}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); if(idx[0]==0) self.min()[idx[1]]=value; else self.max()[idx[1]]=value; }
	Vector2r AlignedBox2r_get_minmax(const AlignedBox2r & self, int idx){ IDX_CHECK(idx,2); if(idx==0) return self.min(); return self.max(); }
	void AlignedBox2r_set_minmax(AlignedBox2r & self, int idx, const Vector2r& value){ IDX_CHECK(idx,2); if(idx==0) self.min()=value; else self.max()=value; }

	/*** I/O including pickling ***/
	string Vector6r_str(const Vector6r & self){ return string("Vector6(")+doubleToShortest(self[0])+","+doubleToShortest(self[1])+","+doubleToShortest(self[2])+", "+doubleToShortest(self[3])+","+doubleToShortest(self[4])+","+doubleToShortest(self[5])+")";}
	string VectorXr_str(const VectorXr & self){ std::ostringstream oss; oss<<"VectorX(["; for(int i=0; i<self.size(); i++) oss<<(i==0?"":((i%3)?",":", "))<<doubleToShortest(self[i]); oss<<"])"; return oss.str(); }
	string Vector6i_str(const Vector6i & self){ return string("Vector6i(")+lexical_cast<string>(self[0])+","+lexical_cast<string>(self[1])+","+lexical_cast<string>(self[2])+", "+lexical_cast<string>(self[3])+","+lexical_cast<string>(self[4])+","+lexical_cast<string>(self[5])+")";}
	string Vector3r_str(const Vector3r & self){ return string("Vector3(")+doubleToShortest(self[0])+","+doubleToShortest(self[1])+","+doubleToShortest(self[2])+")";}
	string Vector3i_str(const Vector3i & self){ return string("Vector3i(")+lexical_cast<string>(self[0])+","+lexical_cast<string>(self[1])+","+lexical_cast<string>(self[2])+")";}
	string Vector2r_str(const Vector2r & self){ return string("Vector2(")+doubleToShortest(self[0])+","+doubleToShortest(self[1])+")";}
	string Vector2i_str(const Vector2i & self){ return string("Vector2i(")+lexical_cast<string>(self[0])+","+lexical_cast<string>(self[1])+")";}
	string AlignedBox3r_str(const AlignedBox3r & self){ return string("AlignedBox3("+Vector3r_str(self.min())+","+Vector3r_str(self.max())+")"); }
	string AlignedBox2r_str(const AlignedBox2r & self){ return string("AlignedBox2("+Vector2r_str(self.min())+","+Vector2r_str(self.max())+")"); }
	string Matrix3r_str(const Matrix3r & self){ std::ostringstream oss; oss<<"Matrix3("; for(int i=0; i<3; i++) for(int j=0; j<3; j++) oss<<doubleToShortest(self(i,j))<<((i==2 && j==2)?")":",")<<((i<2 && j==2)?" ":""); return oss.str(); }
	string Matrix6r_str(const Matrix6r & self){ std::ostringstream oss; oss<<"Matrix6(\n"; for(int i=0; i<6; i++){ oss<<"\t("; for(int j=0; j<6; j++) oss<<doubleToShortest(self(i,j),/*pad*/7)<<(j==2?", ":(j==5?"),\n":",")); } oss<<")"; return oss.str(); }
	string MatrixXr_str(const MatrixXr & self){ std::ostringstream oss; bool wrap=self.rows()>1; oss<<"MatrixX("<<(wrap?"\n":""); for(int i=0; i<self.rows(); i++){ oss<<"\t("; for(int j=0; j<self.cols(); j++) oss<<doubleToShortest(self(i,j),/*pad*/7)<<((((j+1)%3)==0 && j!=self.cols()-1)?", ":(j==(self.cols()-1)?(wrap?"),\n":"),"):",")); } oss<<")"; return oss.str(); }

	int Matrix3r_len(){return 3;} // rows
	int Matrix6r_len(){return 6;} // rows
	int AlignedBox3r_len(){return 2;} // min, max
	int AlignedBox2r_len(){return 2;} // min, max
	int Vector6r_len(){return 6;}
	int Vector6i_len(){return 6;}
	int Vector3r_len(){return 3;}
	int Vector3i_len(){return 3;}
	int Vector2r_len(){return 2;}
	int Vector2i_len(){return 2;}
	static int VectorXr_len(const VectorXr& self){return self.size();}
	static int MatrixXr_len(const MatrixXr& self){return self.rows();} // rows

	// pickling support
	struct Matrix3r_pickle: py::pickle_suite{static py::tuple getinitargs(const Matrix3r& x){ return py::make_tuple(x(0,0),x(0,1),x(0,2),x(1,0),x(1,1),x(1,2),x(2,0),x(2,1),x(2,2));} };
	struct Matrix6r_pickle: py::pickle_suite{static py::tuple getinitargs(const Matrix6r& x){ return py::make_tuple(x.row(0),x.row(1),x.row(2),x.row(3),x.row(4),x.row(5));} };
	struct MatrixXr_pickle: py::pickle_suite{static py::tuple getinitargs(const MatrixXr& x){ std::vector<VectorXr> vv(x.rows()); for(int i=0; i<x.rows(); i++) vv[i]=x.row(i); return py::tuple(vv); } };
	struct AlignedBox3r_pickle: py::pickle_suite{static py::tuple getinitargs(const AlignedBox3r& x){ return py::make_tuple(x.min(),x.max()); } };
	struct AlignedBox2r_pickle: py::pickle_suite{static py::tuple getinitargs(const AlignedBox2r& x){ return py::make_tuple(x.min(),x.max()); } };
	struct Vector6r_pickle: py::pickle_suite{static py::tuple getinitargs(const Vector6r& x){ return py::make_tuple(x[0],x[1],x[2],x[3],x[4],x[5]);} };
	struct Vector6i_pickle: py::pickle_suite{static py::tuple getinitargs(const Vector6i& x){ return py::make_tuple(x[0],x[1],x[2],x[3],x[4],x[5]);} };
	struct Vector3r_pickle: py::pickle_suite{static py::tuple getinitargs(const Vector3r& x){ return py::make_tuple(x[0],x[1],x[2]);} };
	struct Vector3i_pickle: py::pickle_suite{static py::tuple getinitargs(const Vector3i& x){ return py::make_tuple(x[0],x[1],x[2]);} };
	struct Vector2r_pickle: py::pickle_suite{static py::tuple getinitargs(const Vector2r& x){ return py::make_tuple(x[0],x[1]);} };
	struct Vector2i_pickle: py::pickle_suite{static py::tuple getinitargs(const Vector2i& x){ return py::make_tuple(x[0],x[1]);} };
	struct VectorXr_pickle: py::pickle_suite{static py::tuple getinitargs(const VectorXr& x){ return py::make_tuple(py::list(x)); } };


	static Matrix3r* Matrix3r_fromElements(Real m00, Real m01, Real m02, Real m10, Real m11, Real m12, Real m20, Real m21, Real m22){ Matrix3r* m(new Matrix3r); (*m)<<m00,m01,m02,m10,m11,m12,m20,m21,m22; return m; }
	static Matrix3r* Matrix3r_fromRows(const Vector3r& l0, const Vector3r& l1, const Vector3r& l2, bool cols=false){ Matrix3r* m(new Matrix3r); if(cols){m->col(0)=l0; m->col(1)=l1; m->col(2)=l2; } else {m->row(0)=l0; m->row(1)=l1; m->row(2)=l2;} return m; }
	static Matrix3r* Matrix3r_fromDiagonal(const Vector3r& d){ Matrix3r* m(new Matrix3r); *m=d.asDiagonal(); return m; }
	static MatrixXr* MatrixXr_fromDiagonal(const VectorXr& d){ MatrixXr* m(new MatrixXr); *m=d.asDiagonal(); return m;}
	static MatrixXr* MatrixXr_fromRowSeq(const std::vector<VectorXr>& rr,bool setCols){
		int rows=rr.size(),cols=rr.size()>0?rr[0].size():0;
		for(int i=1; i<rows; i++) if(rr[i].size()!=cols) throw std::invalid_argument(("MatrixXr: all rows must have the same length."));
		MatrixXr* m;
		m=setCols?new MatrixXr(cols,rows):new MatrixXr(rows,cols);
		for(int i=0; i<rows; i++){ if(setCols) m->col(i)=rr[i]; else m->row(i)=rr[i]; }
		return m;
	};
	static MatrixXr* MatrixXr_fromRows(const VectorXr& r0, const VectorXr& r1, const VectorXr& r2, const VectorXr& r3, const VectorXr& r4, const VectorXr& r5, const VectorXr& r6, const VectorXr& r7, const VectorXr& r8, const VectorXr& r9, bool setCols){
		/* check vector dimensions */ VectorXr rr[]={r0,r1,r2,r3,r4,r5,r6,r7,r8,r9};
		int cols=-1, rows=-1;
		for(int i=0; i<10; i++){
			if(rows<0 && rr[i].size()==0) rows=i;
			if(rows>=0 && rr[i].size()>0) throw std::invalid_argument("Matrix6r: non-empty rows not allowed after first empty row, which marks end of the matrix.");
		}
		cols=(rows>0?rr[0].size():0);
		for(int i=1; i<rows; i++) if(rr[i].size()!=cols) throw std::invalid_argument(("Matrix6r: all non-empty rows must have the same length (0th row has "+lexical_cast<string>(rr[0].size())+" items, "+lexical_cast<string>(i)+"th row has "+lexical_cast<string>(rr[i].size())+" items)").c_str());
		MatrixXr* m;
		m=setCols?new MatrixXr(cols,rows):new MatrixXr(rows,cols);
		for(int i=0; i<rows; i++){ if(setCols) m->col(i)=rr[i]; else m->row(i)=rr[i]; }
		return m;
	}
	static Matrix6r* Matrix6r_fromBlocks(const Matrix3r& ul, const Matrix3r& ur, const Matrix3r& ll, const Matrix3r& lr){ Matrix6r* m(new Matrix6r); (*m)<<ul,ur,ll,lr; return m; }
	static Matrix6r* Matrix6r_fromRows(const Vector6r& l0, const Vector6r& l1, const Vector6r& l2, const Vector6r& l3, const Vector6r& l4, const Vector6r& l5, bool cols=false){ Matrix6r* m(new Matrix6r); if(cols){ m->col(0)=l0; m->col(1)=l1; m->col(2)=l2; m->col(3)=l3; m->col(4)=l4; m->col(5)=l5; } else { m->row(0)=l0; m->row(1)=l1; m->row(2)=l2; m->row(3)=l3; m->row(4)=l4; m->row(5)=l5; } return m; }
	static Matrix6r* Matrix6r_fromDiagonal(const Vector6r& d){ Matrix6r* m(new Matrix6r); *m=d.asDiagonal(); return m; }

	static Vector6r* Vector6r_fromElements(Real v0, Real v1, Real v2, Real v3, Real v4, Real v5){ Vector6r* v(new Vector6r); (*v)<<v0,v1,v2,v3,v4,v5; return v; }
	static Vector6r* Vector6r_fromHeadTail(const Vector3r& head, const Vector3r& tail){ Vector6r* ret(new Vector6r); ret->segment<3>(0)=head; ret->segment<3>(3)=tail; return ret; }
	static VectorXr* VectorXr_fromList(const std::vector<Real>& ii){ VectorXr* v(new VectorXr(ii.size())); for(size_t i=0; i<ii.size(); i++) (*v)[i]=ii[i]; return v; }
	static Vector6i* Vector6i_fromElements(int v0, int v1, int v2, int v3, int v4, int v5){ Vector6i* v(new Vector6i); (*v)<<v0,v1,v2,v3,v4,v5; return v; }

	static Matrix3r Matrix6r_ul(const Matrix6r& m){ return Matrix3r(m.block<3,3>(0,0)); }
	static Matrix3r Matrix6r_ur(const Matrix6r& m){ return Matrix3r(m.block<3,3>(0,3)); }
	static Matrix3r Matrix6r_ll(const Matrix6r& m){ return Matrix3r(m.block<3,3>(3,0)); }
	static Matrix3r Matrix6r_lr(const Matrix6r& m){ return Matrix3r(m.block<3,3>(3,3)); }
	// swizzles for Vector3r
	static Vector2r Vector3r_xy(const Vector3r& v){ return Vector2r(v[0],v[1]); }
	static Vector2r Vector3r_yx(const Vector3r& v){ return Vector2r(v[1],v[0]); }
	static Vector2r Vector3r_xz(const Vector3r& v){ return Vector2r(v[0],v[2]); }
	static Vector2r Vector3r_zx(const Vector3r& v){ return Vector2r(v[2],v[0]); }
	static Vector2r Vector3r_yz(const Vector3r& v){ return Vector2r(v[1],v[2]); }
	static Vector2r Vector3r_zy(const Vector3r& v){ return Vector2r(v[2],v[1]); }

	static Vector3r Vector3r_cross(const Vector3r& self, const Vector3r& v){ return self.cross(v); }
	static Vector3i Vector3i_cross(const Vector3i& self, const Vector3i& v){ return self.cross(v); }
	static Vector3r Vector6r_head(const Vector6r& self){ return self.head<3>(); }
	static Vector3r Vector6r_tail(const Vector6r& self){ return self.tail<3>(); }
	static Vector3i Vector6i_head(const Vector6i& self){ return self.head<3>(); }
	static Vector3i Vector6i_tail(const Vector6i& self){ return self.tail<3>(); }

	static Vector3r Matrix3r_diagonal(const Matrix3r& m){ return Vector3r(m.diagonal()); }
	static Vector6r Matrix6r_diagonal(const Matrix6r& m){ return Vector6r(m.diagonal()); }
	static VectorXr MatrixXr_diagonal(const MatrixXr& m){ return VectorXr(m.diagonal()); }
	static Vector3r Matrix3r_row(const Matrix3r& m, int ix){ IDX_CHECK(ix,3); return Vector3r(m.row(ix)); }
	static Vector3r Matrix3r_col(const Matrix3r& m, int ix){ IDX_CHECK(ix,3); return Vector3r(m.col(ix)); }
	static Vector6r Matrix6r_row(const Matrix6r& m, int ix){ IDX_CHECK(ix,6); return Vector6r(m.row(ix)); }
	static Vector6r Matrix6r_col(const Matrix6r& m, int ix){ IDX_CHECK(ix,6); return Vector6r(m.col(ix)); }
	static VectorXr MatrixXr_row(const MatrixXr& m, int ix){ IDX_CHECK(ix,m.rows()); return VectorXr(m.row(ix)); }
	static VectorXr MatrixXr_col(const MatrixXr& m, int ix){ IDX_CHECK(ix,m.cols()); return VectorXr(m.col(ix)); }
	static Vector3r AlignedBox3r_min(const AlignedBox3r& self){ return self.min(); }
	static Vector3r AlignedBox3r_max(const AlignedBox3r& self){ return self.max(); }
	static Vector3r AlignedBox3r_center(const AlignedBox3r& self){ return self.center(); }
	static Vector3r AlignedBox3r_size(const AlignedBox3r& self){ return self.sizes(); }

	static Vector2r AlignedBox2r_min(const AlignedBox2r& self){ return self.min(); }
	static Vector2r AlignedBox2r_max(const AlignedBox2r& self){ return self.max(); }
	static Vector2r AlignedBox2r_center(const AlignedBox2r& self){ return self.center(); }
	static Vector2r AlignedBox2r_size(const AlignedBox2r& self){ return self.sizes(); }

	/*** operations not defined in eigen ***/

	// cwise equality
	template<typename Scalar> bool operator==(const Eigen::Matrix<Scalar,2,1>& u, const Eigen::Matrix<Scalar,2,1>& v){ return u.x()==v.x() && u.y()==v.y(); }
	template<typename Scalar> bool operator!=(const Eigen::Matrix<Scalar,2,1>& u, const Eigen::Matrix<Scalar,2,1>& v){ return !(u==v); }
	template<typename Scalar> bool operator==(const Eigen::Matrix<Scalar,3,1>& u, const Eigen::Matrix<Scalar,3,1>& v){ return u.x()==v.x() && u.y()==v.y() && u.z()==v.z(); }
	template<typename Scalar> bool operator!=(const Eigen::Matrix<Scalar,3,1>& u, const Eigen::Matrix<Scalar,3,1>& v){ return !(u==v); }
	template<typename Scalar> bool operator==(const Eigen::Matrix<Scalar,6,1>& u, const Eigen::Matrix<Scalar,6,1>& v){ return u[0]==v[0] && u[1]==v[1] && u[2]==v[2] && u[3]==v[3] && u[4]==v[4] && u[5]==v[5]; }
	template<typename Scalar> bool operator!=(const Eigen::Matrix<Scalar,6,1>& u, const Eigen::Matrix<Scalar,6,1>& v){ return !(u==v); }
	template<typename Scalar> bool operator==(const Eigen::Matrix<Scalar,3,3>& m, const Eigen::Matrix<Scalar,3,3>& n){ for(int i=0;i<3;i++)for(int j=0;j<3;j++)if(m(i,j)!=n(i,j)) return false; return true; }
	template<typename Scalar> bool operator!=(const Eigen::Matrix<Scalar,3,3>& m, const Eigen::Matrix<Scalar,3,3>& n){ return !(m==n); }
	template<typename Scalar> bool operator==(const Eigen::Matrix<Scalar,6,6>& m, const Eigen::Matrix<Scalar,6,6>& n){ for(int i=0;i<6;i++)for(int j=0;j<6;j++)if(m(i,j)!=n(i,j)) return false; return true; }
	template<typename Scalar> bool operator!=(const Eigen::Matrix<Scalar,6,6>& m, const Eigen::Matrix<Scalar,6,6>& n){ return !(m==n); }
	template<typename Scalar> bool operator==(const Eigen::Matrix<Scalar,Eigen::Dynamic,1>& u, const Eigen::Matrix<Scalar,Eigen::Dynamic,1>& v){ if(u.rows()!=v.rows()) return false; for(int i=0; i<u.rows(); i++){ if(u[i]!=v[i]) return false; } return true; }
	template<typename Scalar> bool operator!=(const Eigen::Matrix<Scalar,Eigen::Dynamic,1>& u, const Eigen::Matrix<Scalar,Eigen::Dynamic,1>& v){ return !(u==v); }
	template<typename Scalar> bool operator==(const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& m, const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& n){ if(m.rows()!=n.rows() || m.cols()!=n.cols()) return false; for(int i=0;i<m.rows();i++)for(int j=0;j<m.cols();j++)if(m(i,j)!=n(i,j)) return false; return true; }
	template<typename Scalar> bool operator!=(const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& m, const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& n){ return !(m==n); }

	// remove almost-zero elements
	template<typename MatrixType> MatrixType Matrix_pruned(const MatrixType& obj, typename MatrixType::Scalar absTol=1e-6){ MatrixType ret(MatrixType::Zero(obj.rows(),obj.cols())); for(int i=0;i<obj.rows();i++){ for(int j=0;j<obj.cols();j++){ if(std::abs(obj(i,j))>absTol  && obj(i,j)!=-0) ret(i,j)=obj(i,j); } } return ret; }
	// element sum
	template<typename MatrixType> typename MatrixType::Scalar Matrix_sum(const MatrixType& obj){ return obj.sum(); }
	// outer product of two vectors
	template<typename VecTypeA, typename VecTypeB> typename Eigen::Matrix<typename VecTypeA::Scalar,VecTypeA::RowsAtCompileTime,VecTypeB::RowsAtCompileTime> Vector_outer(const VecTypeA& a, const VecTypeB& b){ return a*b.transpose(); }
	// maximum abs value coefficient
	template<typename MatrixType> typename MatrixType::Scalar Matrix_maxAbsCoeff(const MatrixType& obj){ return Eigen::Array<typename MatrixType::Scalar,MatrixType::RowsAtCompileTime,MatrixType::ColsAtCompileTime>(obj).abs().maxCoeff(); }
	static Real Vector3r_dot(const Vector3r& self, const Vector3r& v){ return self.dot(v); }
	static Real Vector3i_dot(const Vector3i& self, const Vector3i& v){ return self.dot(v); }
	static Real Vector2r_dot(const Vector2r& self, const Vector2r& v){ return self.dot(v); }
	static Real Vector2i_dot(const Vector2i& self, const Vector2i& v){ return self.dot(v); }

	static bool AlignedBox3r_containsPt(const AlignedBox3r& self, const Vector3r& v){ return self.contains(v); }
	static bool AlignedBox3r_containsBox(const AlignedBox3r& self, const AlignedBox3r& b){ return self.contains(b); }
	static void AlignedBox3r_extendPt(AlignedBox3r& self, const Vector3r& v){ self.extend(v); }
	static void AlignedBox3r_extendBox(AlignedBox3r& self, const AlignedBox3r& b){ self.extend(b); }
	static void AlignedBox3r_clamp(AlignedBox3r& self, const AlignedBox3r& b){ self.clamp(b); }

	static bool AlignedBox2r_containsPt(const AlignedBox2r& self, const Vector2r& v){ return self.contains(v); }
	static bool AlignedBox2r_containsBox(const AlignedBox2r& self, const AlignedBox2r& b){ return self.contains(b); }
	static void AlignedBox2r_extendPt(AlignedBox2r& self, const Vector2r& v){ self.extend(v); }
	static void AlignedBox2r_extendBox(AlignedBox2r& self, const AlignedBox2r& b){ self.extend(b); }
	static void AlignedBox2r_clamp(AlignedBox2r& self, const AlignedBox2r& b){ self.clamp(b); }
	template<typename VT> VT Vector_Unit(int ax){ IDX_CHECK(ax,VT::RowsAtCompileTime); return VT::Unit(ax); }
	template<typename VT> Eigen::Matrix<typename VT::Scalar,VT::RowsAtCompileTime,VT::RowsAtCompileTime>
	Vector_asDiagonal(const VT& self){ return self.asDiagonal(); }

	// svd decomposition
	/****** matrix decompositions ******/
	template<typename MatrixType>
	py::tuple Matrix_jacobiSVD(const MatrixType& in) {
		Eigen::JacobiSVD<MatrixType> svd(in, Eigen::ComputeThinU | Eigen::ComputeThinV);
		return py::make_tuple(svd.matrixU(),svd.matrixV(),MatrixType(svd.singularValues().asDiagonal()));
	}

	// polar decomposition
	template<typename MatrixType>
	py::tuple Matrix_computeUnitaryPositive(const MatrixType& in) {
		assert(in.rows()==in.cols());
		MatrixType u, s, v, retU, retP;
		Eigen::JacobiSVD<MatrixType> svd(in, Eigen::ComputeThinU | Eigen::ComputeThinV);
		u = svd.matrixU();
		v = svd.matrixV();
		s = svd.singularValues().asDiagonal();
		retU = u*v.transpose();
		retP = v*s*v.transpose();
		return py::make_tuple(retU,retP);
	}
		
	// eigen decomposition
	template<typename MatrixType>
	py::tuple Matrix_selfAdjointEigenDecomposition(const MatrixType& in) {
		Eigen::SelfAdjointEigenSolver<MatrixType> a(in);
		return py::make_tuple(a.eigenvectors(),a.eigenvalues());
	}

	// those are again needed to wrap the result
	#define EIG_WRAP_METH0_ROWS_COLS(klass,meth) static klass klass##_##meth(int rows, int cols){ return klass::meth(rows,cols); }
	#define EIG_WRAP_METH0_SIZE(klass,meth) static klass klass##_##meth(int size){ return klass::meth(size); }

	EIG_WRAP_METH0_ROWS_COLS(MatrixXr,Zero)
	EIG_WRAP_METH0_ROWS_COLS(MatrixXr,Ones)
	EIG_WRAP_METH0_ROWS_COLS(MatrixXr,Identity)
	EIG_WRAP_METH0_ROWS_COLS(MatrixXr,Random)
	EIG_WRAP_METH0_SIZE(VectorXr,Zero)
	EIG_WRAP_METH0_SIZE(VectorXr,Ones)
	EIG_WRAP_METH0_SIZE(VectorXr,Random)

	#define EIG_WRAP_METH1(klass,meth) static klass klass##_##meth(const klass& self){ return self.meth(); }
	#define EIG_WRAP_METH0(klass,meth) static klass klass##_##meth(){ return klass().meth(); }

	// this one is needed with visitors
	// EIG_WRAP_METH0(Quaternionr,Identity);

	EIG_WRAP_METH1(Matrix3r,transpose);
	EIG_WRAP_METH1(Matrix3r,inverse);
	EIG_WRAP_METH1(Matrix6r,transpose);
	EIG_WRAP_METH1(Matrix6r,inverse);
	EIG_WRAP_METH1(MatrixXr,transpose);
	EIG_WRAP_METH1(MatrixXr,inverse);

	EIG_WRAP_METH0(Matrix3r,Zero);
	EIG_WRAP_METH0(Matrix3r,Identity);
	EIG_WRAP_METH0(Matrix3r,Ones);
	EIG_WRAP_METH0(Matrix3r,Random);
	EIG_WRAP_METH0(Matrix6r,Zero);
	EIG_WRAP_METH0(Matrix6r,Identity);
	EIG_WRAP_METH0(Matrix6r,Ones);
	EIG_WRAP_METH0(Matrix6r,Random);
	EIG_WRAP_METH0(Vector6r,Zero); EIG_WRAP_METH0(Vector6r,Ones); EIG_WRAP_METH0(Vector6r,Random);
	EIG_WRAP_METH0(Vector6i,Zero); EIG_WRAP_METH0(Vector6i,Ones);
	EIG_WRAP_METH0(Vector3r,Zero); EIG_WRAP_METH0(Vector3r,UnitX); EIG_WRAP_METH0(Vector3r,UnitY); EIG_WRAP_METH0(Vector3r,UnitZ); EIG_WRAP_METH0(Vector3r,Ones); EIG_WRAP_METH0(Vector3r,Random);
	EIG_WRAP_METH0(Vector3i,Zero); EIG_WRAP_METH0(Vector3i,UnitX); EIG_WRAP_METH0(Vector3i,UnitY); EIG_WRAP_METH0(Vector3i,UnitZ); EIG_WRAP_METH0(Vector3i,Ones); 
	EIG_WRAP_METH0(Vector2r,Zero); EIG_WRAP_METH0(Vector2r,UnitX); EIG_WRAP_METH0(Vector2r,UnitY); EIG_WRAP_METH0(Vector2r,Ones); EIG_WRAP_METH0(Vector2r,Random);
	EIG_WRAP_METH0(Vector2i,Zero); EIG_WRAP_METH0(Vector2i,UnitX); EIG_WRAP_METH0(Vector2i,UnitY); EIG_WRAP_METH0(Vector2i,Ones);

	/*** expose operators ***/
	#define EIG_OP1_RET(ret,klass,op,sym) ret klass##op(const klass& self){ return (sym self).eval();}
	// for return type same as op1 type
	#define EIG_OP1(klass,op,sym) EIG_OP1_RET(klass,klass,op,sym)
	#define EIG_OP2_RET(ret,klass,op,sym,klass2) ret klass##op##klass2(const klass& self, const klass2& other){ return (self sym other).eval(); }
	// for return type same as op1 type (most of them)
	#define EIG_OP2(klass,op,sym,klass2) EIG_OP2_RET(klass,klass,op,sym,klass2)
	#define EIG_OP2_INPLACE(klass,op,sym,klass2) klass klass##op##klass2(klass& self, const klass2& other){ self sym other; return self; }

	EIG_OP1(Matrix3r,__neg__,-)
	EIG_OP2(Matrix3r,__add__,+,Matrix3r) EIG_OP2_INPLACE(Matrix3r,__iadd__,+=,Matrix3r)
	EIG_OP2(Matrix3r,__sub__,-,Matrix3r) EIG_OP2_INPLACE(Matrix3r,__isub__,-=,Matrix3r)
	EIG_OP2(Matrix3r,__mul__,*,Real) EIG_OP2(Matrix3r,__rmul__,*,Real) EIG_OP2_INPLACE(Matrix3r,__imul__,*=,Real)
	EIG_OP2(Matrix3r,__mul__,*,int) EIG_OP2(Matrix3r,__rmul__,*,int) EIG_OP2_INPLACE(Matrix3r,__imul__,*=,int)
	EIG_OP2_RET(Vector3r,Matrix3r,__mul__,*,Vector3r) EIG_OP2_RET(Vector3r,Matrix3r,__rmul__,*,Vector3r)
	EIG_OP2(Matrix3r,__mul__,*,Matrix3r) EIG_OP2_INPLACE(Matrix3r,__imul__,*=,Matrix3r)
	EIG_OP2(Matrix3r,__div__,/,Real) EIG_OP2_INPLACE(Matrix3r,__idiv__,/=,Real)
	EIG_OP2(Matrix3r,__div__,/,int) EIG_OP2_INPLACE(Matrix3r,__idiv__,/=,int)

	EIG_OP1(Matrix6r,__neg__,-)
	EIG_OP2(Matrix6r,__add__,+,Matrix6r) EIG_OP2_INPLACE(Matrix6r,__iadd__,+=,Matrix6r)
	EIG_OP2(Matrix6r,__sub__,-,Matrix6r) EIG_OP2_INPLACE(Matrix6r,__isub__,-=,Matrix6r)
	EIG_OP2(Matrix6r,__mul__,*,Real) EIG_OP2(Matrix6r,__rmul__,*,Real) EIG_OP2_INPLACE(Matrix6r,__imul__,*=,Real)
	EIG_OP2(Matrix6r,__mul__,*,int) EIG_OP2(Matrix6r,__rmul__,*,int) EIG_OP2_INPLACE(Matrix6r,__imul__,*=,int)
	EIG_OP2_RET(Vector6r,Matrix6r,__mul__,*,Vector6r) EIG_OP2_RET(Vector6r,Matrix6r,__rmul__,*,Vector6r)
	EIG_OP2(Matrix6r,__mul__,*,Matrix6r) EIG_OP2_INPLACE(Matrix6r,__imul__,*=,Matrix6r)
	EIG_OP2(Matrix6r,__div__,/,Real) EIG_OP2_INPLACE(Matrix6r,__idiv__,/=,Real)
	EIG_OP2(Matrix6r,__div__,/,int) EIG_OP2_INPLACE(Matrix6r,__idiv__,/=,int)

	EIG_OP1(MatrixXr,__neg__,-)
	EIG_OP2(MatrixXr,__add__,+,MatrixXr) EIG_OP2_INPLACE(MatrixXr,__iadd__,+=,MatrixXr)
	EIG_OP2(MatrixXr,__sub__,-,MatrixXr) EIG_OP2_INPLACE(MatrixXr,__isub__,-=,MatrixXr)
	EIG_OP2(MatrixXr,__mul__,*,Real) EIG_OP2(MatrixXr,__rmul__,*,Real) EIG_OP2_INPLACE(MatrixXr,__imul__,*=,Real)
	EIG_OP2(MatrixXr,__mul__,*,int) EIG_OP2(MatrixXr,__rmul__,*,int) EIG_OP2_INPLACE(MatrixXr,__imul__,*=,int)
	EIG_OP2_RET(VectorXr,MatrixXr,__mul__,*,VectorXr) EIG_OP2_RET(VectorXr,MatrixXr,__rmul__,*,VectorXr)
	EIG_OP2(MatrixXr,__mul__,*,MatrixXr) EIG_OP2_INPLACE(MatrixXr,__imul__,*=,MatrixXr)
	EIG_OP2(MatrixXr,__div__,/,Real) EIG_OP2_INPLACE(MatrixXr,__idiv__,/=,Real)
	EIG_OP2(MatrixXr,__div__,/,int) EIG_OP2_INPLACE(MatrixXr,__idiv__,/=,int)

	EIG_OP1(VectorXr,__neg__,-);
	EIG_OP2(VectorXr,__add__,+,VectorXr); EIG_OP2_INPLACE(VectorXr,__iadd__,+=,VectorXr)
	EIG_OP2(VectorXr,__sub__,-,VectorXr); EIG_OP2_INPLACE(VectorXr,__isub__,-=,VectorXr)
	EIG_OP2(VectorXr,__mul__,*,Real) EIG_OP2(VectorXr,__rmul__,*,Real) EIG_OP2_INPLACE(VectorXr,__imul__,*=,Real) EIG_OP2(VectorXr,__div__,/,Real) EIG_OP2_INPLACE(VectorXr,__idiv__,/=,Real)
	EIG_OP2(VectorXr,__mul__,*,int) EIG_OP2(VectorXr,__rmul__,*,int) EIG_OP2_INPLACE(VectorXr,__imul__,*=,int) EIG_OP2(VectorXr,__div__,/,int) EIG_OP2_INPLACE(VectorXr,__idiv__,/=,int)


	EIG_OP1(Vector6r,__neg__,-);
	EIG_OP2(Vector6r,__add__,+,Vector6r); EIG_OP2_INPLACE(Vector6r,__iadd__,+=,Vector6r)
	EIG_OP2(Vector6r,__sub__,-,Vector6r); EIG_OP2_INPLACE(Vector6r,__isub__,-=,Vector6r)
	EIG_OP2(Vector6r,__mul__,*,Real) EIG_OP2(Vector6r,__rmul__,*,Real) EIG_OP2_INPLACE(Vector6r,__imul__,*=,Real) EIG_OP2(Vector6r,__div__,/,Real) EIG_OP2_INPLACE(Vector6r,__idiv__,/=,Real)
	EIG_OP2(Vector6r,__mul__,*,int) EIG_OP2(Vector6r,__rmul__,*,int) EIG_OP2_INPLACE(Vector6r,__imul__,*=,int) EIG_OP2(Vector6r,__div__,/,int) EIG_OP2_INPLACE(Vector6r,__idiv__,/=,int)

	EIG_OP1(Vector6i,__neg__,-);
	EIG_OP2(Vector6i,__add__,+,Vector6i); EIG_OP2_INPLACE(Vector6i,__iadd__,+=,Vector6i)
	EIG_OP2(Vector6i,__sub__,-,Vector6i); EIG_OP2_INPLACE(Vector6i,__isub__,-=,Vector6i)
	EIG_OP2(Vector6i,__mul__,*,int) EIG_OP2(Vector6i,__rmul__,*,int) EIG_OP2_INPLACE(Vector6i,__imul__,*=,int) EIG_OP2(Vector6i,__div__,/,int) EIG_OP2_INPLACE(Vector6i,__idiv__,/=,int)

	EIG_OP1(Vector3r,__neg__,-);
	EIG_OP2(Vector3r,__add__,+,Vector3r); EIG_OP2_INPLACE(Vector3r,__iadd__,+=,Vector3r)
	EIG_OP2(Vector3r,__sub__,-,Vector3r); EIG_OP2_INPLACE(Vector3r,__isub__,-=,Vector3r)
	EIG_OP2(Vector3r,__mul__,*,Real) EIG_OP2(Vector3r,__rmul__,*,Real) EIG_OP2_INPLACE(Vector3r,__imul__,*=,Real) EIG_OP2(Vector3r,__div__,/,Real) EIG_OP2_INPLACE(Vector3r,__idiv__,/=,Real)
	EIG_OP2(Vector3r,__mul__,*,int) EIG_OP2(Vector3r,__rmul__,*,int) EIG_OP2_INPLACE(Vector3r,__imul__,*=,int) EIG_OP2(Vector3r,__div__,/,int) EIG_OP2_INPLACE(Vector3r,__idiv__,/=,int)

	EIG_OP1(Vector3i,__neg__,-);
	EIG_OP2(Vector3i,__add__,+,Vector3i); EIG_OP2_INPLACE(Vector3i,__iadd__,+=,Vector3i)
	EIG_OP2(Vector3i,__sub__,-,Vector3i); EIG_OP2_INPLACE(Vector3i,__isub__,-=,Vector3i)
	EIG_OP2(Vector3i,__mul__,*,int) EIG_OP2(Vector3i,__rmul__,*,int)  EIG_OP2_INPLACE(Vector3i,__imul__,*=,int)

	EIG_OP1(Vector2r,__neg__,-);
	EIG_OP2(Vector2r,__add__,+,Vector2r); EIG_OP2_INPLACE(Vector2r,__iadd__,+=,Vector2r)
	EIG_OP2(Vector2r,__sub__,-,Vector2r); EIG_OP2_INPLACE(Vector2r,__isub__,-=,Vector2r)
	EIG_OP2(Vector2r,__mul__,*,Real) EIG_OP2(Vector2r,__rmul__,*,Real) EIG_OP2_INPLACE(Vector2r,__imul__,*=,Real) EIG_OP2(Vector2r,__div__,/,Real) EIG_OP2_INPLACE(Vector2r,__idiv__,/=,Real)
	EIG_OP2(Vector2r,__mul__,*,int) EIG_OP2(Vector2r,__rmul__,*,int) EIG_OP2_INPLACE(Vector2r,__imul__,*=,int) EIG_OP2(Vector2r,__div__,/,int) EIG_OP2_INPLACE(Vector2r,__idiv__,/=,int)

	EIG_OP1(Vector2i,__neg__,-);
	EIG_OP2(Vector2i,__add__,+,Vector2i); EIG_OP2_INPLACE(Vector2i,__iadd__,+=,Vector2i)
	EIG_OP2(Vector2i,__sub__,-,Vector2i); EIG_OP2_INPLACE(Vector2i,__isub__,-=,Vector2i)
	EIG_OP2(Vector2i,__mul__,*,int)  EIG_OP2_INPLACE(Vector2i,__imul__,*=,int) EIG_OP2(Vector2i,__rmul__,*,int)
#endif /* _NEW_VISITORS */

BOOST_PYTHON_MODULE(minieigen){
	py::scope().attr("__doc__")="miniEigen is wrapper for a small part of the `Eigen <http://eigen.tuxfamily.org>`_ library. Refer to its documentation for details. All classes in this module support pickling.";

	py::docstring_options docopt;
	docopt.enable_all();
	docopt.disable_cpp_signatures();

	custom_VectorAnyAny_from_sequence<VectorXr>();
	custom_VectorAnyAny_from_sequence<Vector6r>();
	custom_VectorAnyAny_from_sequence<Vector6i>();
	custom_VectorAnyAny_from_sequence<Vector3r>();
	custom_VectorAnyAny_from_sequence<Vector3i>();
	custom_VectorAnyAny_from_sequence<Vector2r>();
	custom_VectorAnyAny_from_sequence<Vector2i>();
	custom_alignedBoxNr_from_seq<2>();
	custom_alignedBoxNr_from_seq<3>();
	custom_Quaternionr_from_axisAngle_or_angleAxis();

	custom_MatrixAnyAny_from_sequence<Matrix3r>();
	custom_MatrixAnyAny_from_sequence<Matrix6r>();
	custom_MatrixAnyAny_from_sequence<MatrixXr>();

	#ifdef _COMPLEX_SUPPORT
		custom_VectorAnyAny_from_sequence<Vector2cr>();
		custom_VectorAnyAny_from_sequence<Vector3cr>();
		custom_VectorAnyAny_from_sequence<Vector6cr>();
		custom_VectorAnyAny_from_sequence<VectorXcr>();
		custom_MatrixAnyAny_from_sequence<Matrix3cr>();
		custom_MatrixAnyAny_from_sequence<Matrix6cr>();
		custom_MatrixAnyAny_from_sequence<MatrixXcr>();

		py::class_<Vector2cr>("Vector2c","/*TODO*/",py::init<>()).def(VectorVisitor<Vector2cr>());
		py::class_<Vector3cr>("Vector3c","/*TODO*/",py::init<>()).def(VectorVisitor<Vector3cr>());
		py::class_<Vector6cr>("Vector6c","/*TODO*/",py::init<>()).def(VectorVisitor<Vector6cr>());
		py::class_<VectorXcr>("VectorXc","/*TODO*/",py::init<>()).def(VectorVisitor<VectorXcr>());

		py::class_<Matrix3cr>("Matrix3c","/*TODO*/",py::init<>()).def(MatrixVisitor<Matrix3cr>());
		py::class_<Matrix6cr>("Matrix6c","/*TODO*/",py::init<>()).def(MatrixVisitor<Matrix6cr>());
		py::class_<MatrixXcr>("MatrixXc","/*TODO*/",py::init<>()).def(MatrixVisitor<MatrixXcr>());
		;
	#endif

	py::class_<Quaternionr>("Quaternion","Quaternion representing rotation.\n\nSupported operations (``q`` is a Quaternion, ``v`` is a Vector3): ``q*q`` (rotation composition), ``q*=q``, ``q*v`` (rotating ``v`` by ``q``), ``q==q``, ``q!=q``.\n\nStatic attributes: ``Identity``.",py::init<>())
		#ifdef _NEW_VISITORS
			.def(QuaternionVisitor<Quaternionr>())
		#else
			.def("__init__",py::make_constructor(&Quaternionr_fromAxisAngle,py::default_call_policies(),(py::arg("axis"),py::arg("angle"))))
			.def("__init__",py::make_constructor(&Quaternionr_fromAngleAxis,py::default_call_policies(),(py::arg("angle"),py::arg("axis"))))
			.def(py::init<Real,Real,Real,Real>((py::arg("w"),py::arg("x"),py::arg("y"),py::arg("z")),"Initialize from coefficients.\n\n.. note:: The order of coefficients is *w*, *x*, *y*, *z*. The [] operator numbers them differently, 0...4 for *x* *y* *z* *w*!"))
			.def(py::init<Matrix3r>((py::arg("rotMatrix")))) //,"Initialize from given rotation matrix.")
			.def(py::init<Quaternionr>((py::arg("other"))))
			.def_pickle(Quaternionr_pickle())
			// properties
			.add_static_property("Identity",&Quaternionr_Identity)
			// methods
			.def("setFromTwoVectors",&Quaternionr_setFromTwoVectors,((py::arg("u"),py::arg("v"))))
			.def("conjugate",&Quaternionr::conjugate)
			.def("toAxisAngle",&Quaternionr_toAxisAngle).def("toAngleAxis",&Quaternionr_toAngleAxis)
			.def("toRotationMatrix",&Quaternionr::toRotationMatrix)
			.def("toRotationVector",&Quaternionr_toRotationVector)
			.def("Rotate",&Quaternionr_Rotate,((py::arg("v"))))
			.def("inverse",&Quaternionr::inverse)
			.def("norm",&Quaternionr::norm)
			.def("normalize",&Quaternionr::normalize)
			.def("normalized",&Quaternionr::normalized)
			#if 0
				.def("random",&Quaternionr_random,"Assign random orientation to the quaternion.")
			#endif
			// operators
			.def(py::self * py::self)
			.def(py::self *= py::self)
			.def(py::self * py::other<Vector3r>())
			.def("__eq__",&Quaternionr__eq__).def("__ne__",&Quaternionr__ne__)
			.def("__sub__",&Quaternionr__sub__Quaternionr) 
			// specials
			.def("__abs__",&Quaternionr::norm)
			.def("__len__",&Quaternionr_len).staticmethod("__len__")
			.def("__setitem__",&Quaternionr_set_item).def("__getitem__",&Quaternionr_get_item)
			.def("__str__",&Quaternionr_str).def("__repr__",&Quaternionr_str)
		#endif
	;



	py::class_<Matrix3r>("Matrix3","3x3 float matrix.\n\nSupported operations (``m`` is a Matrix3, ``f`` if a float/int, ``v`` is a Vector3): ``-m``, ``m+m``, ``m+=m``, ``m-m``, ``m-=m``, ``m*f``, ``f*m``, ``m*=f``, ``m/f``, ``m/=f``, ``m*m``, ``m*=m``, ``m*v``, ``v*m``, ``m==m``, ``m!=m``.\n\nStatic attributes: ``Zero``, ``Ones``, ``Identity``.",py::init<>())
		.def(py::init<Quaternionr const &>((py::arg("q"))))
		#ifdef _NEW_VISITORS
			.def(MatrixVisitor<Matrix3r>())
		#else
			.def("__init__",py::make_constructor(&Matrix3r_fromElements,py::default_call_policies(),(py::arg("m00"),py::arg("m01"),py::arg("m02"),py::arg("m10"),py::arg("m11"),py::arg("m12"),py::arg("m20"),py::arg("m21"),py::arg("m22"))))
			.def("__init__",py::make_constructor(&Matrix3r_fromRows,py::default_call_policies(),(py::arg("r0"),py::arg("r1"),py::arg("r2"),py::arg("cols")=false)))
			.def("__init__",py::make_constructor(&Matrix3r_fromDiagonal,py::default_call_policies(),(py::arg("diag"))))
			.def(py::init<Matrix3r const &>((py::arg("m"))))
			.def_pickle(Matrix3r_pickle())
			//
			.def("__len__",&::Matrix3r_len).staticmethod("__len__").def("__setitem__",&::Matrix3r_set_item).def("__getitem__",&::Matrix3r_get_item)
			.def("__str__",&::Matrix3r_str).def("__repr__",&::Matrix3r_str)
			.def("__setitem__",&::Matrix3r_set_row).def("__getitem__",&::Matrix3r_get_row)
			.def("determinant",&Matrix3r::determinant)
			.def("trace",&Matrix3r::trace)
			.def("norm",&Matrix3r::norm)
			.def("inverse",&Matrix3r_inverse)
			.def("transpose",&Matrix3r_transpose)
			.def("diagonal",&Matrix3r_diagonal)
			.def("row",&Matrix3r_row)
			.def("col",&Matrix3r_col)
			.def("pruned",&Matrix_pruned<Matrix3r>,py::arg("absTol")=1e-6)
			.def("sum",&Matrix_sum<Matrix3r>)
			.def("maxAbsCoeff",&Matrix_maxAbsCoeff<Matrix3r>)

			//
			.def("__neg__",&Matrix3r__neg__)
			.def("__add__",&Matrix3r__add__Matrix3r).def("__iadd__",&Matrix3r__iadd__Matrix3r)
			.def("__sub__",&Matrix3r__sub__Matrix3r).def("__isub__",&Matrix3r__isub__Matrix3r)
			.def("__mul__",&Matrix3r__mul__Real).def("__rmul__",&Matrix3r__rmul__Real).def("__imul__",&Matrix3r__imul__Real)
			.def("__mul__",&Matrix3r__mul__int).def("__rmul__",&Matrix3r__rmul__int).def("__imul__",&Matrix3r__imul__int)
			.def("__mul__",&Matrix3r__mul__Vector3r).def("__rmul__",&Matrix3r__rmul__Vector3r)
			.def("__mul__",&Matrix3r__mul__Matrix3r).def("__imul__",&Matrix3r__imul__Matrix3r)
			.def("__div__",&Matrix3r__div__Real).def("__idiv__",&Matrix3r__idiv__Real)
			.def("__div__",&Matrix3r__div__int).def("__idiv__",&Matrix3r__idiv__int)
			.def(py::self == py::self)
			.def(py::self != py::self)
			//
			/* extras for matrices */
			.add_static_property("Identity",&Matrix3r_Identity)
			.add_static_property("Zero",&Matrix3r_Zero)
			.add_static_property("Ones",&Matrix3r_Ones)
			.def("Random",&Matrix3r_Random).staticmethod("Random")
			.def("jacobiSVD",&Matrix_jacobiSVD<Matrix3r>,"Compute SVD decomposition of matrix, retuns (U,S,V) such that self=U*S*V.transpose()")
			.def("svd",&Matrix_jacobiSVD<Matrix3r>,"Shortcut for :obj:`jacobiSVD`.")
			.def("computeUnitaryPositive",&Matrix_computeUnitaryPositive<Matrix3r>,"Compute polar decomposition (unitary matrix U and positive semi-definite symmetric matrix P such that self=U*P.")
			.def("polarDecomposition",&Matrix_computeUnitaryPositive<Matrix3r>,"Shortcut for :obj:`computeUnitaryPositive`.")
			.def("selfAdjointEigenDecomposition",&Matrix_selfAdjointEigenDecomposition<Matrix3r>,"Compute eigen (spectral) decomposition of symmetric matrix, returns (eigVecs,eigVals). eigVecs is orthogonal Matrix3 with columns ar normalized eigenvectors, eigVals is Vector3 with corresponding eigenvalues. self=eigVecs*diag(eigVals)*eigVecs.transpose()")
			.def("spectralDecomposition",&Matrix_selfAdjointEigenDecomposition<Matrix3r>,"Shortcut for selfAdjointEigenDecomposition")
		#endif
	;

	py::class_<Matrix6r>("Matrix6","6x6 float matrix. Constructed from 4 3x3 sub-matrices, from 6xVector6 (rows).\n\nSupported operations (``m`` is a Matrix6, ``f`` if a float/int, ``v`` is a Vector6): ``-m``, ``m+m``, ``m+=m``, ``m-m``, ``m-=m``, ``m*f``, ``f*m``, ``m*=f``, ``m/f``, ``m/=f``, ``m*m``, ``m*=m``, ``m*v``, ``v*m``, ``m==m``, ``m!=m``.\n\nStatic attributes: ``Zero``, ``Ones``, ``Identity``.",py::init<>())
		#ifdef _NEW_VISITORS
			.def(MatrixVisitor<Matrix6r>())
		#else
			.def("__init__",py::make_constructor(&Matrix6r_fromBlocks,py::default_call_policies(),(py::arg("ul"),py::arg("ur"),py::arg("ll"),py::arg("lr"))))
			.def("__init__",py::make_constructor(&Matrix6r_fromRows,py::default_call_policies(),(py::arg("l0"),py::arg("l1"),py::arg("l2"),py::arg("l3"),py::arg("l4"),py::arg("l5"),py::arg("cols")=false)))

			/* extras for matrices */
			.def("ul",&Matrix6r_ul)
			.def("ur",&Matrix6r_ur)
			.def("ll",&Matrix6r_ll)
			.def("lr",&Matrix6r_lr)
			.def("__init__",py::make_constructor(&Matrix6r_fromDiagonal,py::default_call_policies(),(py::arg("diag"))))
			.def_pickle(Matrix6r_pickle())
			.def(py::init<Matrix6r const &>((py::arg("m"))))
			//
			.def("__str__",&::Matrix6r_str).def("__repr__",&::Matrix6r_str)
			.def("__len__",&::Matrix6r_len).staticmethod("__len__").def("__setitem__",&::Matrix6r_set_item).def("__getitem__",&::Matrix6r_get_item)
			.def("__setitem__",&::Matrix6r_set_row).def("__getitem__",&::Matrix6r_get_row)
			.def("determinant",&Matrix6r::determinant)
			.def("trace",&Matrix6r::trace)
			.def("inverse",&Matrix6r_inverse)
			.def("transpose",&Matrix6r_transpose)
			//.def("polarDecomposition",&Matrix6r_polarDecomposition)
			.def("diagonal",&Matrix6r_diagonal)
			.def("sum",&Matrix_sum<Matrix6r>)
			.def("row",&Matrix6r_row)
			.def("col",&Matrix6r_col)
			//
			.def("jacobiSVD",&Matrix_jacobiSVD<Matrix6r>,"Compute SVD decomposition of matrix, returns (U,S,V) such that self=U*S*V.transpose()")
			.def("svd",&Matrix_jacobiSVD<Matrix6r>,"Shortcut for :obj:`jacobiSVD`.")
			.def("computeUnitaryPositive",&Matrix_computeUnitaryPositive<Matrix6r>,"Compute polar decomposition (unitary matrix U and positive semi-definite symmetric matrix P such that self=U*P")
			.def("polarDecomposition",&Matrix_computeUnitaryPositive<Matrix6r>,"Shortcut for computeUnitaryPositive")
			.def("selfAdjointEigenDecomposition",&Matrix_selfAdjointEigenDecomposition<Matrix6r>,"Compute eigen (spectral) decomposition of symmetric matrix, returns (eigVecs,eigVals). eigVecs is orthogonal Matrix6 with columns ar normalized eigenvectors, eigVals is Vector6 with corresponding eigenvalues. self=eigVecs*diag(eigVals)*eigVecs.transpose()")
			.def("spectralDecomposition",&Matrix_selfAdjointEigenDecomposition<Matrix6r>,"Shortcut for :obj:`selfAdjointEigenDecomposition`.")

			//
			.def("__neg__",&Matrix6r__neg__)
			.def("__add__",&Matrix6r__add__Matrix6r).def("__iadd__",&Matrix6r__iadd__Matrix6r)
			.def("__sub__",&Matrix6r__sub__Matrix6r).def("__isub__",&Matrix6r__isub__Matrix6r)
			.def("__mul__",&Matrix6r__mul__Real).def("__rmul__",&Matrix6r__rmul__Real).def("__imul__",&Matrix6r__imul__Real)
			.def("__mul__",&Matrix6r__mul__Matrix6r).def("__imul__",&Matrix6r__imul__Matrix6r)
			.def("__mul__",&Matrix6r__mul__int).def("__rmul__",&Matrix6r__rmul__int).def("__imul__",&Matrix6r__imul__int)
			.def("__mul__",&Matrix6r__mul__Vector6r).def("__rmul__",&Matrix6r__rmul__Vector6r)
			.def("__div__",&Matrix6r__div__Real).def("__idiv__",&Matrix6r__idiv__Real)
			.def("__div__",&Matrix6r__div__int).def("__idiv__",&Matrix6r__idiv__int)
			.def(py::self == py::self)
			.def(py::self != py::self)
			//
			.add_static_property("Identity",&Matrix6r_Identity)
			.add_static_property("Zero",&Matrix6r_Zero)
			.add_static_property("Ones",&Matrix6r_Ones)
			.def("Random",&Matrix6r_Random).staticmethod("Random")
		#endif
	;

	py::class_<VectorXr>("VectorX","Dynamic-sized float vector.\n\nSupported operations (``f`` if a float/int, ``v`` is a VectorX): ``-v``, ``v+v``, ``v+=v``, ``v-v``, ``v-=v``, ``v*f``, ``f*v``, ``v*=f``, ``v/f``, ``v/=f``, ``v==v``, ``v!=v``.\n\nImplicit conversion from sequence (list, tuple, ...) of X floats.",py::init<>())
		#ifdef _NEW_VISITORS
			.def(VectorVisitor<VectorXr>())
		#else
			// properties
			// methods
			.def("__init__",py::make_constructor(&VectorXr_fromList,py::default_call_policies(),(py::arg("vv"))))
			.def(py::init<VectorXr>((py::arg("other"))))
			.add_static_property("Ones",&VectorXr_Ones).add_static_property("Zero",&VectorXr_Zero)
			.def_pickle(VectorXr_pickle())
			.def("resize",&VectorXr_resize)
			.def("norm",&VectorXr::norm).def("squaredNorm",&VectorXr::squaredNorm).def("normalize",&VectorXr::normalize).def("normalized",&VectorXr::normalized)
			.def("asDiagonal",&Vector_asDiagonal<VectorXr>)
			.def("size",&VectorXr::size)
			.def("sum",&Matrix_sum<VectorXr>)
			.def("Random",&VectorXr_Random).staticmethod("Random")
			// operators
			.def("__neg__",&VectorXr__neg__) // -v
			.def("__add__",&VectorXr__add__VectorXr).def("__iadd__",&VectorXr__iadd__VectorXr) // +, +=
			.def("__sub__",&VectorXr__sub__VectorXr).def("__isub__",&VectorXr__isub__VectorXr) // -, -=
			.def("__mul__",&VectorXr__mul__Real).def("__rmul__",&VectorXr__rmul__Real) // f*v, v*f
			.def("__div__",&VectorXr__div__Real).def("__idiv__",&VectorXr__idiv__Real) // v/f, v/=f
			.def("__mul__",&VectorXr__mul__int).def("__rmul__",&VectorXr__rmul__int) // f*v, v*f
			.def("__div__",&VectorXr__div__int).def("__idiv__",&VectorXr__idiv__int) // v/f, v/=f
			.def(py::self != py::self).def(py::self == py::self)
			// specials
			.def("__abs__",&VectorXr::norm)
			.def("__len__",&::VectorXr_len)
			.def("__setitem__",&::VectorXr_set_item).def("__getitem__",&::VectorXr_get_item)
			.def("__str__",&::VectorXr_str).def("__repr__",&::VectorXr_str)
		#endif
	;


	py::class_<MatrixXr>("MatrixX","XxX (dynamic-sized) float matrix. Constructed from list of rows (as VectorX).\n\nSupported operations (``m`` is a MatrixX, ``f`` if a float/int, ``v`` is a VectorX): ``-m``, ``m+m``, ``m+=m``, ``m-m``, ``m-=m``, ``m*f``, ``f*m``, ``m*=f``, ``m/f``, ``m/=f``, ``m*m``, ``m*=m``, ``m*v``, ``v*m``, ``m==m``, ``m!=m``.",py::init<>())
		/* extras for matrices */
		#ifdef _NEW_VISITORS
			.def(MatrixVisitor<MatrixXr>())
		#else
			.def("__init__",py::make_constructor(&MatrixXr_fromRows,py::default_call_policies(),(
				py::arg("r0")=VectorXr(),
				py::arg("r1")=VectorXr(),
				py::arg("r2")=VectorXr(),
				py::arg("r3")=VectorXr(),
				py::arg("r4")=VectorXr(),
				py::arg("r5")=VectorXr(),
				py::arg("r6")=VectorXr(),
				py::arg("r7")=VectorXr(),
				py::arg("r8")=VectorXr(),
				py::arg("r9")=VectorXr(),
				py::arg("cols")=false)
			))
			.def("__init__",py::make_constructor(&MatrixXr_fromDiagonal,py::default_call_policies(),(py::arg("diag"))))
			.def("__init__",py::make_constructor(&MatrixXr_fromRowSeq,py::default_call_policies(),(py::arg("rows"),py::arg("cols")=false)))
			.def(py::init<MatrixXr const &>((py::arg("m"))))
			.def_pickle(MatrixXr_pickle())
			.def("rows",&MatrixXr::rows)
			.def("cols",&MatrixXr::cols)
			.def("resize",&MatrixXr_resize)
			//
	 		.def("__len__",&::MatrixXr_len).def("__setitem__",&::MatrixXr_set_item).def("__getitem__",&::MatrixXr_get_item)
			.def("__str__",&::MatrixXr_str).def("__repr__",&::MatrixXr_str)
			.def("__setitem__",&::MatrixXr_set_row).def("__getitem__",&::MatrixXr_get_row)
			.def("determinant",&MatrixXr::determinant)
			.def("trace",&MatrixXr::trace)
			.def("inverse",&MatrixXr_inverse)
			.def("transpose",&MatrixXr_transpose)
			#if 0
				.def("pinv",&MatrixXr_pseudoInverse_wrapper)
			#endif
			.def("diagonal",&MatrixXr_diagonal)
			.def("row",&MatrixXr_row)
			.def("col",&MatrixXr_col)
			.def("pruned",&Matrix_pruned<MatrixXr>,py::arg("absTol")=1e-6)
			.def("sum",&Matrix_sum<MatrixXr>)
			.def("maxAbsCoeff",&Matrix_maxAbsCoeff<MatrixXr>)
			//
			.def("jacobiSVD",&Matrix_jacobiSVD<MatrixXr>,"Compute SVD decomposition of matrix, returns (U,S,V) such that self=U*S*V.transpose()")
			.def("svd",&Matrix_jacobiSVD<MatrixXr>,"Shortcut for :obj:`jacobiSVD`.")
			.def("computeUnitaryPositive",&Matrix_computeUnitaryPositive<MatrixXr>,"Compute polar decomposition (unitary matrix U and positive semi-definite symmetric matrix P such that self=U*P")
			.def("polarDecomposition",&Matrix_computeUnitaryPositive<MatrixXr>,"Shortcut for computeUnitaryPositive")
			.def("selfAdjointEigenDecomposition",&Matrix_selfAdjointEigenDecomposition<MatrixXr>,"Compute eigen (spectral) decomposition of symmetric matrix, returns (eigVecs,eigVals). eigVecs is orthogonal MatrixX with columns ar normalized eigenvectors, eigVals is VectorX with corresponding eigenvalues. self=eigVecs*diag(eigVals)*eigVecs.transpose()")
			.def("spectralDecomposition",&Matrix_selfAdjointEigenDecomposition<MatrixXr>,"Shortcut for selfAdjointEigenDecomposition")

			//
			.def("__neg__",&MatrixXr__neg__)
			.def("__add__",&MatrixXr__add__MatrixXr).def("__iadd__",&MatrixXr__iadd__MatrixXr)
			.def("__sub__",&MatrixXr__sub__MatrixXr).def("__isub__",&MatrixXr__isub__MatrixXr)
			.def("__mul__",&MatrixXr__mul__Real).def("__rmul__",&MatrixXr__rmul__Real).def("__imul__",&MatrixXr__imul__Real)
			.def("__mul__",&MatrixXr__mul__int).def("__rmul__",&MatrixXr__rmul__int).def("__imul__",&MatrixXr__imul__int)
			.def("__mul__",&MatrixXr__mul__MatrixXr).def("__imul__",&MatrixXr__imul__MatrixXr)
			.def("__div__",&MatrixXr__div__Real).def("__idiv__",&MatrixXr__idiv__Real)
			.def("__div__",&MatrixXr__div__int).def("__idiv__",&MatrixXr__idiv__int)
			.def("__mul__",&MatrixXr__mul__VectorXr).def("__rmul__",&MatrixXr__rmul__VectorXr)
			.def(py::self == py::self)
			.def(py::self != py::self)
			//
			.add_static_property("Identity",&MatrixXr_Identity)
			.add_static_property("Zero",&MatrixXr_Zero)
			.add_static_property("Ones",&MatrixXr_Ones)
			.def("Random",&MatrixXr_Random).staticmethod("Random")
		#endif
	;

	py::class_<Vector6r>("Vector6","6-dimensional float vector.\n\nSupported operations (``f`` if a float/int, ``v`` is a Vector6): ``-v``, ``v+v``, ``v+=v``, ``v-v``, ``v-=v``, ``v*f``, ``f*v``, ``v*=f``, ``v/f``, ``v/=f``, ``v==v``, ``v!=v``.\n\nImplicit conversion from sequence (list, tuple, ...) of 6 floats.\n\nStatic attributes: ``Zero``, ``Ones``.",py::init<>())
		#ifdef _NEW_VISITORS
			.def(VectorVisitor<Vector6r>())
		#else
			.def("__init__",py::make_constructor(&Vector6r_fromElements,py::default_call_policies(),(py::arg("v0"),py::arg("v1"),py::arg("v2"),py::arg("v3"),py::arg("v4"),py::arg("v5"))))
			.def("__init__",py::make_constructor(&Vector6r_fromHeadTail,py::default_call_policies(),(py::arg("head"),py::arg("tail"))))
			.def("head",&Vector6r_head).def("tail",&Vector6r_tail)
			.def(py::init<Vector6r>((py::arg("other"))))
			.def_pickle(Vector6r_pickle())
			// properties
			.add_static_property("Ones",&Vector6r_Ones).add_static_property("Zero",&Vector6r_Zero)
			//.add_static_property("UnitX",&Vector6r_UnitX).add_static_property("UnitY",&Vector6r_UnitY).add_static_property("UnitZ",&Vector6r_UnitZ)
			// methods
			//.def("dot",&Vector6r_dot).def("cross",&Vector6r_cross)
			.def("Random",&Vector6r_Random).staticmethod("Random")
			.def("norm",&Vector6r::norm).def("squaredNorm",&Vector6r::squaredNorm).def("normalize",&Vector6r::normalize).def("normalized",&Vector6r::normalized)
			.def("asDiagonal",&Vector_asDiagonal<Vector6r>)
			.def("sum",&Matrix_sum<Vector6r>)
			// operators
			.def("__neg__",&Vector6r__neg__) // -v
			.def("__add__",&Vector6r__add__Vector6r).def("__iadd__",&Vector6r__iadd__Vector6r) // +, +=
			.def("__sub__",&Vector6r__sub__Vector6r).def("__isub__",&Vector6r__isub__Vector6r) // -, -=
			.def("__mul__",&Vector6r__mul__Real).def("__rmul__",&Vector6r__rmul__Real) // f*v, v*f
			.def("__div__",&Vector6r__div__Real).def("__idiv__",&Vector6r__idiv__Real) // v/f, v/=f
			.def("__mul__",&Vector6r__mul__int).def("__rmul__",&Vector6r__rmul__int) // f*v, v*f
			.def("__div__",&Vector6r__div__int).def("__idiv__",&Vector6r__idiv__int) // v/f, v/=f
			.def(py::self != py::self).def(py::self == py::self)
			// specials
			.def("__abs__",&Vector6r::norm)
			.def("__len__",&::Vector6r_len).staticmethod("__len__")
			.def("__setitem__",&::Vector6r_set_item).def("__getitem__",&::Vector6r_get_item)
			.def("__str__",&::Vector6r_str).def("__repr__",&::Vector6r_str)
		#endif
		#if 0
			// specials
			.def("toSymmTensor",&Vector6r_toSymmTensor,(py::args("strain")=false),"Convert Vector6 in the Voigt notation to the corresponding 2nd order symmetric tensor (as Matrix3); if *strain* is ``True``, multiply non-diagonal components by .5")
		#endif
	;

	py::class_<Vector6i>("Vector6i","6-dimensional float vector.\n\nSupported operations (``f`` if a float/int, ``v`` is a Vector6): ``-v``, ``v+v``, ``v+=v``, ``v-v``, ``v-=v``, ``v*f``, ``f*v``, ``v*=f``, ``v/f``, ``v/=f``, ``v==v``, ``v!=v``.\n\nImplicit conversion from sequence (list, tuple, ...) of 6 floats.\n\nStatic attributes: ``Zero``, ``Ones``.",py::init<>())
		#ifdef _NEW_VISITORS
			.def(VectorVisitor<Vector6i>())
		#else
			.def("__init__",py::make_constructor(&Vector6i_fromElements,py::default_call_policies(),(py::arg("v0"),py::arg("v1"),py::arg("v2"),py::arg("v3"),py::arg("v4"),py::arg("v5"))))
			.def("head",&Vector6i_head).def("tail",&Vector6i_tail)
			.def(py::init<Vector6i>((py::arg("other"))))
			.def_pickle(Vector6i_pickle())
			// properties
			.add_static_property("Ones",&Vector6i_Ones).add_static_property("Zero",&Vector6i_Zero)
			//.add_static_property("UnitX",&Vector6i_UnitX).add_static_property("UnitY",&Vector6i_UnitY).add_static_property("UnitZ",&Vector6i_UnitZ)
			// methods
			//.def("dot",&Vector6i_dot).def("cross",&Vector6i_cross)
			//.def("norm",&Vector6i::norm).def("squaredNorm",&Vector6i::squaredNorm).def("normalize",&Vector6i::normalize).def("normalized",&Vector6i::normalized)
			.def("sum",&Matrix_sum<Vector6i>)
			// operators
			.def("__neg__",&Vector6i__neg__) // -v
			.def("__add__",&Vector6i__add__Vector6i).def("__iadd__",&Vector6i__iadd__Vector6i) // +, +=
			.def("__sub__",&Vector6i__sub__Vector6i).def("__isub__",&Vector6i__isub__Vector6i) // -, -=
			.def("__mul__",&Vector6i__mul__int).def("__rmul__",&Vector6i__rmul__int) // f*v, v*f
			.def("__div__",&Vector6i__div__int).def("__idiv__",&Vector6i__idiv__int) // v/f, v/=f
			.def(py::self != py::self).def(py::self == py::self)
			// specials
			// .def("__abs__",&Vector6i::norm)
			.def("__len__",&::Vector6i_len).staticmethod("__len__")
			.def("__setitem__",&::Vector6i_set_item).def("__getitem__",&::Vector6i_get_item)
			.def("__str__",&::Vector6i_str).def("__repr__",&::Vector6i_str)
		#endif
	;

	py::class_<Vector3r>("Vector3","3-dimensional float vector.\n\nSupported operations (``f`` if a float/int, ``v`` is a Vector3): ``-v``, ``v+v``, ``v+=v``, ``v-v``, ``v-=v``, ``v*f``, ``f*v``, ``v*=f``, ``v/f``, ``v/=f``, ``v==v``, ``v!=v``, plus operations with ``Matrix3`` and ``Quaternion``.\n\nImplicit conversion from sequence (list, tuple, ...) of 3 floats.\n\nStatic attributes: ``Zero``, ``Ones``, ``UnitX``, ``UnitY``, ``UnitZ``.",py::init<>())
		#ifdef _NEW_VISITORS
			.def(VectorVisitor<Vector3r>())
		#else
			// swizzles
			.def("xy",&Vector3r_xy).def("yx",&Vector3r_yx).def("xz",&Vector3r_xz).def("zx",&Vector3r_zx).def("yz",&Vector3r_yz).def("zy",&Vector3r_zy)
			.def(py::init<Real,Real,Real>((py::arg("x"),py::arg("y"),py::arg("z"))))
			.def_pickle(Vector3r_pickle())
			.def(py::init<Vector3r>((py::arg("other"))))
			// properties
			.add_static_property("Ones",&Vector3r_Ones).add_static_property("Zero",&Vector3r_Zero)
			.add_static_property("UnitX",&Vector3r_UnitX).add_static_property("UnitY",&Vector3r_UnitY).add_static_property("UnitZ",&Vector3r_UnitZ)
			// methods
			.def("Random",&Vector3r_Random).staticmethod("Random")
			.def("dot",&Vector3r_dot).def("cross",&Vector3r_cross)
			.def("norm",&Vector3r::norm).def("squaredNorm",&Vector3r::squaredNorm).def("normalize",&Vector3r::normalize).def("normalized",&Vector3r::normalized)
			.def("asDiagonal",&Vector_asDiagonal<Vector3r>)
			.def("Unit",&Vector_Unit<Vector3r>).staticmethod("Unit")
			.def("pruned",&Matrix_pruned<Vector3r>,py::arg("absTol")=1e-6)
			.def("maxAbsCoeff",&Matrix_maxAbsCoeff<Vector3r>)
			.def("sum",&Matrix_sum<Vector3r>)
			.def("outer",&Vector_outer<Vector3r,Vector3r>)
			// operators
			.def("__neg__",&Vector3r__neg__) // -v
			.def("__add__",&Vector3r__add__Vector3r).def("__iadd__",&Vector3r__iadd__Vector3r) // +, +=
			.def("__sub__",&Vector3r__sub__Vector3r).def("__isub__",&Vector3r__isub__Vector3r) // -, -=
			.def("__mul__",&Vector3r__mul__Real).def("__rmul__",&Vector3r__rmul__Real) // f*v, v*f
			.def("__div__",&Vector3r__div__Real).def("__idiv__",&Vector3r__idiv__Real) // v/f, v/=f
			.def("__mul__",&Vector3r__mul__int).def("__rmul__",&Vector3r__rmul__int) // f*v, v*f
			.def("__div__",&Vector3r__div__int).def("__idiv__",&Vector3r__idiv__int) // v/f, v/=f
			.def(py::self != py::self).def(py::self == py::self)
			// specials
			.def("__abs__",&Vector3r::norm)
			.def("__len__",&::Vector3r_len).staticmethod("__len__")
			.def("__setitem__",&::Vector3r_set_item).def("__getitem__",&::Vector3r_get_item)
			.def("__str__",&::Vector3r_str).def("__repr__",&::Vector3r_str)
		#endif
	;	
	py::class_<Vector3i>("Vector3i","3-dimensional integer vector.\n\nSupported operations (``i`` if an int, ``v`` is a Vector3i): ``-v``, ``v+v``, ``v+=v``, ``v-v``, ``v-=v``, ``v*i``, ``i*v``, ``v*=i``, ``v==v``, ``v!=v``.\n\nImplicit conversion from sequence  (list, tuple, ...) of 3 integers.\n\nStatic attributes: ``Zero``, ``Ones``, ``UnitX``, ``UnitY``, ``UnitZ``.",py::init<>())
		#ifdef _NEW_VISITORS
			.def(VectorVisitor<Vector3i>())
		#else
			.def(py::init<int,int,int>((py::arg("x"),py::arg("y"),py::arg("z"))))
			.def("cross",&Vector3i_cross)
			.def(py::init<Vector3i>((py::arg("other"))))
			.def_pickle(Vector3i_pickle())
			.def("dot",&Vector3i_dot)
			.add_static_property("UnitX",&Vector3i_UnitX).add_static_property("UnitY",&Vector3i_UnitY).add_static_property("UnitZ",&Vector3i_UnitZ)
			// properties
			.add_static_property("Ones",&Vector3i_Ones).add_static_property("Zero",&Vector3i_Zero)
			.def("Unit",&Vector_Unit<Vector3i>).staticmethod("Unit")
			// methods
			.def("sum",&Matrix_sum<Vector3i>)
			//.def("norm",&Vector3i::norm).def("squaredNorm",&Vector3i::squaredNorm)
			// operators
			.def("__neg__",&Vector3i__neg__) // -v
			.def("__add__",&Vector3i__add__Vector3i).def("__iadd__",&Vector3i__iadd__Vector3i) // +, +=
			.def("__sub__",&Vector3i__sub__Vector3i).def("__isub__",&Vector3i__isub__Vector3i) // -, -=
			.def("__mul__",&Vector3i__mul__int).def("__rmul__",&Vector3i__rmul__int) // f*v, v*f
			.def(py::self != py::self).def(py::self == py::self)
			// specials
			//.def("__abs__",&Vector3i::norm)
			.def("__len__",&::Vector3i_len).staticmethod("__len__")
			.def("__setitem__",&::Vector3i_set_item).def("__getitem__",&::Vector3i_get_item)
			.def("__str__",&::Vector3i_str).def("__repr__",&::Vector3i_str)
		#endif
	;	
	py::class_<Vector2r>("Vector2","3-dimensional float vector.\n\nSupported operations (``f`` if a float/int, ``v`` is a Vector3): ``-v``, ``v+v``, ``v+=v``, ``v-v``, ``v-=v``, ``v*f``, ``f*v``, ``v*=f``, ``v/f``, ``v/=f``, ``v==v``, ``v!=v``.\n\nImplicit conversion from sequence (list, tuple, ...) of 2 floats.\n\nStatic attributes: ``Zero``, ``Ones``, ``UnitX``, ``UnitY``.",py::init<>())
		#ifdef _NEW_VISITORS
			.def(VectorVisitor<Vector2r>())
		#else
			.def(py::init<Real,Real>((py::arg("x"),py::arg("y"))))
			.def("Unit",&Vector_Unit<Vector2r>).staticmethod("Unit")
			.add_static_property("UnitX",&Vector2r_UnitX).add_static_property("UnitY",&Vector2r_UnitY)
			.def("Random",&Vector2r_Random).staticmethod("Random")
			.def("dot",&Vector2r_dot)
			.def(py::init<Vector2r>((py::arg("other"))))
			.def_pickle(Vector2r_pickle())
			// properties
			.add_static_property("Ones",&Vector2r_Ones).add_static_property("Zero",&Vector2r_Zero)
			// methods
			.def("sum",&Matrix_sum<Vector2r>)
			// operators
			.def("__neg__",&Vector2r__neg__) // -v
			.def("__add__",&Vector2r__add__Vector2r).def("__iadd__",&Vector2r__iadd__Vector2r) // +, +=
			.def("__sub__",&Vector2r__sub__Vector2r).def("__isub__",&Vector2r__isub__Vector2r) // -, -=
			.def("__mul__",&Vector2r__mul__Real).def("__rmul__",&Vector2r__rmul__Real) // f*v, v*f
			.def("__div__",&Vector2r__div__Real).def("__idiv__",&Vector2r__idiv__Real) // v/f, v/=f
			.def("__mul__",&Vector2r__mul__int).def("__rmul__",&Vector2r__rmul__int) // f*v, v*f
			.def("__div__",&Vector2r__div__int).def("__idiv__",&Vector2r__idiv__int) // v/f, v/=f
			.def(py::self != py::self).def(py::self == py::self)
			// specials
			.def("__len__",&::Vector2r_len).staticmethod("__len__")
			.def("__setitem__",&::Vector2r_set_item).def("__getitem__",&::Vector2r_get_item)
			.def("__str__",&::Vector2r_str).def("__repr__",&::Vector2r_str)
			.def("norm",&Vector2r::norm).def("squaredNorm",&Vector2r::squaredNorm).def("normalize",&Vector2r::normalize).def("normalized",&Vector2r::normalized)
			.def("__abs__",&Vector2r::norm)
		#endif
	;	
	py::class_<Vector2i>("Vector2i","2-dimensional integer vector.\n\nSupported operations (``i`` if an int, ``v`` is a Vector2i): ``-v``, ``v+v``, ``v+=v``, ``v-v``, ``v-=v``, ``v*i``, ``i*v``, ``v*=i``, ``v==v``, ``v!=v``.\n\nImplicit conversion from sequence (list, tuple, ...) of 2 integers.\n\nStatic attributes: ``Zero``, ``Ones``, ``UnitX``, ``UnitY``.",py::init<>())
		#ifdef _NEW_VISITORS
			.def(VectorVisitor<Vector2i>())
		#else
			.def(py::init<int,int>((py::arg("x"),py::arg("y"))))
			.def(py::init<Vector2i>((py::arg("other"))))
			.def_pickle(Vector2i_pickle())
			// properties
			.add_static_property("Ones",&Vector2i_Ones).add_static_property("Zero",&Vector2i_Zero)
			.add_static_property("UnitX",&Vector2i_UnitX).add_static_property("UnitY",&Vector2i_UnitY)
			// methods
			.def("dot",&Vector2i_dot)
			.def("sum",&Matrix_sum<Vector2i>)
			//.def("norm",&Vector2i::norm).def("squaredNorm",&Vector2i::squaredNorm).def("normalize",&Vector2i::normalize)
			.def("Unit",&Vector_Unit<Vector2i>).staticmethod("Unit")
			// operators
			.def("__neg__",&Vector2i__neg__) // -v
			.def("__add__",&Vector2i__add__Vector2i).def("__iadd__",&Vector2i__iadd__Vector2i) // +, +=
			.def("__sub__",&Vector2i__sub__Vector2i).def("__isub__",&Vector2i__isub__Vector2i) // -, -=
			.def("__mul__",&Vector2i__mul__int).def("__rmul__",&Vector2i__rmul__int) // f*v, v*f
			.def(py::self != py::self).def(py::self == py::self)
			// specials
			//.def("__abs__",&Vector2i::norm)
			.def("__len__",&::Vector2i_len).staticmethod("__len__")
			.def("__setitem__",&::Vector2i_set_item).def("__getitem__",&::Vector2i_get_item)
			.def("__str__",&::Vector2i_str).def("__repr__",&::Vector2i_str)
		#endif
	;	
	

	py::class_<AlignedBox3r>("AlignedBox3","Axis-aligned box object, defined by its minimum and maximum corners",py::init<>())
		#ifdef _NEW_VISITORS
			.def(AabbVisitor<AlignedBox3r>())
		#else
			.def(py::init<AlignedBox3r>((py::arg("other"))))
			.def(py::init<Vector3r,Vector3r>((py::arg("min"),py::arg("max"))))
			.def_pickle(AlignedBox3r_pickle())
			.def("volume",&AlignedBox3r::volume)
			.def("empty",&::AlignedBox3r::isEmpty)
			.def("center",&::AlignedBox3r_center)
			.def("sizes",&::AlignedBox3r_size)
			.def("contains",&::AlignedBox3r_containsPt)
			.def("contains",&::AlignedBox3r_containsBox)
			// for the "in" operator
			.def("__contains__",&::AlignedBox3r_containsPt) 
			.def("__contains__",&::AlignedBox3r_containsBox)
			.def("extend",&::AlignedBox3r_extendPt)
			.def("extend",&::AlignedBox3r_extendBox)
			.def("clamp",&::AlignedBox3r_clamp)
			// return new objects
			.def("intersection",&AlignedBox3r::intersection)
			.def("merged",&AlignedBox3r::merged)
			// those return internal references, which is what we want
			.add_property("min",&::AlignedBox3r_min) 
			.add_property("max",&::AlignedBox3r_max)
			.def("__len__",&::AlignedBox3r_len).staticmethod("__len__")
			.def("__setitem__",&::AlignedBox3r_set_item).def("__getitem__",&::AlignedBox3r_get_item)
			.def("__setitem__",&::AlignedBox3r_set_minmax).def("__getitem__",&::AlignedBox3r_get_minmax)
			.def("__str__",&::AlignedBox3r_str).def("__repr__",&::AlignedBox3r_str)
		#endif
	;
	py::class_<AlignedBox2r>("AlignedBox2","Axis-aligned box object in 2d, defined by its minimum and maximum corners",py::init<>())
		#ifdef _NEW_VISITORS
			.def(AabbVisitor<AlignedBox2r>())
		#else
			.def(py::init<AlignedBox2r>((py::arg("other"))))
			.def(py::init<Vector2r,Vector2r>((py::arg("min"),py::arg("max"))))
			.def_pickle(AlignedBox2r_pickle())
			.def("volume",&AlignedBox2r::volume)
			.def("empty",&::AlignedBox2r::isEmpty)
			.def("center",&::AlignedBox2r_center)
			.def("sizes",&::AlignedBox2r_size)
			.def("area",&AlignedBox2r::volume)
			.def("contains",&::AlignedBox2r_containsPt)
			.def("contains",&::AlignedBox2r_containsBox)
			.def("extend",&::AlignedBox2r_extendPt)
			.def("extend",&::AlignedBox2r_extendBox)
			.def("clamp",&::AlignedBox2r_clamp)
			// for the "in" operator
			.def("__contains__",&::AlignedBox2r_containsPt) 
			.def("__contains__",&::AlignedBox2r_containsBox)
			.def("intersection",&AlignedBox2r::intersection)
			.def("merged",&AlignedBox2r::merged)
			// those return internal references, which is what we want
			.add_property("min",&::AlignedBox2r_min) 
			.add_property("max",&::AlignedBox2r_max)
			.def("__len__",&::AlignedBox2r_len).staticmethod("__len__")
			.def("__setitem__",&::AlignedBox2r_set_item).def("__getitem__",&::AlignedBox2r_get_item)
			.def("__setitem__",&::AlignedBox2r_set_minmax).def("__getitem__",&::AlignedBox2r_get_minmax)
			.def("__str__",&::AlignedBox2r_str).def("__repr__",&::AlignedBox2r_str)
		#endif
	;
};








