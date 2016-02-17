#pragma once
#include"common.hpp"
// classes dealing with exposing actual types, with many switches inside depending on template arguments

// methods common for vectors and matrices
template<typename MatrixBaseT>
class MatrixBaseVisitor: public py::def_visitor<MatrixBaseVisitor<MatrixBaseT> >{
	typedef typename MatrixBaseT::Scalar Scalar; // could be complex number
	typedef typename MatrixBaseT::RealScalar RealScalar; // this is the "real" (math) scalar
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
		.def("isApprox",&MatrixBaseVisitor::isApprox,(py::arg("other"),py::arg("prec")=Eigen::NumTraits<Scalar>::dummy_precision()),"Approximate comparison with precision *prec*.")
		.def("rows",&MatrixBaseT::rows,"Number of rows.")
		.def("cols",&MatrixBaseT::cols,"Number of columns.")
		;
		visit_if_float<Scalar,PyClass>(cl);
		visit_fixed_or_dynamic<MatrixBaseT,PyClass>(cl);

		// reductions
		cl
		.def("sum",&MatrixBaseT::sum,"Sum of all elements.")
		.def("prod",&MatrixBaseT::prod,"Product of all elements.")
		.def("mean",&MatrixBaseT::mean,"Mean value over all elements.")
		.def("maxAbsCoeff",&MatrixBaseVisitor::maxAbsCoeff,"Maximum absolute value over all elements.")
		;
		// reductions only meaningful for non-complex matrices (e.g. maxCoeff, minCoeff)
		visit_reductions_noncomplex<Scalar,PyClass>(cl);
		
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
		// __trudiv__ is for py3k
		cl
		.def("__mul__",&MatrixBaseVisitor::__mul__scalar<Scalar>)
		.def("__rmul__",&MatrixBaseVisitor::__rmul__scalar<Scalar>)
		.def("__imul__",&MatrixBaseVisitor::__imul__scalar<Scalar>)
		.def("__div__",&MatrixBaseVisitor::__div__scalar<long>)
		.def("__truediv__",&MatrixBaseVisitor::__div__scalar<long>)
		.def("__idiv__",&MatrixBaseVisitor::__idiv__scalar<long>)
		.def("__itruediv__",&MatrixBaseVisitor::__div__scalar<long>)
		.def("__div__",&MatrixBaseVisitor::__div__scalar<Scalar>)
		.def("__truediv__",&MatrixBaseVisitor::__div__scalar<Scalar>)
		.def("__idiv__",&MatrixBaseVisitor::__idiv__scalar<Scalar>)
		.def("__itruediv__",&MatrixBaseVisitor::__idiv__scalar<Scalar>)
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
	static RealScalar maxAbsCoeff(const MatrixBaseT& m){ return m.array().abs().maxCoeff(); }
	static MatrixBaseT Ones(){ return MatrixBaseT::Ones(); }
	static MatrixBaseT Zero(){ return MatrixBaseT::Zero(); }
	static MatrixBaseT Random(){ return MatrixBaseT::Random(); }
	static MatrixBaseT Identity(){ return MatrixBaseT::Identity(); }

	static bool __eq__(const MatrixBaseT& a, const MatrixBaseT& b){
		if(a.rows()!=b.rows() || a.cols()!=b.cols()) return false;
		return a.cwiseEqual(b).all();
	}
	static bool isApprox(const MatrixBaseT& a, const MatrixBaseT& b, const RealScalar& prec){ return a.isApprox(b,prec); }
	static bool __ne__(const MatrixBaseT& a, const MatrixBaseT& b){ return !__eq__(a,b); }
	static MatrixBaseT __neg__(const MatrixBaseT& a){ return -a; };
	static MatrixBaseT __add__(const MatrixBaseT& a, const MatrixBaseT& b){ return a+b; }
	static MatrixBaseT __sub__(const MatrixBaseT& a, const MatrixBaseT& b){ return a-b; }
	static MatrixBaseT __iadd__(MatrixBaseT& a, const MatrixBaseT& b){ a+=b; return a; };
	static MatrixBaseT __isub__(MatrixBaseT& a, const MatrixBaseT& b){ a-=b; return a; };

	template<typename Scalar2> static MatrixBaseT __mul__scalar(const MatrixBaseT& a, const Scalar2& scalar){ return a*scalar; }
	template<typename Scalar2> static MatrixBaseT __imul__scalar(MatrixBaseT& a, const Scalar2& scalar){ a*=scalar; return a; }
	template<typename Scalar2> static MatrixBaseT __rmul__scalar(const MatrixBaseT& a, const Scalar2& scalar){ return a*scalar; }
	template<typename Scalar2> static MatrixBaseT __div__scalar(const MatrixBaseT& a, const Scalar2& scalar){ return a/scalar; }
	template<typename Scalar2> static MatrixBaseT __idiv__scalar(MatrixBaseT& a, const Scalar2& scalar){ a/=scalar; return a; }

	template<typename Scalar, class PyClass> static	void visit_reductions_noncomplex(PyClass& cl, typename boost::enable_if_c<Eigen::NumTraits<Scalar>::IsComplex >::type* dummy = 0){ /* do nothing*/ }
	template<typename Scalar, class PyClass> static	void visit_reductions_noncomplex(PyClass& cl, typename boost::disable_if_c<Eigen::NumTraits<Scalar>::IsComplex >::type* dummy = 0){
		// must be wrapped since there are overloads:
		//   maxCoeff(), maxCoeff(IndexType*), maxCoeff(IndexType*,IndexType*)
		cl
		.def("maxCoeff",&MatrixBaseVisitor::maxCoeff0,"Maximum value over all elements.")
		.def("minCoeff",&MatrixBaseVisitor::minCoeff0,"Minimum value over all elements.")
		;
	}
	static Scalar maxCoeff0(const MatrixBaseT& m){ return m.array().maxCoeff(); }
	static Scalar minCoeff0(const MatrixBaseT& m){ return m.array().minCoeff(); }

	// we want to keep -0 (rather than replacing it by 0), but that does not work for complex numbers
	// hence two versions
	template<typename Scalar> static bool prune_element(const Scalar& num, RealScalar absTol, typename boost::disable_if_c<Eigen::NumTraits<Scalar>::IsComplex >::type* dummy=0){using namespace std; return abs(num)<=absTol || num!=-0; }
	template<typename Scalar> static bool prune_element(const Scalar& num, RealScalar absTol, typename boost::enable_if_c<Eigen::NumTraits<Scalar>::IsComplex >::type* dummy=0){using namespace std; return abs(num)<=absTol; }
	
	static MatrixBaseT pruned(const MatrixBaseT& a, double absTol=1e-6){ // typename MatrixBaseT::Scalar absTol=1e-6){
		MatrixBaseT ret(MatrixBaseT::Zero(a.rows(),a.cols()));
		for(Index c=0;c<a.cols();c++){ for(Index r=0;r<a.rows();r++){ if(!prune_element(a(c,r),absTol)) ret(c,r)=a(c,r); } }
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
	
	// 4-vector
	template<typename VectorT2, class PyClass> static void visit_special_sizes(PyClass& cl, typename boost::enable_if_c<VectorT2::RowsAtCompileTime==4>::type* dummy=0){
		cl
		.def(py::init<typename VectorT2::Scalar,typename VectorT2::Scalar,typename VectorT2::Scalar>((py::arg("v0"),py::arg("v1"),py::arg("v2"),py::arg("v3"))))
		;
	}

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


	static VectorT dyn_Ones(Index size){ return VectorT::Ones(size); }
	static VectorT dyn_Zero(Index size){ return VectorT::Zero(size); }
	static VectorT dyn_Random(Index size){ return VectorT::Random(size); }
	static VectorT Unit(Index ix){ IDX_CHECK(ix,(Index)Dim); return VectorT::Unit(ix); }
	static VectorT dyn_Unit(Index size,Index ix){ IDX_CHECK(ix,size); return VectorT::Unit(size,ix); }
	static bool dyn(){ return Dim==Eigen::Dynamic; }
	static Index __len__(){ assert(!dyn()); return Dim; }
	static Index dyn__len__(const VectorT& self){ return self.size(); }
	static void resize(VectorT& self, Index size){ self.resize(size); }
	static Scalar dot(const VectorT& self, const VectorT& other){ return self.dot(other); }
	static CompatMatrixT outer(const VectorT& self, const VectorT& other){ return self*other.transpose(); }
	static CompatMatrixT asDiagonal(const VectorT& self){ return self.asDiagonal(); }
	static Scalar get_item(const VectorT& self, Index ix){ IDX_CHECK(ix,dyn()?(Index)self.size():(Index)Dim); return self[ix]; }
	static void set_item(VectorT& self, Index ix, Scalar value){ IDX_CHECK(ix,dyn()?(Index)self.size():(Index)Dim); self[ix]=value; }
	struct VectorPickle: py::pickle_suite{
		static py::tuple getinitargs(const VectorT& x){
			// if this fails, add supported size to the switch below
			BOOST_STATIC_ASSERT(Dim==2 || Dim==3 || Dim==4 || Dim==6 || Dim==Eigen::Dynamic);
			switch((Index)Dim){
				case 2: return py::make_tuple(x[0],x[1]);
				case 3: return py::make_tuple(x[0],x[1],x[2]);
				case 4: return py::make_tuple(x[0],x[1],x[2],x[3]);
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
		for(Index i=0; i<self.size(); i++) oss<<(i==0?"":(((i%3)!=0 || pad>0)?",":", "))<<num_to_string(self.row(i/self.cols())[i%self.cols()],/*pad*/pad);
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
		// matrix*matrix product
		.def("__mul__",&MatrixVisitor::__mul__).def("__imul__",&MatrixVisitor::__imul__)
		// matrix*vector product
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
		.def("__len__",&MatrixVisitor::dyn__len__)
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
		.def("inverse",&MatrixVisitor::inverse,"Return inverted matrix.");
		// decompositions are only meaningful on non-complex numbers
		visit_if_decompositions_meaningful<Scalar,PyClass>(cl);
	}
	// for complex numbers, do nothing
	template<typename Scalar, class PyClass> static	void visit_if_decompositions_meaningful(PyClass& cl, typename boost::enable_if_c<Eigen::NumTraits<Scalar>::IsComplex >::type* dummy = 0){ /* do nothing */ }
	// for non-complex numbers, define decompositions
	template<typename Scalar, class PyClass> static	void visit_if_decompositions_meaningful(PyClass& cl, typename boost::disable_if_c<Eigen::NumTraits<Scalar>::IsComplex >::type* dummy = 0){
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


	static MatrixT dyn_Ones(Index rows, Index cols){ return MatrixT::Ones(rows,cols); }
	static MatrixT dyn_Zero(Index rows, Index cols){ return MatrixT::Zero(rows,cols); }
	static MatrixT dyn_Random(Index rows, Index cols){ return MatrixT::Random(rows,cols); }
	static MatrixT dyn_Identity(Index rows, Index cols){ return MatrixT::Identity(rows,cols); }
	static typename MatrixT::Index dyn__len__(MatrixT& a){ return a.rows(); }
	static typename MatrixT::Index __len__(){ return MatrixT::RowsAtCompileTime; }
	static MatrixT Identity(){ return MatrixT::Identity(); }
	static MatrixT transpose(const MatrixT& m){ return m.transpose(); }
	static CompatVectorT diagonal(const MatrixT& m){ return m.diagonal(); }
	static MatrixT* fromDiagonal(const CompatVectorT& d){ MatrixT* m(new MatrixT); *m=d.asDiagonal(); return m; }
	static void resize(MatrixT& self, Index rows, Index cols){ self.resize(rows,cols); }
	static CompatVectorT get_row(const MatrixT& a, Index ix){ IDX_CHECK(ix,a.rows()); return a.row(ix); }
	static void set_row(MatrixT& a, Index ix, const CompatVectorT& r){ IDX_CHECK(ix,a.rows()); a.row(ix)=r; }
	static Scalar get_item(const MatrixT& a, py::tuple _idx){ Index idx[2]; Index mx[2]={a.rows(),a.cols()}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); return a(idx[0],idx[1]); }
	static void set_item(MatrixT& a, py::tuple _idx, const Scalar& value){ Index idx[2]; Index mx[2]={a.rows(),a.cols()}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); a(idx[0],idx[1])=value; }

	static MatrixT __imul__(MatrixT& a, const MatrixT& b){ a*=b; return a; };
	static MatrixT __mul__(const MatrixT& a, const MatrixT& b){ return a*b; }
	static CompatVectorT __mul__vec(const MatrixT& m, const CompatVectorT& v){ return m*v; }
	// float matrices only
	static MatrixT inverse(const MatrixT& m){ return m.inverse(); }
	static MatrixT __div__(const MatrixT& a, const MatrixT& b){ return a/b; }
	// static void __idiv__(MatrixT& a, const MatrixT& b){ a/=b; };
	static CompatVectorT row(const MatrixT& m, Index ix){ IDX_CHECK(ix,m.rows()); return m.row(ix); }
	static CompatVectorT col(const MatrixT& m, Index ix){ IDX_CHECK(ix,m.cols()); return m.col(ix); }

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
			for(Index r=0; r<m.rows(); r++){
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
			switch((Index)Dim){
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
		// those return internal references, which is what we want
		.add_property("min",py::make_function(&AabbVisitor::min,py::return_internal_reference<>())) 
		.add_property("max",py::make_function(&AabbVisitor::max,py::return_internal_reference<>()))
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
	static VectorType& min(Box& self){ return self.min(); }
	static VectorType& max(Box& self){ return self.max(); }
	static VectorType center(const Box& self){ return self.center(); }
	static VectorType sizes(const Box& self){ return self.sizes(); }
	struct BoxPickle: py::pickle_suite{
		static py::tuple getinitargs(const Box& x){ return py::make_tuple(x.min(),x.max()); }
	};
	static Index len(){ return 2; }
	// getters and setters 
	static Scalar get_item(const Box& self, py::tuple _idx){ Index idx[2]; Index mx[2]={2,Box::AmbientDimAtCompileTime}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); if(idx[0]==0) return self.min()[idx[1]]; return self.max()[idx[1]]; }
	static void set_item(Box& self, py::tuple _idx, Scalar value){ Index idx[2]; Index mx[2]={2,Box::AmbientDimAtCompileTime}; IDX2_CHECKED_TUPLE_INTS(_idx,mx,idx); if(idx[0]==0) self.min()[idx[1]]=value; else self.max()[idx[1]]=value; }
	static VectorType get_minmax(const Box& self, Index idx){ IDX_CHECK(idx,2); if(idx==0) return self.min(); return self.max(); }
	static void set_minmax(Box& self, Index idx, const VectorType& value){ IDX_CHECK(idx,2); if(idx==0) self.min()=value; else self.max()=value; }
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
		.def("angularDistance",&QuaternionVisitor::angularDistance)
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
		.def("slerp",&QuaternionVisitor::slerp,(py::arg("t"),py::arg("other")))
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
	static QuaternionT* fromAxisAngle(const CompatVec3& axis, const Scalar& angle){ QuaternionT* ret=new QuaternionT(AngleAxisT(angle,axis)); ret->normalize();  return ret; }
	static QuaternionT* fromAngleAxis(const Scalar& angle, const CompatVec3& axis){ QuaternionT* ret=new QuaternionT(AngleAxisT(angle,axis)); ret->normalize(); return ret; }
	static QuaternionT* fromTwoVectors(const CompatVec3& u, const CompatVec3& v){ QuaternionT* q(new QuaternionT); q->setFromTwoVectors(u,v); return q; }

	// those must be wrapped since "other" is declared as QuaternionBase<OtherDerived>; the type is then not inferred when using .def
	static QuaternionT slerp(const QuaternionT& self, const Real& t, const QuaternionT& other){ return self.slerp(t,other); }
	static Real angularDistance(const QuaternionT& self, const QuaternionT& other){ return self.angularDistance(other); }

	struct QuaternionPickle: py::pickle_suite{static py::tuple getinitargs(const QuaternionT& x){ return py::make_tuple(x.w(),x.x(),x.y(),x.z());} };
	static QuaternionT Identity(){ return QuaternionT::Identity(); }
	static Vector3r Rotate(const QuaternionT& self, const Vector3r& u){ return self*u; }
	static py::tuple toAxisAngle(const QuaternionT& self){ AngleAxisT aa(self); return py::make_tuple(aa.axis(),aa.angle());}
	static py::tuple toAngleAxis(const QuaternionT& self){ AngleAxisT aa(self); return py::make_tuple(aa.angle(),aa.axis());}
	static CompatVec3 toRotationVector(const QuaternionT& self){ AngleAxisT aa(self); return aa.angle()*aa.axis();}
	static void setFromTwoVectors(QuaternionT& self, const Vector3r& u, const Vector3r& v){ self.setFromTwoVectors(u,v); /*return self;*/ }

	static bool __eq__(const QuaternionT& u, const QuaternionT& v){ return u.x()==v.x() && u.y()==v.y() && u.z()==v.z() && u.w()==v.w(); }
	static bool __ne__(const QuaternionT& u, const QuaternionT& v){ return !__eq__(u,v); }
	static CompatVecX __sub__(const QuaternionT& a, const QuaternionT& b){ CompatVecX r(4); r<<a.w()-b.w(),a.x()-b.x(),a.y()-b.y(),a.z()-b.z(); return r; }

	static Scalar __getitem__(const QuaternionT & self, Index idx){ IDX_CHECK(idx,4); if(idx==0) return self.x(); if(idx==1) return self.y(); if(idx==2) return self.z(); return self.w(); }
	static void __setitem__(QuaternionT& self, Index idx, Real value){ IDX_CHECK(idx,4); if(idx==0) self.x()=value; else if(idx==1) self.y()=value; else if(idx==2) self.z()=value; else if(idx==3) self.w()=value; }
	static string __str__(const py::object& obj){
		const QuaternionT& self=py::extract<QuaternionT>(obj)();
		AngleAxisT aa(self);
		return string(object_class_name(obj)+"((")+num_to_string(aa.axis()[0])+","+num_to_string(aa.axis()[1])+","+num_to_string(aa.axis()[2])+"),"+num_to_string(aa.angle())+")";
	}
	static Index __len__(){return 4;}
};


