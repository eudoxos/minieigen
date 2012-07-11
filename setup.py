#from setuptools import setup,Extension
from distutils.core import setup,Extension
setup(name='minieigen',
	version='0.2',
	description='',
	ext_modules=[Extension('miniEigen',
		sources=['src/miniEigen.cpp',
			'src/double-conversion/bignum.cc',
			'src/double-conversion/bignum-dtoa.cc',
			'src/double-conversion/cached-powers.cc',
			'src/double-conversion/diy-fp.cc',
			'src/double-conversion/double-conversion.cc',
			'src/double-conversion/fast-dtoa.cc',
			'src/double-conversion/fixed-dtoa.cc',
			'src/double-conversion/strtod.cc'
		],
		libraries=['boost_python'],
		include_dirs=['/usr/include/eigen3'],
		define_macros=[('EIGEN_DONT_ALIGN',None)]
	)],
	install_requires=['distribute']
)

