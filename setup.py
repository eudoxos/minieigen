# encoding: utf-8
from distutils.core import setup,Extension
import glob
setup(name='minieigen',
	version='0.3-2',
	author='Václav Šmilauer',
	author_email='eu@doxos.eu',
	url='http://www.launchpad.net/minieigen',
	description='Wrap parts of Eigen3, c++ library for basic math and geometry.',
	long_description='''
A small wrapper for core parts of EIgen, c++ library for linear algebra. It is mainly useful for inspecting c++ code which already uses eigen and boost::python. Supported types are Vectors (2,3,6 and dynamic-sized with integer and floating-point values), Matrices (3x3, 6x6 and dynamic-sized with floating-point values), Quaternions and 2d and 3d AlignedBox's. Numerous methods are wrapped and the original API of Eigen is followed. The code compiles with a c++99 compiler.
''',
	classifiers=[
		'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
		'Programming Language :: Python',
		'Programming Language :: Python :: 3',
		'Operating System :: OS Independent',
		'Topic :: Scientific/Engineering :: Mathematics',
		'Intended Audience :: Science/Research',
		'Development Status :: 4 - Beta'
	],
	ext_modules=[Extension('minieigen',
		sources=['minieigen/minieigen.cpp',
			'minieigen/double-conversion/bignum.cc',
			'minieigen/double-conversion/bignum-dtoa.cc',
			'minieigen/double-conversion/cached-powers.cc',
			'minieigen/double-conversion/diy-fp.cc',
			'minieigen/double-conversion/double-conversion.cc',
			'minieigen/double-conversion/fast-dtoa.cc',
			'minieigen/double-conversion/fixed-dtoa.cc',
			'minieigen/double-conversion/strtod.cc'
		],
		libraries=['boost_python'],
		include_dirs=['/usr/include/eigen3','minieigen'],
		define_macros=[('EIGEN_DONT_ALIGN',None)]
	)],
	data_files=[('headers',glob.glob('minieigen/double-conversion/*.h'))],
	#install_requires=['distribute'],
)

