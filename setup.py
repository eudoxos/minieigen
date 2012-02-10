#from setuptools import setup,Extension
from distutils.core import setup,Extension
setup(name='minieigen',
	version='0.1',
	description='',
	ext_modules=[Extension('miniEigen',
		sources=['src/miniEigen.cpp'],
		libraries=['boost_python'],
		include_dirs=['/usr/include/eigen3','/usr/include/eigen2'],
	)],
	install_requires=['distribute']
)

