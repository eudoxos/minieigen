
minieigen documentation
========================

Overview
--------

.. todo:: Something concise here.

Examples
---------

.. todo:: Some examples of what can be done with :obj:`minieigen`.

Naming conventions
------------------

* Classes are suffixed with number indicating size where it makes sense (it does not make sense for :obj:`minieigen.Quaternion`):

  - :obj:`minieigen.Vector3` is a 3-vector (column vector);
  - :obj:`minieigen.Matrix3` is a 3×3 matrix;
  - :obj:`minieigen.AlignedBox3` is aligned box in 3d;
  - ``X`` indicates dynamic-sized types, such as :obj:`minieigen.VectorX` or :obj:`minieigen.MatrixX`.

* Scalar (element) type is suffixed at the end:

  - nothing is suffixed for floats (:obj:`minieigen.Matrix3`);
  - ``i`` indicates integers (:obj:`minieigen.Matrix3i`);
  - ``c`` indicates complex numbers (:obj:`minieigen.Matrix3c`).

* Methods are named as follows:

  - static methods are upper-case (as in c++), e.g. :obj:`minieigen.Matrix3.Random`;

    - nullary static methods are exposed as properties, if they return a constant (e.g. :obj:`minieigen.Matrix3.Identity`); if they don't, they are exposed as methods (:obj:`minieigen.Matrix3.Random`); the idea is that the necessity to call the method (``Matrix3.Random()``) singifies that there is some computation going on, whereas constants behave like immutable singletons.

  - non-static methods are lower-case (as in c++), e.g. :obj:`minieigen.Matrix3.inverse`.

* Return types:

  - methods modifying the instance in-place return ``None`` (e.g. :obj:`minieigen.Vector3.normalize`); some methods in c++ (e.g. `Quaternion::setFromTwoVectors <http://eigen.tuxfamily.org/dox-devel/classEigen_1_1QuaternionBase.html#ac35460294d855096e9b687cadf821452>`__) both modify the instance *and* return the reference to it, which we don't want to do in Python (:obj:`minieigen.Quaternion.setFromTwoVectors`);
  - methods returning another object (e.g. :obj:`minieigen.Vector3.normalized`) do not modify the instance;
  - methods returning (non-const) references return by value in python

Limitations
-----------
* Type conversions (e.g. float to complex) are not supported.
* Methods returning references in c++ return values in Python (so e.g. ``Matrix3().diagonal()[2]=0`` would zero the last diagonal element in c++ but not in Python).
* Many methods are not wrapped, though they are fairly easy to add.
* Conversion from 1-column :obj:`MatrixX` to :obj:`VectorX` is not automatic in places where the algebra requires it.
* Alignment of matrices is not supported (therefore Eigen cannot vectorize the code well); it might be a performance issue in some cases; c++ code interfacing with minieigen (in a way that c++ values can be *set* from Python) **must** compile with ``EIGEN_DONT_ALIGN``, otherwise there might be crashes at runtime when vector instructions receive unaligned data. It seems that alignment is difficult to do with boost::python.
* Proper automatic tests are missing.

Links
------
* http://eigen.tuxfamily.org (Eigen itself)
* http://www.launchpad.net/minieigen (upstream repository, bug reports, answers)
* https://pypi.python.org/pypi/minieigen (Python package index page, used by ``easy_install``)
* packages:

  - `Debian <http://packages.debian.org/search?keywords=minieigen>`__
  - Ubuntu: `distribution <http://packages.ubuntu.com/search?keywords=minieigen&searchon=names&suite=all&section=all>`__, `PPA <https://code.launchpad.net/~eudoxos/+archive/minieigen>`__

Documentation
-------------

* :ref:`genindex`
* :ref:`search`

.. automodule:: minieigen
  :members:
  :undoc-members:
