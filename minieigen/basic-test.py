import sys
sys.path=['.']+sys.path
from minieigen import *
print Vector3(0,0,0)==Vector3.Zero
print Vector3(0,0,0)!=Vector3.Zero
print Vector3(0,0,0)!=Vector3(0,0,0)
print Vector3(1,2,3).sum()
print Vector3.Random()
m=Matrix3.Random()
m.inverse()


box=AlignedBox3((1,2,3),(4,5,6))
box2=AlignedBox3((0,0,0),(5,5,5))
print box, box.volume()
print box.intersection(box2)
print box.contains((3,3,3))

m3=Matrix3(0,1,2, 3,4,5, 6,7,8)
v3=Vector3(1,2,3)
v3c=Vector3c(1,5j,1-3j)
v6=Vector6(0,1,2, 3,4,5)
v6_=Vector6((0,1,2),(3,4,5))
print v3*1.22
print 1.22*v3
print m3,v3,v3c,v6,v6_
print Vector6c.Random
print Matrix6.Random
print VectorX.Random(7)
print MatrixX.Random(4,4)
print Vector2.UnitY, Vector3.UnitZ, Vector3.UnitX
print Vector3.Unit(2)
print MatrixXc.Random(7,7)
print VectorX.Ones(11)
