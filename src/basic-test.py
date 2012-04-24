import sys
sys.path=['.']+sys.path
from miniEigen import *
box=AlignedBox3((1,2,3),(4,5,6))
box2=AlignedBox3((0,0,0),(5,5,5))
print box, box.volume()
print box.intersection(box2)
print box.contains((3,3,3))
