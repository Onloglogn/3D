from numpy import array, cos, sin, floor, ceil
from compare import *
from scipy.spatial.transform import Rotation as rot
fname = input()

eps = float(input())
x, y, z = tuple(map(float, input().split()))
faces, vidx, idxv, vf, Ns = read_binary(fname)
R = rot.from_rotvec([x, y, z]).as_matrix()
idxv = array(idxv).T
idxv = R.dot(idxv)
idxv = idxv.T.tolist()
x=1e6
y=1e6
z=1e6
for e in idxv:
    x=min(e[0], x)
    y=min(e[1], y)
    z=min(e[2], z)
idxv = [[e[0]-x, e[1]-y, e[2]-z] for e in idxv]
u=dict()
d=dict()
for e in idxv:
    x, y= int(floor(e[0]/eps)), int(floor(e[1]/eps))
    if u.get((x, y), None)==None:
        u[(x, y)]=e[2]
        d[(x, y)]=e[2]
    else:
        u[(x, y)]=max(u[(x, y)], e[2])
        d[(x, y)]=min(d[(x, y)], e[2])
a=0
b=0
c=1e6
for e in d:
    if d[e]<c:
        a=e[0]
        b=e[1]
        c=d[e]
print(a, b)
K = len(list(d.keys()))
print(K)
for e in u:
    print(e[0], e[1], round(u[e], 6))
print(K)
for e in d:
    print(e[0], e[1], round(d[e], 6))
