cdef extern from "math.h":
    double sin(double x)
    double cos(double x)
    double floor(double arg)
    double ceil(double arg)
    double sqrt(double arg)
    double atan2(double x, double y)
cdef double N=1<<15
cdef double pi=3.141592653589793
cpdef double dot(list a, list b):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
cpdef list sub(list a, list b):
    return [a[i]-b[i] for i in range(len(a))]
cpdef list add(list a, list b):
    return [a[i]+b[i] for i in range(len(a))]
cpdef list mul(list a, double b):
    return [el*b for el in a]
cpdef list div(list a, double b):
    return [el/b for el in a]
cpdef double dott(tuple a, tuple b):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
cpdef tuple cross(tuple a, tuple b):
    return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
cpdef tuple subt(tuple a, tuple b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
cpdef tuple addt(tuple a, tuple b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
cpdef tuple mult(tuple a, double b):
    return (a[0]*b, a[1]*b, a[2]*b)
cpdef tuple divt(tuple a, double b):
    return (a[0]/b, a[1]/b, a[2]/b)
cpdef list translate(list idxv, tuple v):
    return [addt(idxv[k], v) for k in range(len(idxv))]
cpdef double vnorm(tuple a):
    return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
cpdef double norm(tuple a):
    return a[0]*a[0]+a[1]*a[1]+a[2]*a[2]
cpdef int lg(double n):#there are faster ways with binary search, but unnecessary
    s=0
    while (n>(1<<s)):
        s+=1
    return s
cpdef list mxmn(list idxv):
    cdef float xmn=N, xmx=-N
    cdef float ymn=N, ymx=-N
    cdef float zmn=N, zmx=-N
    for v in idxv:
        xmx = max(xmx, v[0])
        xmn = min(xmn, v[0])
        ymx = max(ymx, v[1])
        ymn = min(ymn, v[1])
        zmx = max(zmx, v[2])
        zmn = min(zmn, v[2])
    return [xmx, xmn, ymx, ymn, zmx, zmn]
cpdef align(list idxv):
    cdef list l = mxmn(idxv)
    return translate(idxv, (-l[1], -l[3], -l[5]))

cpdef tuple normal(list f):
    return cross(subt(f[2], f[0]), subt(f[1], f[0]))
cpdef tuple matrixmul(tuple vertex, list m):
    return (dott(m[0], vertex), dott(m[1], vertex), dott(m[2], vertex))
cpdef list matrixmulv(list idxv, list m):
    return [matrixmul(idxv[k], m) for k in range(len(idxv))]
cpdef tuple R3(tuple vertex, tuple xyz):#something is wrong here idk what
    cdef double x= xyz[0]
    cdef double y= xyz[1]
    cdef double z= xyz[2]
    cdef tuple R0 = (cos(x)*cos(y)*cos(z)-sin(x)*sin(z),\
            -cos(x)*cos(y)*sin(z)-sin(x)*cos(z), cos(x)*sin(y))
    cdef tuple R1=(sin(x)*cos(y)*cos(z)+cos(x)*sin(z),\
            -sin(x)*cos(y)*sin(z)+cos(x)*cos(z),\
            sin(x)*sin(y))
    cdef tuple R2 = (-sin(y)*sin(z), sin(y)*sin(z), cos(y))
    return (dott(R0, vertex), dott(R1, vertex), dott(R2, vertex))
cpdef tuple R3inv(tuple vertex, tuple xyz):
    cdef double x= xyz[0]
    cdef double y= xyz[1]
    cdef double z= xyz[2]
    cdef tuple R0 = (cos(x)*cos(y)*cos(z)-sin(x)*sin(z),\
            sin(x)*cos(y)*cos(z)+cos(x)*sin(z), -sin(y)*cos(z))
    cdef tuple R1=(-cos(x)*cos(y)*sin(z)-sin(x)*cos(z),\
            -sin(x)*cos(y)*sin(z)+cos(x)*cos(z),\
            sin(y)*sin(z))
    cdef tuple R2 = (cos(x)*sin(y), sin(x)*sin(y), cos(y))
    return (dott(R0, vertex), dott(R1, vertex), dott(R2, vertex))
cpdef list R3v(list idxv, tuple xyz):
    if xyz[0]==0 and xyz[1]==0 and xyz[2]==0:
        return idxv[:]
    return [R3(idxv[k], xyz) for k in range(len(idxv))]
cpdef list R3invv(list idxv, tuple xyz):
    if xyz[0]==0 and xyz[1]==0 and xyz[2]==0:
        return idxv[:]
    return [R3inv(idxv[k], xyz) for k in range(len(idxv))]
cpdef list perm(list idxv, tuple abc):
    return [(v[abc[0]], v[abc[1]], v[abc[2]]) for v in idxv]
cpdef double det(list vertex):
    cdef tuple a = vertex[0]
    cdef tuple b = vertex[1]
    cdef tuple c = vertex[2]
    return (a[0]*b[1]*c[2]+b[0]*c[1]*a[2]+c[0]*a[1]*b[2]\
            -c[0]*b[1]*a[2]-b[0]*a[1]*c[2]-a[0]*c[1]*b[2])
cpdef double faceVol(list faces, list idxv):
    cdef double s=0
    for f in faces:
        s+=det([idxv[f[0]], idxv[f[1]], idxv[f[2]]])
    return abs(s/6)

cpdef int inside_tetra(p, t):
    cdef tuple n1= cross(sub(t[2],t[0]), sub(t[1],t[0]))
    cdef tuple n2= cross(sub(t[1],t[0]), sub(t[3],t[0]))
    cdef tuple n3= cross(sub(t[3], t[0]), sub(t[2], t[0]))
    cdef tuple n4= cross(sub(t[2], t[1]), sub(t[3], t[1]))
    cdef double a= dott(n1, sub(p, t[0]))
    cdef double b= dott(n2, sub(p, t[0]))
    cdef double c= dott(n3, sub(p, t[0]))
    cdef double d= dott(n4, sub(p, t[1]))
    return (a>0 and b>0 and c>0 and d>0) or (a<0 and b<0 and c<0 and d<0)
cpdef list circumsphere(list t):
    cdef list ff = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]
    cdef list n = [norm(t[0]),norm(t[1]),norm(t[2]), norm(t[3])]
    cdef double dx = sum([det([(n[ff[i][0]], n[ff[i][1]], n[ff[i][2]]),\
            (t[ff[i][0]][1], t[ff[i][1]][1], t[ff[i][2]][1]), \
            (t[ff[i][0]][2], t[ff[i][1]][2], t[ff[i][2]][2])])*(1 if i%2 else -1) for i in range(4)])
    cdef double dy = -sum([det([(n[ff[i][0]], n[ff[i][1]], n[ff[i][2]]),\
            (t[ff[i][0]][0], t[ff[i][1]][0], t[ff[i][2]][0]), \
            (t[ff[i][0]][2], t[ff[i][1]][2], t[ff[i][2]][2])])*(1 if i%2 else -1) for i in range(4)])
    cdef double dz = sum([det([(n[ff[i][0]], n[ff[i][1]], n[ff[i][2]]),\
            (t[ff[i][0]][0], t[ff[i][1]][0], t[ff[i][2]][0]), \
            (t[ff[i][0]][1], t[ff[i][1]][1], t[ff[i][2]][1])])*(1 if i%2 else -1) for i in range(4)])
    cdef double u = det([t[1], t[2], t[3]])
    cdef double v = det([t[0], t[2], t[3]])
    cdef double w = det([t[0], t[1], t[3]])
    cdef double x = det([t[0], t[1], t[2]])
    cdef double a = -u+v-w+x
    cdef double c = n[0]*u-n[1]*v+n[2]*w-n[3]*x
    return [divt((dx, dy, dz), 2*a), abs(sqrt(dx*dx+dy*dy+dz*dz-4*a*c)/(2*a))]
