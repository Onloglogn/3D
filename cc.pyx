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
cpdef list mxmnB(list idxv, list b):
    cdef float xmn=N, xmx=-N
    cdef float ymn=N, ymx=-N
    cdef float zmn=N, zmx=-N
    for v in idxv:
        m = dott(v, b[0])
        n = dott(v, b[1])
        o = dott(v, b[2])
        xmx = max(xmx, m)
        xmn = min(xmn, m)
        ymx = max(ymx, n)
        ymn = min(ymn, n)
        zmx = max(zmx, o)
        zmn = min(zmn, o)
    return [xmx, xmn, ymx, ymn, zmx, zmn]
cpdef double D(list a, list b):
    cdef float r=0
    for i in range(len(a)):
        r+=(a[i]-b[i])*(a[i]-b[i])
    return r
cpdef double __sum__(list l):
    cdef double s=0
    for i in range(len(l)):
        s+=l[i]
    return s
cpdef tuple normal(list f):
    return cross(subt(f[2], f[0]), subt(f[1], f[0]))
cpdef tuple matrixmul(tuple vertex, list m):
    return (dott(m[0], vertex), dott(m[1], vertex), dott(m[2], vertex))
cpdef list matrixmulv(list idxv, list m):
    return [matrixmul(idxv[k], m) for k in range(len(idxv))]
cpdef tuple R3(tuple vertex, tuple xyz):
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
cpdef list bar(list faces, list idxv):
    return [divt(addt(idxv[f[0]], addt(idxv[f[1]], idxv[f[2]])), 3) for f in faces]
cpdef tuple mean(list faces, list idxv):
    cdef tuple mu= (0.0, 0.0, 0.0)
    for u in range(0, len(faces)):
        for v in range(3):
            mu = addt(mu, idxv[faces[u][v]])
    mu = divt(mu, 3*len(faces))
    return mu
cpdef tuple pcmean(list idxv):
    cdef tuple mu = (0.0, 0.0, 0.0)
    for i in range(len(idxv)):
        mu = addt(mu, idxv[i])
    return divt(mu, len(idxv))
cpdef list covariance(list faces, list idxv):
    cdef tuple mu= (0.0, 0.0, 0.0)
    cdef list C = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    cdef double a, b, c, d, e, f, g, h, i
    for u in range(0, len(faces)):
        for v in range(3):
            mu = addt(mu, idxv[faces[u][v]])
    mu = divt(mu, 3*len(faces))
    for u in range(0, len(faces)):
        a = idxv[faces[u][0]][0]-mu[0]#
        b = idxv[faces[u][0]][1]-mu[1]
        c = idxv[faces[u][0]][2]-mu[2]
        d = idxv[faces[u][1]][0]-mu[0]#
        e = idxv[faces[u][1]][1]-mu[1]
        f = idxv[faces[u][1]][2]-mu[2]
        g = idxv[faces[u][2]][0]-mu[0]#
        h = idxv[faces[u][2]][1]-mu[1]
        i = idxv[faces[u][2]][2]-mu[2]
        C[0][0]+=a*a+d*d+g*g
        C[1][1]+=b*b+e*e+h*h
        C[2][2]+=c*c+f*f+i*i
        C[0][1]+=a*b+d*e+g*h
        C[0][2]+=a*c+d*f+g*i
        C[1][2]+=b*c+e*f+h*i
    for u in range(3):
        for v in range(3):
            C[u][v] = C[u][v]/(3*len(faces))
    C[1][0]=C[0][1]
    C[2][0]=C[0][2]
    C[2][1]=C[1][2]
    return [mu, C]
cpdef list pccovariance(list idxv):
    cdef tuple mu = pcmean(idxv)
    cdef list C = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    cdef double a, b, c
    for u in range(len(idxv)):
        a = idxv[u][0]
        b = idxv[u][1]
        c = idxv[u][2]
        C[0][0]+=a*a
        C[1][1]+=b*b
        C[2][2]+=c*c
        C[0][1]+=a*b
        C[0][2]+=a*c
        C[1][2]+=b*c
    for u in range(3):
        for v in range(3):
            C[u][v] = C[u][v]/(len(idxv))
    C[1][0]=C[0][1]
    C[2][0]=C[0][2]
    C[2][1]=C[1][2]
    return [mu, C]
cpdef list BBS(vertex, c, n, h):#bounding box search
    cdef list vc = [subt(vertex[i], c) for i in range(len(vertex))]
    return [vertex[i] for i in range(len(vc)) if abs(dott(vc[i], n[0]))<h[0] and abs(dott(vc[i], n[1]))<h[1] and abs(dott(vc[i], n[2]))<h[2]]
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
cpdef list BBT(list idxv, int D):
    cdef list tr = [0]
    cdef list vn = [[idxv[i], 1] for i in range(len(idxv))]
    cdef int n=1
    cdef list l
    
    for i in range(D):
        for j in range(1<<i):
            l = mxmn([e[0] for e in vn if e[1]==n])
            tr.append([[((l[0]+l[1])/2, (l[2]+l[3])/2, )]])
            n+=1
    return []
