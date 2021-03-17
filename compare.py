from struct import pack, unpack
from cc import *
from numpy import linalg as LA
from numpy import ceil, floor
from random import random
def read(file_name, sort=True):
    f = open(file_name, 'r')
    r=0
    faces = list()
    face=list()
    N = None
    vidx = dict()
    idxv = list()
    vf = dict()# list of vertices connecting to faces...
    Ns = list()
    idx=0
    for line in f:
        if "normal" in line:
            N = tuple(map(float, [el for el in line.split(' ') if el != ''][2:5]))
        if "vertex" in line:
            if r==0:
                r=1
            vertex = list(map(float,  \
                    [el for el in line.split(" ") if el!=''][1:4]))
            vertex = [round(e, 6) for e in vertex]
            vertex = tuple(vertex)
            if vidx.get(vertex, None)==None:
                vidx[vertex]=idx
                idxv.append(vertex)
                face.append(idx)
                idx+=1
            else:
                face.append(vidx[vertex])
        else:
            if r==1:
                faces.append(face)
                Ns.append(N)
                for v in face:
                    if vf.get(v, None)==None:
                        vf[v] = [len(faces)-1]
                    else:
                        vf[v].append(len(faces)-1)
                face = list()
                r=0
    return faces, vidx, idxv, vf, Ns
def read_binary(filename):
    f = open(filename, 'rb')
    s = f.read()
    f.close()
    faces = list()
    n = unpack('i', s[80:84])[0]
    vidx = dict()
    idxv = list()
    vf = dict()
    Ns = list()
    idx = 0
    for i in range(84, 84+50*n, 50):
        face = list()
        x, y, z = unpack("3f", s[i:i+12])
        normal = (x, y, z)
        for j in range(1, 4):
            x, y, z = unpack("3f", s[i+12*j:i+12*(j+1)])
            x, y, z = round(x, 6), round(y, 6), round(z, 6)
            if vidx.get((x, y, z), None)==None:
                vidx[(x, y, z)]=idx
                idxv.append((x, y, z))
                face.append(idx)
                idx+=1
            else:
                face.append(vidx[(x, y, z)])
        faces.append(face)
        Ns.append(normal)
        for v in face:
            if vf.get(v, None)==None:
                vf[v] = [len(faces)-1]
            else:
                vf[v].append(len(faces)-1)
    xyz = mxmn(idxv)
    abc = [(xyz[0]-xyz[1], 0), (xyz[2]-xyz[3], 1), (xyz[4]-xyz[5], 2)]
    abc = sorted(abc, key=lambda x:-x[0])
    abc = tuple([e[1] for e in abc])
    idxv = perm(idxv, abc)
    return faces, vidx, idxv, vf, Ns
def orient_faces(faces, idxv, Ns):
    for i in range(len(faces)):
        v = [idxv[faces[i][j]] for j in range(3)]
        if len(Ns[i])>0 and dott(normal(v), Ns[i])<0:
            faces[i] = (faces[i][0], faces[i][2], faces[i][1])
    return faces
def create_fake(filename, faces, idxv):
    O = b'\x00'
    f = open(filename, 'wb+')
    f.write(80*O)
    f.write(pack('i', len(faces)))
    for e in faces:
        f.write(pack('f', e[1][0]))
        f.write(pack('f', e[1][1]))
        f.write(pack('f', e[1][2]))
        for e2 in e[0]:
            for i in range(3):
                f.write(pack('f', idxv[e2][i]))
        f.write(O*2)
    f.close()
    return True
def mxid(v):
    return sorted([[v[i], i] for i in range(len(v))])[-1][1]
class BoundingBox:
    def __init__(self, center, h, v):
        self.center = center
        self.h = h
        self.v = v
def obb(faces, idxv):
    m, C = covariance(faces, idxv)
    u, v = LA.eig(C)
    v = [tuple(e) for e in v]
    v = [divt(e, vnorm(e)) for e in v]
    mm = mxmnB(idxv, v)#max and min in each direction
    return BoundingBox(addt(mult(v[0], (mm[0]+mm[1])/2), addt(mult(v[1], \
            (mm[2]+mm[3])/2), mult(v[2], (mm[4]+mm[5])/2))), \
            ((mm[0]-mm[1])/2, (mm[2]-mm[3])/2, (mm[4]-mm[5])/2), v), u#position of center 
            #and height of box in the direction of respective basis and eigen values and vectors
def obbtree(faces, idxv, depth):
    b, u = obb(faces, idxv)
    tree = {1:b}
    gtei = {1:mxid(u)}#greatest eigenvalue
    #to the left add a 0, else add one to go to child
    ba = bar(faces, idxv)
    v = [[ba[i], 1] for i in range(len(ba))]
    for i in range(1, depth+1):
        for j in range(len(v)):
            r=0
            if dott(subt(v[j][0], tree[v[j][1]].center), tree[v[j][1]].v[gtei[v[j][1]]])>0:
                r=1
            v[j][1]=v[j][1]*2+r
        nd = {j:list() for j in range(1<<i, 1<<(i+1))}
        for j in range(len(v)):
            nd[v[j][1]].append(j)
        for j in range(1<<i, 1<<(i+1)):
            face = [faces[e] for e in nd[j]]
            if len(face)>0:
                b, u = obb(face, idxv)
                gtei[j] = mxid(b.h)
                tree[j] = b
            else:
                tree[j] = None
                gtei[j] = 0
    return tree
def build_graph(faces):
    edg2f = dict()
    graph = {i:list() for i in range(len(faces))}
    for i in range(len(faces)):
        per = [(j, k) for k in range(len(faces[i])) for j in range(len(faces[i])) if j<k]
        for p in per:
            ed = (faces[i][p[0]], faces[i][p[1]]) if faces[i][p[0]]<faces[i][p[1]] else (faces[i][p[1]], faces[i][p[0]])
            if edg2f.get(ed, None)==None:
               edg2f[ed]=[i]
            else:
                edg2f[ed].append(i)
    for ed in edg2f.keys():
        for i in range(len(edg2f[ed])):
            for j in range(i+1, len(edg2f[ed])):
                graph[edg2f[ed][i]].append([edg2f[ed][j], ed])
                graph[edg2f[ed][j]].append([edg2f[ed][i], ed])
    return graph
def HM(idxv, eps):#height matrix
    dc1, dc2= dict(), dict()
    for v in idxv:
        a, b = int(floor(v[0]/eps)), int(floor(v[1]/eps))
        if dc1.get((a, b), None)==None:
            dc1[(a, b)] = v[2]
        else:
            dc1[(a, b)] = max(v[2], dc1[(a, b)])
        if dc2.get((a, b), None)==None:
            dc2[(a, b)] = v[2]
        else:
            dc2[(a, b)] = min(v[2], dc2[(a, b)])
    z = dc2[min(list(dc2.keys()), key=lambda x:dc2[x])]
    for k in dc2:
        dc1[k]-=z
        dc2[k]-=z
    return dc1, dc2
