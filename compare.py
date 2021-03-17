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
