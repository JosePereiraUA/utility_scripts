from pymol import cmd
from pymol.cgo import *
import numpy as np



def hyperballs(objname='hyperballs', selection='sele',
               scale='0.2', radius='1.0', resolution='41', export=False):
    '''
    
    hyperballs objname='hyperballs', selection='sele',
                scale=1.0, radius=1.0, resolution=41, export=False

                
    '''

    def rotation_about_vector(matrix, angle, axis):
        theta = 0.5*angle
        q0 = np.cos(theta)
        q1,q2,q3 = np.sin(theta)*axis/np.linalg.norm(axis)
        Rmat = np.array([[1-2*q2*q2-2*q3*q3,   2*q1*q2-2*q0*q3,   2*q1*q3+2*q0*q2],
                        [   2*q2*q1+2*q0*q3, 1-2*q3*q3-2*q1*q1,   2*q2*q3-2*q0*q1],
                        [   2*q3*q1-2*q0*q2,   2*q3*q2+2*q0*q1, 1-2*q1*q1-2*q2*q2]], float)
        return np.dot(Rmat, matrix.T).T

    def fitCircle(a,b,c,d,n):
        Ry = (-a**2 + b**2 - 2*a*c + 2*b*c + d**2)/(2*d)
        Rx = np.sqrt((b+c)**2 - Ry**2)
        Xa = d+a*(Ry-d)/(a+c)
        Xb = b*Ry/(b+c)
        Y = np.linspace(Xa,Xb,n)
        R = -np.sqrt(c*c - (Y-Ry)**2) + Rx
        return Y,R

    def getCone(xbase,xtip,rbase,rtip,cbase,ctip):
        return [CONE,
                xbase[0], xbase[1], xbase[2],
                xtip[0],  xtip[1],  xtip[2],
                rbase, rtip,
                cbase[0], cbase[1], cbase[2],
                ctip[0],  ctip[1],  ctip[2], 1.0, 1.0]

    def getSphere(x,r,c):
        return [COLOR, c[0], c[1], c[2], SPHERE, x[0], x[1], x[2], r]

    def interpolateRGB(rgb0, rgb1, n):
        drgb = rgb1-rgb0
        x = np.linspace(0.0,1.0,n)
        return [rgb0 + xi*drgb for xi in x]


    # validate input
    #-------------------------------------------------------
    RR = float(radius)
    SS = float(scale)
    export=bool(export)
    RESOLUTION = int(resolution)
    if not RESOLUTION%2:
        RESOLUTION+=1
    model  = cmd.get_model(selection, cmd.get_state())

    # get colors
    #-------------------------------------------------------
    pymol.clist = []
    cmd.iterate(selection, 'pymol.clist.append(color)')
    colors = np.array([cmd.get_color_tuple(color) for color in pymol.clist])
    
    # get atomic information
    #-------------------------------------------------------
    coord = np.array([atom.coord for atom in model.atom])
    radii = np.array([atom.vdw   for atom in model.atom])*SS
    bonds = np.array([bond.index for bond in model.bond])

    # get bond vector
    #-------------------------------------------------------
    Rij = coord[bonds[:,0],:]-coord[bonds[:,1],:]
    dij = np.sqrt(np.sum(Rij**2, axis=1))
    Rijn = Rij/dij.repeat(3).reshape(-1,3)

    #-------------------------------------------------------
    # builg CGO object
    #-------------------------------------------------------
    obj = []
    
    # -> spheres
    for x1, rad, color in zip(coord, radii, colors):
        obj.extend(getSphere(x1, rad, color))
        
    # -> hyperboloids
    for n,(i,j) in enumerate(bonds):
        col = interpolateRGB(colors[i], colors[j], RESOLUTION)
        k,r = fitCircle(radii[i], radii[j], RR, dij[n], RESOLUTION)
        for m in range(RESOLUTION-1):
            x1 = coord[j] + k[m  ]*Rijn[n]
            x2 = coord[j] + k[m+1]*Rijn[n]
            obj.extend(getCone(x1,x2,r[m],r[m+1], col[m],col[m+1]))

    # -> load object
    currentview = cmd.get_view(output=1, quiet=1)
    cmd.delete(objname)
    cmd.load_cgo(obj, objname, 1)
    cmd.set_view(currentview)


    if export:
        
        radialRESOLUTION = 20
        Npoints = RESOLUTION*radialRESOLUTION
        theta  = np.linspace(0., 2*np.pi, radialRESOLUTION, False)
        stheta = np.tile(np.sin(theta), RESOLUTION)
        ctheta = np.tile(np.cos(theta), RESOLUTION)

        l1 = np.arange(radialRESOLUTION+1) % radialRESOLUTION
        l2 = l1 + radialRESOLUTION
        faceLoop = np.c_[l1[:-1],l1[1:],l2[1:],l2[:-1]]
        faces = 1 + np.vstack([n+faceLoop for n in range(0,Npoints-radialRESOLUTION,radialRESOLUTION)])

        vtexList = []   # vertex list
        faceList = []   # Face list

        cg = coord.mean(0)
        for n,(i,j) in enumerate(bonds):
            z,r = fitCircle(radii[i], radii[j], RR, dij[n], RESOLUTION)
            zz = np.repeat(z, radialRESOLUTION)
            rr = np.repeat(r, radialRESOLUTION)

            xyz = np.c_[rr*ctheta, rr*stheta, zz]
            rax = np.cross([0.,0.,1.0], Rijn[n])
            angle = np.arccos(Rijn[n,2])
            vtexList.append(coord[j]+rotation_about_vector(xyz, angle,rax))
            faceList.append(faces+n*Npoints)
        
        vtexList = np.reshape(vtexList, (-1,3))-cg
        faceList = np.reshape(faceList, (-1,4))

        form  = 'Transform { translation %.6f %.6f %.6f\n'
        form += '            children Shape { geometry Sphere { radius %.3f } }}'

        vtexStr = ['v %12.6f %12.6f %12.6f'%(x,y,z) for x,y,z in vtexList]
        faceStr = ['f %d %d %d %d'%(i,j,k,l) for i,j,k,l in faceList]
        sphereStr = [form%(x,y,z,r) for r,(x,y,z) in zip(radii,coord-cg)]

        # export *.obj file
        #fout = open('/Users/SergioSantos/hyperballs.obj','w')
        fout = open('hyperballs.obj','w')
        fout.write('\n'.join(vtexStr+faceStr))
        fout.close()

        # export *.wrl file
        #fout = open('/Users/SergioSantos/hyperballs.wrl','w')
        fout = open('hyperballs.wrl','w')
        fout.write('\n'.join(sphereStr))
        fout.close()
        
    return


cmd.extend('hyperballs', hyperballs)




