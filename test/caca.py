import numpy as np



def XYfromNormalVector(n):
    
    siny = -n[0]
    #Pathological case
    if siny*siny > 1.0:
        if siny > 0.0:
            return np.asarray([0.0, np.pi/2.0, 0.0])
        else:
            return np.asarray([0.0, -np.pi/2.0, 0.0])
    y = np.arcsin(siny)
    cosy = np.sqrt(1.0 - siny*siny)    
    cosx = n[2]/cosy
    #We set the range of cosx
    if np.abs(cosx) > 1.0:
        cosxp = 1.0
    y = np.arcsin(siny)  
    x = np.arccos(cosx)
    sinxp = np.sin(x)
    sinxm = -np.sin(x)
    cosxp = np.cos(x)
    cosxm = -np.cos(x)
    cosym = np.cos(y)
    cosyp = -np.cos(y)
    
    if siny > 0:          
        if normalVectorMatches(n, sinxp, siny, cosxp, cosyp):
            x = x
        elif normalVectorMatches(n, sinxp, siny, cosxp, cosym):
            y = np.pi-y
        elif normalVectorMatches(n, sinxp, siny, cosxm, cosyp):
            x = np.pi-x
        elif normalVectorMatches(n, sinxp, siny, cosxm, cosym):
            x = np.pi-x
            y = np.pi-y
        elif normalVectorMatches(n, sinxm, siny, cosxp, cosyp):
            x = -x
        elif normalVectorMatches(n, sinxm, siny, cosxp, cosym):
            x = -x
            y = np.pi-y
        elif normalVectorMatches(n, sinxm, siny, cosxm, cosyp):
            x = np.pi + x
        else:
            x = np.pi + x
            y = np.pi-y
    else:
        if normalVectorMatches(n, sinxp, siny, cosxp, cosyp):
            y = -y
        elif normalVectorMatches(n, sinxp, siny, cosxp, cosym):
            y = np.pi+y
        elif normalVectorMatches(n, sinxp, siny, cosxm, cosyp):
            x = np.pi-x
            y = -y
        elif normalVectorMatches(n, sinxp, siny, cosxm, cosym):
            x = np.pi-x
            y = np.pi+y
        elif normalVectorMatches(n, sinxm, siny, cosxp, cosyp):
            x = -x
            y = -y
        elif normalVectorMatches(n, sinxm, siny, cosxp, cosym):
            x = -x
            y = np.pi+y
        elif normalVectorMatches(n, sinxm, siny, cosxm, cosyp):
            x = np.pi + x
            y = -y
        else:
            x = np.pi + x
            y = np.pi+y
    return np.asarray([x,y,0.0])


def normalVectorMatches(n, sinxm, siny, cosxp, cosyp):
        
    tol = 0.001
    A = makeXYspecialMatrix(sinxm, siny, cosxp, cosyp)
    z = np.asarray([0.0, 0.0, 1.0])
    vz = np.asarray(A.dot(z))[0]
    if np.abs(vz[0]-n[0]) < tol and np.abs(vz[1]-n[1]) < tol and np.abs(vz[2]-n[2]) < tol:   
        return True
    return False


def makeXYspecialMatrix(sinx, siny, cosx, cosy):

    A = np.asmatrix([[cosy, 0, -siny],
                    [-sinx*siny, cosx, -sinx*cosy],
                    [cosx*siny, sinx, cosx*cosy]])
    return A

   

if __name__=='__main__':

    for i in range(0, 5):
        for j in range(0, 5):
            theta = np.pi/5.0 * i
            phi = 2.0*np.pi/5.0 * j
            n = np.asarray([np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)])
            anglenom = XYfromNormalVector(n)       
            print(anglenom)
            matxnom_ = [[1.0, 0, 0],
                        [0.0, np.cos(anglenom[0]), -np.sin(anglenom[0])],
                        [0.0, np.sin(anglenom[0]), np.cos(anglenom[0])]]
            matxnom = np.asmatrix(matxnom_)
            #RotY
            matynom_ = [[np.cos(anglenom[1]), 0.0, -np.sin(anglenom[1])],
                        [0.0, 1.0, 0.0],
                        [np.sin(anglenom[1]), 0.0, np.cos(anglenom[1])]]
            matynom = np.asmatrix(matynom_)
            #RotZ
            matznom_ = [[np.cos(anglenom[2]), -np.sin(anglenom[2]), 0.0],
                        [np.sin(anglenom[2]), np.cos(anglenom[2]), 0.0],
                        [0.0, 0.0, 1.0]]
            matznom = np.asmatrix(matznom_)
            rotnom = matxnom.dot(matynom.dot(matznom))
            znom = np.asarray([0.0, 0.0, 1.0])
            nnew = np.asarray(rotnom.dot(znom))[0]

            print('Real', n, 'Other', nnew)   