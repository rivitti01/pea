import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg

def OPENmcssSolve(NStations, NClasses, D, c, l):
    Rkc  = np.zeros((NStations, NClasses))
    Rk  = np.zeros(NStations)
    Nk  = np.zeros(NStations)
    Uk  = np.zeros(NStations)
    
    R = 0
    X = sum(l)
    
    for k in range(0, NStations):
        U = 0
        for ic in range(0, NClasses):
            U = U + l[ic] * D[k, ic]
    
        Uk[k] = U
        
        if (c[k] <= 0.0):
            RF = 1
        else:
            RF = 1 / (1 - U)
    
        for ic in range(0, NClasses):
            Rkc[k, ic] = D[k, ic] * RF
            Rk[k] = Rk[k] + Rkc[k, ic] * l[ic] / X
        R = R + Rk[k]
    
    for k in range(0, NStations):
        Nk[k] = X * Rk[k]
    
    #print(X, R, Rk, Nk, Uk)
    #print(Rkc)
    return {'X':X, 'R':R, 'Nk':Nk, 'Rk':Rk, 'Uk':Uk, 'Rkc':Rkc}