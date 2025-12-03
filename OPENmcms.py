import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg

def OPENmcmsSolve(NStations, NClasses, D, c, l):
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
    
        if(c[k] > 0):
            Uk[k] = U / c[k]
        else:
            Uk[k] = U
        
        if (c[k] <= 0.0):
            RF = 1
        elif (c[k] >= 2.0):
            Fc = 1
            Sm = 1
            Tr = 1
            for j in range(1, int(c[k])):
                Tr = Tr * U / j
                Fc = Fc * j / U
                Sm = Sm + Tr
            RF = (1 + 1 / (c[k] - U) / (1 + (c[k] - U) / U * Fc * Sm))
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
    return {'X':X, 'R':R, 'Nk':Nk, 'Rk':Rk, 'Uk':Uk * c, 'Rkc':Rkc}