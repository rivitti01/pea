import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg

def OPENscssSolve(NStations, D, c, l):
    Rk  = np.zeros(NStations)
    Nk  = np.zeros(NStations)
    Uk  = np.zeros(NStations)

    R = 0

    for k in range(0, NStations):
        U = l * D[k]
        Uk[k] = U
        if (c[k] <= 0.0):
            Rk[k] = D[k]
        else:
            Rk[k] = D[k] / (1 - U)
        R = R + Rk[k]

    X = l

    for k in range(0, NStations):
        Nk[k] = X * Rk[k]
    
    #print(X, R, Rk, Nk)
    return {'X':X, 'R':R, 'Nk':Nk, 'Rk':Rk, 'Uk':Uk}
