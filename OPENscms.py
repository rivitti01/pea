import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg

def OPENscmsSolve(NStations, D, c, l):
    Rk  = np.zeros(NStations)
    Nk  = np.zeros(NStations)
    Uk  = np.zeros(NStations)

    R = 0

    for k in range(0, NStations):
        U = l * D[k]
        if(c[k] > 0):
            Uk[k] = U / c[k]
        else:
            Uk[k] = U

        if (c[k] <= 0.0):
            Rk[k] = D[k]
        elif (c[k] >= 2.0):
            Fc = 1
            Sm = 1
            Tr = 1
            for j in range(1, int(c[k])):
                Tr = Tr * U / j
                Fc = Fc * j / U
                Sm = Sm + Tr
            Rk[k] = D[k] * (1 + 1 / (c[k] - U) / (1 + (c[k] - U) / U * Fc * Sm))
        else:
            Rk[k] = D[k] / (1 - U)
        R = R + Rk[k]

    X = l

    for k in range(0, NStations):
        Nk[k] = X * Rk[k]

    #print(X, R, Rk, Nk)
    return {'X':X, 'R':R, 'Nk':Nk, 'Rk':Rk, 'Uk':Uk * c}
