import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg

def MVAscssSolve(NStations, D, c, N, Z):
    Rk  = np.zeros(NStations)
    Nk  = np.zeros(NStations)

    for n in range(1, N+1):
        R = 0

        for k in range(0, NStations):
            if (c[k] <= 0.0):
                Rk[k] = D[k]
            else:
                Rk[k] = D[k] * (1 + Nk[k])
            R = R + Rk[k]

        X = n / (R + Z)

        for k in range(0, NStations):
            Nk[k] = X * Rk[k]

    return {'X':X, 'R':R, 'Nk':Nk, 'Rk':Rk}
    