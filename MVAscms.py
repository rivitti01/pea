import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg

def MVAscmsSolve(NStations, D, c, N, Z):
    Rk  = np.zeros(NStations)
    Nk  = np.zeros(NStations)
    Pki = np.zeros((NStations, N+1))
    for k in range(0, NStations):
        Pki[k, 0] = 1

    for n in range(1, N+1):
        R = 0

        for k in range(0, NStations):
            if (c[k] <= 0.0):
                Rk[k] = D[k]
            elif (c[k] >= 2.0):
                Rk[k] = 0
                for j in range(1, n+1):
                    cd = min(float(j), c[k])
                    Rk[k] = Rk[k] + float(j) / cd * Pki[k, j-1]
                Rk[k] = Rk[k] * D[k]
            else:
                Rk[k] = D[k] * (1 + Nk[k])
            R = R + Rk[k]

        X = n / (R + Z)

        for k in range(0, NStations):
            if (c[k] >= 2.0):
                spk = 0
                for j in range(n, 0, -1):
                    cd = min(float(j), c[k])
                    Pki[k, j] = X * D[k] / cd * Pki[k, j-1]
                    spk = spk + Pki[k, j]
                Pki[k, 0] = 1 - spk
            Nk[k] = X * Rk[k]
    
    return {'X':X, 'R':R, 'Nk':Nk, 'Rk':Rk}
