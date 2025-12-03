import numpy as np

def MMcKSolve(D, c, K, l, pMax = -1):
    if (pMax == -1) or (pMax > K):
        pMax = K
    
    rho = l * D 

    Pn = np.zeros(K + 1)
    sm = 0.0
    for i in range(0, K + 1):
        if i==0:
            Pn[i] = 1.0
            sm = 1.0
        elif i < c:
            Pn[i] = Pn[i - 1] * rho / i
            sm = sm + Pn[i]
        else:
            Pn[i] = Pn[i - 1] * rho / c
            sm = sm + Pn[i]

    for i in range(0, K + 1):
        Pn[i] = Pn[i] / sm

    N = 0
    U = 0
    for i in range(1, K + 1):
        N = N + i * Pn[i]
        U = U + min(i, c) * Pn[i]
    Uave = U / c
    pB = Pn[K]
    Dr = pB  * l
    X = l - Dr
    R = N / X
    Qn = N - U
    Qtt = R - D

    return {'type': 'M/M/c/K', 'X':X, 'R':R, 'N':N, 'U':U, 'Uave': Uave, 'Qn':Qn, 'Qtt': Qtt, 'Pn': Pn[0:pMax+1],
            'pB': pB, 'Dr': Dr}
