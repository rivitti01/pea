import numpy as np

def MMcSolve(D, c, l, pMax = 10):
    X = l
    U = l * D
    Uave = U / c

    
    Fc = 1
    Sm = 1
    Tr = 1
    R = 0
    
    for j in range(1, int(c)):
        Tr = Tr * U / j
        Fc = Fc * j / U
        Sm = Sm + Tr
    R = D * (1 + 1 / (c - U) / (1 + (c - U) / U * Fc * Sm))

    N = X * R

    Qn = N - U
    Qtt = R - D

    
    Pn = np.zeros(pMax + 1)
    for i in range(0, pMax + 1):
        if i==0:
            Pn[i] = 1.0 / (Sm + Tr * Uave / (1.0 - Uave))
        elif i < c:
            Pn[i] = Pn[i-1] * U / i
        else:
            Pn[i] = Pn[i - 1] * Uave
    
    return {'type': 'M/M/c', 'X':X, 'R':R, 'N':N, 'U':U, 'Uave': Uave, 'Qn':Qn, 'Qtt': Qtt, 'Pn': Pn}
