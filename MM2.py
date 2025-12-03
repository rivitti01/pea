import numpy as np

def MM2Solve(D, l, pMax = 10):
    X = l
    U = l * D
    Uave = U / 2.0
    R = D / (1 - Uave * Uave)
    N = U / (1 - Uave * Uave)
    Qn = N - U
    Qtt = R - D
    
    Pn = np.zeros(pMax + 1)
    for i in range(0, pMax + 1):
        if i==0:
            Pn[i] = (1 - Uave) / (1 + Uave)
        elif i==1:
            Pn[i] = 2.0 * Pn[i-1] * Uave
        else:
            Pn[i] = Pn[i - 1] * Uave
    
    return {'type': 'M/M/2', 'X':X, 'R':R, 'N':N, 'U':U, 'Uave': Uave, 'Qn':Qn, 'Qtt': Qtt, 'Pn': Pn}
