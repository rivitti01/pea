import numpy as np

def MMinfSolve(D, l, pMax = 10):
    X = l
    U = l * D
    
    R = D
    N = U
    Qn = 0
    Qtt = 0
    
    Pn = np.zeros(pMax + 1)
    for i in range(0, pMax + 1):
        if i==0:
            Pn[i] = np.exp(-U)
        else:
            Pn[i] = Pn[i - 1] * U / i
    

    return {'type': 'M/M/inf', 'X':X, 'R':R, 'N':N, 'U':U, 'Qn':Qn, 'Qtt': Qtt, 'Pn': Pn}
