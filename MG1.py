import numpy as np

def MG1Solve(D, cv, l):
    X = l
    U = l * D
    Qn = U * U * (1 + cv * cv) / (2.0 * (1.0 - U))
    Qtt = U * D * (1 + cv * cv) / (2.0 * (1.0 - U))
    R = D + Qtt
    N = U + Qn
    
    return {'type': 'M/G/1', 'X':X, 'R':R, 'N':N, 'U':U, 'Qn':Qn, 'Qtt': Qtt}
