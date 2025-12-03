import numpy as np

def MM1KSolve(D, K, l, pMax = -1):
    if (pMax == -1) or (pMax > K):
        pMax = K

    rho = l * D      
    
    if rho == 1:
        U = K / (K + 1)
        N = K / 2.0
        pB = 1 / (K + 1)
    else:
        U = rho * (1.0 - np.pow(rho, K)) / (1.0 - np.pow(rho, K+1))
        N = rho / (1 - rho) - (K + 1) * np.pow(rho, K + 1) / (1.0 - np.pow(rho, K + 1))
        pB = np.pow(rho, K) * (1.0 - rho) / (1.0 - np.pow(rho, K+1))

    Dr = pB  * l
    X = l - Dr
    R = N / X
    Qn = N - U
    Qtt = R - D
    
    Pn = np.zeros(pMax + 1)
    for i in range(0, pMax + 1):
        if rho == 1:
            Pn[i] = 1 / (K + 1)
        else:
            if i==0:
                Pn[i] = (1 - rho) / (1.0 - np.pow(rho, K+1))
            else:
                Pn[i] = Pn[i - 1] * rho
    
    return {'type': 'M/M/1/K', 'X':X, 'R':R, 'N':N, 'U':U, 'Qn':Qn, 'Qtt': Qtt, 'Pn': Pn,
            'pB': pB, 'Dr': Dr}
