import numpy as np

def MM1Solve(D, l, pMax = 10, perc = [25.0, 50.0, 75.0]):
    X = l
    U = l * D
    R = D / (1 - U)
    N = U / (1 - U)
    Qn = N - U
    Qtt = R - D
    
    Pn = np.zeros(pMax + 1)
    for i in range(0, pMax + 1):
        if i==0:
            Pn[i] = 1 - U
        else:
            Pn[i] = Pn[i - 1] * U
    
    Rperc = np.zeros(len(perc))
    i = 0
    for p in perc:
        Rperc[i] = - np.log(1.0 - p / 100.0) * R
        i = i + 1

    return {'type': 'M/M/1', 'X':X, 'R':R, 'N':N, 'U':U, 'Qn':Qn, 'Qtt': Qtt, 'Pn': Pn,
            'Rperc': Rperc}
