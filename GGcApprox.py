import numpy as np
from MMc import MMcSolve


def GGcApproxSolve(D, cv, c, l, ca):
    resMMc = MMcSolve(D, c, l)
    
    X = l
    U = l * D
    Qtt = (cv * cv + ca * ca) / 2.0 * resMMc['Qtt']
    Qn = Qtt * l
    R = D + Qtt
    N = U + Qn

    return {'type': 'G/G/c (Approx.)', 'X':X, 'R':R, 'N':N, 'U':U, 'Qn':Qn, 'Qtt': Qtt}
