import numpy as np
from MMinf import MMinfSolve


def MGinfSolve(D, l, pMax = 10):
    res = MMinfSolve(D, l, pMax)
    res['type'] = 'M/G/inf'
    
    return res