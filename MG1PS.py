import numpy as np
from MM1 import MM1Solve


def MG1PSSolve(D, l, pMax = 10):
    res = MM1Solve(D, l, pMax)
    
    # removes percentiles
    del res['Rperc']
    res['type'] = 'M/G/1/PS'

    return res
