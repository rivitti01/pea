import numpy as np
import matplotlib.pyplot as plt
from OPENmcss import OPENmcssSolve
from OPENmcms import OPENmcmsSolve


D = np.array([[1.0, 2.0], [1.8, 0.5], [0.25, 0.25]])
c = [1.0, 1.0, 1.0]


beta = np.zeros(101)
R = np.zeros((101, 4))

lT = 0.4

for n in range(0,101):
    beta[n] = float(n) / 100.0
    l = [lT * beta[n], lT * (1 - beta[n])]
    res = OPENmcssSolve(3, 2, D, c, l)
    R[n,0] = res['Rk'][0]
    R[n,1] = res['Rk'][1]
    R[n,2] = res['Rk'][2]
    R[n,3] = res['R']

plt.plot(beta, R)
plt.show()