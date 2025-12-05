import numpy as np
import matplotlib.pyplot as plt
from MVAmcss import MVAmcssSolve
from MVAmcms import MVAmcmsSolve


D = np.array([[1.0, 2.0], [1.8, 0.5], [0.25, 0.25]])
c = [1.0, 1.0, 1.0]
Z = [10.0, 15.0]

beta = np.zeros(101)
R = np.zeros((101, 3))

N = 100

for n in range(0,N+1):
    beta[n] = float(n) / N
    nc = [n, N-n]
    res = MVAmcssSolve(3, 2, D, c, nc, Z)
    R[n,0] = res['Rc'][0]
    R[n,1] = res['Rc'][1]
    R[n,2] = (res['Rc'][0] * res['Xc'][0] + res['Rc'][1] * res['Xc'][1]) / (res['Xc'][0] + res['Xc'][1])
plt.plot(beta, R)
plt.show()