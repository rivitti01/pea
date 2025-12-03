import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg


ref = 0

P = np.array([[0.0, 1.0, 0.0, 0.0],
              [0.2, 0.0, 0.3, 0.5],
              [0.0, 1.0, 0.0, 0.0],
              [0.0, 1.0, 0.0, 0.0]])

P[:, ref] = np.zeros(4)

Q = np.eye(4) - P
l = np.zeros(4)
l[ref] = 1.0

v = linalg.solve(Q.T, l)

print("Visits: ", v)