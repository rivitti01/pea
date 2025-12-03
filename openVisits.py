import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg



P = np.array([[0.0, 0.5, 0.3],
              [1.0, 0.0, 0.0],
              [1.0, 0.0, 0.0]])

Q = np.eye(3) - P
l = np.array([1.0, 0.0, 0.0])

v = linalg.solve(Q.T, l)

print("Visits: ", v)