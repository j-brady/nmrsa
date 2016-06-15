import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d



x = np.zeros((6,6))
y = np.zeros((6,6))
z = np.zeros((6,6))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(x, y, z, rstride=10, cstride=10)
plt.show()
