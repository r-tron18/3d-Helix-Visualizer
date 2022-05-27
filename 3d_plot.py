from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import json
import numpy as np
import matplotlib.lines as mlines


fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')


X,Y,Z=[],[],[]
ax.scatter([0], [0], [0], c='r', marker='*')
with open('point.txt') as fil:
    data = fil.read()
    data = data.split('\n')
    for line in data:
        if not line:
            break
        x,y,z = line.split(',')
        x,y,z = float(x),float(y),float(z)
        X.append(x)
        Y.append(y)
        Z.append(z)

#ax.set_yticks(np.arange(-4,4,1))
#ax.set_xticks(np.arange(-15,-10,1))
#ax.set_zticks(np.arange(3,6,1))

ax.scatter(X[0], Y[0], Z[0], c='b', marker='o')
ax.scatter(X[1:-1], Y[1:-1], Z[1:-1], c='r', marker='o')
ax.scatter(X[-1], Y[-1], Z[-1], c='g', marker='o')

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()
