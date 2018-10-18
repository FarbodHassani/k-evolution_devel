#! /home/gf/pakages/miniconda2/bin/python2.7
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation


data = np.loadtxt('path')[:,2:5].T

def update_lines(num, data, line):
    line.set_data(data[0:2, num-1:num])
    line.set_3d_properties(data[2, num-1:num])
    return line

# Attaching 3D axis to the figure
fig = plt.figure()
ax = p3.Axes3D(fig)

# NOTE: Can't pass empty arrays into 3d version of plot()
lines = ax.plot(data[0, 0:1], data[1, 0:1], data[2, 0:1], 'k-o')[0]

# Setting the axes properties
ax.set_xlim3d([0.0, 1.0])
ax.set_xlabel('X')

ax.set_ylim3d([0.0, 1.0])
ax.set_ylabel('Y')

ax.set_zlim3d([0.0, 1.0])
ax.set_zlabel('Z')

ax.set_title('')

# Creating the Animation object
anim = animation.FuncAnimation(fig, update_lines, 55, fargs=(data, lines),interval=50, blit=False)

anim.save('part_path.gif', writer='imagemagick')
#plt.show()

