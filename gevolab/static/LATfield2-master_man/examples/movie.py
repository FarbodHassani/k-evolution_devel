#! /home/gf/pakages/miniconda2/bin/python2.7

"""
===========
MovieWriter
===========

This example uses a MovieWriter directly to grab individual frames and write
them to a file. This avoids any event loop integration, but has the advantage
of working with even the Agg backend. This is not recommended for use in an
interactive setting.

"""
# -*- noplot -*-

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Movie Test', artist='Matplotlib',
                comment='Movie support!')
writer = FFMpegWriter(fps=15, metadata=metadata)

fig = plt.figure()
l, = plt.plot([], [], 'k-o')

plt.xlim(0, 1)
plt.ylim(0, 1)

x0, y0 = 0, 0

data = np.loadtxt('path')
#data = data[data[:,0].argsort()]

#p = data[data[:,0]==0][:,2:5]
p = data[data[:,0]==1][:,2:5]

x = p[:,0]
y = p[:,1]

with writer.saving(fig, "writer_test.mp4", data.shape[0]):
    for i in range(p.shape[0]):
        x0 = x[i]
        y0 = y[i]
        l.set_data(x0, y0)
        writer.grab_frame()
