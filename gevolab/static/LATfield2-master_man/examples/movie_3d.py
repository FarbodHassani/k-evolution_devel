#! /home/gf/pakages/miniconda2/bin/python2.7
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from moviepy.editor import VideoClip

import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

clr=['r','b','g','k','orange']

data = np.loadtxt('path')

n_part = int(np.max(data[:,0])+1)

parts = [data[data[:,0]==i][:,2:5] for i in range(n_part)]
p_dot = [None for i in range(n_part)]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

ax.set_xlim(0,1)
ax.set_ylim(0,1)
ax.set_zlim(0,1)
#plt.show()

outf = 'test.avi'
rate = 10

cmdstring = ('ffmpeg',
             '-r', '%d' % rate,
             '-f','image2pipe',
             '-vcodec', 'png',
             '-i', 'pipe:', outf
             )
p = subprocess.Popen(cmdstring, stdin=subprocess.PIPE)

frames = parts[0].shape[0]

for t in range(frames):
	for i in range(n_part):
		path = parts[i]
		p_dot[i] = ax.scatter(path[t,0], path[t,1], path[t,2], c=clr[i], marker='o')

	ax.view_init(60, 20)

	plt.savefig(p.stdin, format='png')

	for i in range(n_part):
		p_dot[i].remove()

p.stdin.close()



#import numpy as np
#import matplotlib.pyplot as plt
#from sklearn import svm # sklearn = scikit-learn
#from sklearn.datasets import make_moons
#from moviepy.editor import VideoClip
#from moviepy.video.io.bindings import mplfig_to_npimage

#X, Y = make_moons(50, noise=0.1, random_state=2) # semi-random data

#fig, ax = plt.subplots(1, figsize=(4, 4), facecolor=(1,1,1))
#fig.subplots_adjust(left=0, right=1, bottom=0)
#xx, yy = np.meshgrid(np.linspace(-2,3,500), np.linspace(-1,2,500))

#def make_frame(t):
#    ax.clear()
#    ax.axis('off')
#    ax.set_title("SVC classification", fontsize=16)

#    classifier = svm.SVC(gamma=2, C=1)
#    # the varying weights make the points appear one after the other
#    weights = np.minimum(1, np.maximum(0, t**2+10-np.arange(50)))
#    classifier.fit(X, Y, sample_weight=weights)
#    Z = classifier.decision_function(np.c_[xx.ravel(), yy.ravel()])
#    Z = Z.reshape(xx.shape)
#    ax.contourf(xx, yy, Z, cmap=plt.cm.bone, alpha=0.8,
#                vmin=-2.5, vmax=2.5, levels=np.linspace(-2,2,20))
#    ax.scatter(X[:,0], X[:,1], c=Y, s=50*weights, cmap=plt.cm.bone)

#    return mplfig_to_npimage(fig)

#animation = VideoClip(make_frame, duration = 7)
#animation.write_gif("svm.gif", fps=15)

