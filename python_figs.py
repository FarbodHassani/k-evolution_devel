import sys
import os
sys.path.insert(0,'../')

#from readgadget import *
#from pygadgetreader import *
#from pylab import *
import numpy as np

import h5py
import matplotlib.pylab as plt
#%matplotlib inline
##%matplotlib notebook
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import matplotlib.pylab as plt
import matplotlib.animation as animation
import subprocess
im=listofzeros = [0] * 40
#,'cdm'
snaps=['pi_k']
for i in range(29):
    f = h5py.File('output/pi_k_'+str(i)+'.h5', "r")
    s = np.array(list(f['data']))
    im[i]=plt.imshow(s[:,:,1])
    plt.savefig('Plots/Field/kess'+str(i)+'.jpg')
    plt.close()
    del f
im=listofzeros = [0] * 40
#,'cdm'
snaps=['pi_k']
for i in range(29):
    f = h5py.File('output/pi_v_'+str(i)+'.h5', "r")
    s = np.array(list(f['data']))
    im[i]=plt.imshow(s[:,:,1])
    plt.savefig('Plots/Velocity/vkess'+str(i)+'.jpg')
    plt.close()
    del f
import os
os.system("ffmpeg -f image2 -r 1 -i ./Plots/Field/kess%01d.jpg -vcodec mpeg4 -y ./Plots/Movie_kess.mp4")
import os
os.system("ffmpeg -f image2 -r 1 -i ./Plots/Velocity/vkess%01d.jpg -vcodec mpeg4 -y ./Plots/Movie_vkess.mp4")
