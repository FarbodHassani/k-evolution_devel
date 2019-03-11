
import numpy as np
import matplotlib.pylab as plt
def makeplot(file_1,file_2,color_fh,label_1,label_2):
    # ColorsI = ["red","blueviolet","olive","darkblue"]
    # ColorsII = ['darkred','purple','green','blue']
    plt.figure(figsize=(13,8))
    ax = plt.gca()
    ax.tick_params(axis = 'both', which = 'major', labelsize = 24)
    ax.tick_params(axis = 'both', which = 'minor', labelsize = 18)
    plt.loglog(file_1[:,0],file_1[:,1],color=color_fh[0],linestyle='solid',lw=2.5,label=label_1)
    plt.loglog(file_2[:,0], file_2[:,1],color=color_fh[1],linestyle='dashed',lw=2.5, label=label_2 )
    # plt.loglog(matter_int_2[:,0],matter_int_2[:,1],color="green",linestyle='solid',lw=2.5,label="interpolated")
    plt.legend()
    plt.show()

def makediffplot(file_1,file_2,color_fh,label_1):
    # ColorsI = ["red","blueviolet","olive","darkblue"]
    # ColorsII = ['darkred','purple','green','blue']
    plt.figure(figsize=(13,5))
    ax = plt.gca()
    ax.tick_params(axis = 'both', which = 'major', labelsize = 24)
    ax.tick_params(axis = 'both', which = 'minor', labelsize = 18)
    plt.loglog(file_1[:,0],np.abs(file_1[:,1]-file_2[:,1])/file_2[:,1],color=color_fh[0],linestyle='solid',lw=2.5,label=label_1)
    # plt.axhline(y=0.05,color="black",lw=1.5,label=r"$5\%$ error")
    # plt.ylabel(r"$\Delta\mathcal{P}_{\delta/\mathcal{P}_{\delta }$",fontsize=25)

    # plt.loglog(file_2[:,0], file_2[:,1],color=color_fh[1],linestyle='dashed',lw=2.5, label=label_2 )
    # plt.loglog(matter_int_2[:,0],matter_int_2[:,1],color="green",linestyle='solid',lw=2.5,label="interpolated")
    plt.legend()
    plt.show()
