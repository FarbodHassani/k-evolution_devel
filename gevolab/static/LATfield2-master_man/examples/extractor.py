import subprocess
import numpy as np
import h5py

import matplotlib.pylab as plt
import mpl_toolkits.mplot3d.axes3d as p3

def extractor(exe,file_name,id_add,outname,function='None'):
    f_snap0 = h5py.File(file_name+'.h5', 'r')
    snap0 = f_snap0['data']

    IDs = np.loadtxt(id_add)
    mask = np.in1d(snap0['ID'],IDs)

    subprocess.call(['mpirun', '-np', '4', exe, '-s', file_name, '-i', str(IDs.shape[0])])

    file_name = 'pre-out.h5'
    f1 = h5py.File(file_name, 'r+')
#    f1['data']['positionX'] =snap0[mask]

    if function==None:
        f1['data']['positionX'] =snap0[mask]['positionX']
        f1['data']['positionY'] =snap0[mask]['positionY']
        f1['data']['positionZ'] =snap0[mask]['positionZ']
    else:
        f1['data']['positionX'] =function(snap0[mask]['positionX'])
#        f1['data']['positionX'] =snap0[mask]['positionX']
        f1['data']['positionY'] =function(snap0[mask]['positionY'])
#        f1['data']['positionY'] =snap0[mask]['positionY']
        f1['data']['positionZ'] =function(snap0[mask]['positionZ'])
#        f1['data']['positionZ'] =snap0[mask]['positionZ']

        n_out = np.any((f1['data']['positionX']<=0,\
        f1['data']['positionX']>=1,\
        f1['data']['positionY']<=0,\
        f1['data']['positionY']>=1,\
        f1['data']['positionZ']<=0,\
        f1['data']['positionZ']>=1), axis=0).sum()
        if n_out!=0 :
            print 'Warning!'
            print '%d particles went out by transformation!' %n_out

    f1['data']['velocityX'] =snap0[mask]['velocityX']
    f1['data']['velocityY'] =snap0[mask]['velocityY']
    f1['data']['velocityZ'] =snap0[mask]['velocityZ']
    f1['data']['ID'] =snap0[mask]['ID']        

    poc_0 = np.logical_and(f1['data']['positionY']<0.5,f1['data']['positionZ']<0.5)
    poc_1 = np.logical_and(f1['data']['positionY']<0.5,f1['data']['positionZ']>0.5)
    poc_2 = np.logical_and(f1['data']['positionY']>0.5,f1['data']['positionZ']<0.5)
    poc_3 = np.logical_and(f1['data']['positionY']>0.5,f1['data']['positionZ']>0.5)

    n1 = poc_0.sum()
    n2 = poc_1.sum()
    n3 = poc_2.sum()
    n4 = poc_3.sum()

    intype = type(f1['numParts'][0])
    n1 = np.array([n1], dtype=intype)[0]
    n2 = np.array([n2], dtype=intype)[0]
    n3 = np.array([n3], dtype=intype)[0]
    n4 = np.array([n4], dtype=intype)[0]
    n_tot = n1+n2+n3+n4

    proc_list = np.zeros(n_tot)
    proc_list[poc_0] = 0
    proc_list[poc_1] = 1
    proc_list[poc_2] = 2
    proc_list[poc_3] = 3

    sorted_proc_list = proc_list.argsort()
    f1['data']['ID'] = f1['data']['ID'][sorted_proc_list]
    f1['data']['positionX'] = f1['data']['positionX'][sorted_proc_list]
    f1['data']['positionY'] = f1['data']['positionY'][sorted_proc_list]
    f1['data']['positionZ'] = f1['data']['positionZ'][sorted_proc_list]
    f1['data']['velocityX'] = f1['data']['velocityX'][sorted_proc_list]
    f1['data']['velocityY'] = f1['data']['velocityY'][sorted_proc_list]
    f1['data']['velocityZ'] = f1['data']['velocityZ'][sorted_proc_list]

    f1['numParts'][0] = n1
    f1['numParts'][1] = n2
    f1['numParts'][2] = n3
    f1['numParts'][3] = n4

    f1.close()

    subprocess.call(['mv', 'pre-out.h5',outname+'.h5'])

    print('Done!')

# Attaching 3D axis to the figure
def plotter(sn,outname):
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    # NOTE: Can't pass empty arrays into 3d version of plot()
    ax.plot(sn['positionX'], sn['positionY'], sn['positionZ'], 'ko')

    # Setting the axes properties
    ax.set_xlim3d([0.0, 1.0])
    ax.set_xlabel('X')

    ax.set_ylim3d([0.0, 1.0])
    ax.set_ylabel('Y')

    ax.set_zlim3d([0.0, 1.0])
    ax.set_zlabel('Z')

    ax.set_title('')

#    plt.savefig(outname+'.jpg')
    plt.show()

file_name = '/home/gf/work/forsat/geneva_works/gevolution/code/LATfield2-master/examples/files/lcdm_snap000_cdm'
outname = '/home/gf/work/forsat/geneva_works/gevolution/code/LATfield2-master/examples/files/ic'
id_add = '/home/gf/work/forsat/geneva_works/gevolution/code/LATfield2-master/examples/files/ids'
exe = '/home/gf/work/forsat/geneva_works/gevolution/code/LATfield2-master/examples/pre-ext'

def f(x):
	x_mean = np.mean(x)
	xmin = np.min(x)
	xmax = np.max(x)
	if xmax-xmin>0.8:
		x = x+0.5
	x = np.mod(x, 1.)

	xmin = np.min(x)
	xmax = np.max(x)
	x = 1.7*(x-xmin)
	x_mean = np.mean(x)

	return x-x_mean+0.5

#f=None

extractor(exe,file_name,id_add,outname,f)

#outname = '/home/gf/work/forsat/geneva_works/gevolution/code/LATfield2-master/examples/files/out2'
f1 = h5py.File(outname+'.h5', 'r')
#print f1['numParts'][0],f1['numParts'][1],f1['numParts'][2],f1['numParts'][3]
#print f1['numParts'][0]+f1['numParts'][1]+f1['numParts'][2]+f1['numParts'][3],

ot = np.array(list(f1['data']))
plotter(ot,'salam')

f1.close()


