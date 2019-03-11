import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
def interp_files(Gev_lin_deltakess_z_all, z_index_interp,zrange, kstep, kini, k_fin, cosmic_var_skip,nyqvist_skip_l,nyqvist_skip_m,nyqvist_skip_s):
    k=kini;
    # zrange=30;
    k_range_full=[]
    while k<k_fin:
        k=k+kstep;
        k_range_full.append(round(k, 3))

    kini_larg = round(0.002*cosmic_var_skip,3);
    k_fin_larg = round(1.340/nyqvist_skip_l,3) #2.320;

    k=kini_larg;
    k_range_l=[]
    while k<k_fin_larg:
        k=k+kstep;
        k_range_l.append(round(k, 3))
        ###########
    # Interpolation of power and count at ks
    ###########

    #Redshifts are =[zlist.index(2.0),zlist.index(1.0),zlist.index(0.5),zlist.index(0.0)];
    # counter = 0
    Interp_delta_kess_kev_z_all_Large = []
    for i in z_index_interp:
    # for i in [8,19]:
        Interp_delta_kess_kev_z=np.zeros((np.shape(k_range_full)[0],3))
        interp_deltakess_pow_l=InterpolatedUnivariateSpline( Gev_lin_deltakess_z_all[i][:,0],Gev_lin_deltakess_z_all[i][:,1],k=5)
        interp_deltakess_count_l=(InterpolatedUnivariateSpline( Gev_lin_deltakess_z_all[i][:,0],Gev_lin_deltakess_z_all[i][:,4],k=5))
        Interp_delta_kess_kev_z[:np.shape(k_range_l)[0],0]= k_range_l
        Interp_delta_kess_kev_z[:np.shape(k_range_l)[0],1]= interp_deltakess_pow_l(k_range_l)
        Interp_delta_kess_kev_z[:np.shape(k_range_l)[0],2]= interp_deltakess_count_l(k_range_l)
        Interp_delta_kess_kev_z_all_Large.append(Interp_delta_kess_kev_z);
    #     print("test:",Interp_delta_kess_kev_z[:np.shape(k_range_l)[0]])
    #     print("append list:",Interp_delta_kess_kev_z_all_Large)
    #     print("\n")
    #     counter = counter +1
    ########## mid ###
    kini_mid = round(0.013*cosmic_var_skip,3);
    k_fin_mid = round(9.500/nyqvist_skip_m,3);
    k=kini_mid;
    k_range_m=[]
    while k<k_fin_mid:
        k=k+kstep;
        k_range_m.append(round(k, 3))
    # # ###########
    # # # Interpolation of power and count at ks
    # # ###########

    Interp_delta_kess_kev_z_all_Mid = []

    # # # #Redshifts are =[zlist.index(2.0),zlist.index(1.0),zlist.index(0.5),zlist.index(0.0)];
    for i in z_index_interp:
        Interp_delta_kess_kev_z=np.zeros((np.shape(k_range_full)[0],3))

        interp_deltakess_pow_m=InterpolatedUnivariateSpline( Gev_lin_deltakess_z_all[i+zrange][:,0],Gev_lin_deltakess_z_all[i+zrange][:,1],k=5)
        interp_deltakess_count_m=(InterpolatedUnivariateSpline( Gev_lin_deltakess_z_all[i+zrange][:,0],Gev_lin_deltakess_z_all[i+zrange][:,4],k=5))
        Interp_delta_kess_kev_z[:np.shape(k_range_m)[0],0]= k_range_m
        Interp_delta_kess_kev_z[:np.shape(k_range_m)[0],1]= interp_deltakess_pow_m(k_range_m)
        Interp_delta_kess_kev_z[:np.shape(k_range_m)[0],2]= interp_deltakess_count_m(k_range_m)
        Interp_delta_kess_kev_z_all_Mid.append(Interp_delta_kess_kev_z)

    kini_small = round(0.030*cosmic_var_skip,3);
    k_fin_small = round(24.500/nyqvist_skip_s,3);
    k=kini_small;
    k_range_s=[]
    while k<k_fin_small:
        k=k+kstep;
        k_range_s.append(round(k, 3))

    # # ###########
    # # # Interpolation of power and count at ks
    # # ###########
    Interp_delta_kess_kev_z_all_Small = []

    # #Redshifts are =[zlist.index(2.0),zlist.index(1.0),zlist.index(0.5),zlist.index(0.0)];
    for i in z_index_interp:
        Interp_delta_kess_kev_z=np.zeros((np.shape(k_range_full)[0],3))
    #     print(i,zlist_class[i])
        interp_deltakess_pow_s=InterpolatedUnivariateSpline( Gev_lin_deltakess_z_all[i+zrange*2][:,0],Gev_lin_deltakess_z_all[i+zrange*2][:,1],k=5)
        interp_deltakess_count_s=(InterpolatedUnivariateSpline( Gev_lin_deltakess_z_all[i+zrange*2][:,0],Gev_lin_deltakess_z_all[i+zrange*2][:,4],k=5))
        Interp_delta_kess_kev_z[:np.shape(k_range_s)[0],0]= k_range_s
        Interp_delta_kess_kev_z[:np.shape(k_range_s)[0],1]= interp_deltakess_pow_s(k_range_s)
        Interp_delta_kess_kev_z[:np.shape(k_range_s)[0],2]= interp_deltakess_count_s(k_range_s)
        Interp_delta_kess_kev_z_all_Small.append(Interp_delta_kess_kev_z)

    return Interp_delta_kess_kev_z_all_Large, Interp_delta_kess_kev_z_all_Mid,Interp_delta_kess_kev_z_all_Small
