import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
def extrapolation_z0 (z_ini, z_ini_index, order, Gev_lin_delta_m_zall ):
# order=3
    Gev_deltam_z0_extpl=np.zeros((np.shape(Gev_lin_delta_m_zall[18][:,0])[0],2))
    for i in range(np.shape(Gev_lin_delta_m_zall[z_ini_index[0]][:,1])[0]):
        y_data=[Gev_lin_delta_m_zall[z_ini_index[0]][i,1],Gev_lin_delta_m_zall[z_ini_index[1]][i,1],Gev_lin_delta_m_zall[z_ini_index[2]][i,1]]
        s = InterpolatedUnivariateSpline(z_ini, y_data, k=order)
        Gev_deltam_z0_extpl[i,0] = Gev_lin_delta_m_zall[z_ini_index[0]][i,0]
        Gev_deltam_z0_extpl[i,1] = s(0.0)

    return Gev_deltam_z0_extpl
