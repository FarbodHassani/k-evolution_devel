import numpy as np
def wighted_power(arr_1, arr_2, arr_3, number_binning, kmin, kmax,nyqvist_skip,cosmicvar_skip):
    du  = (np.log10(kmax)-np.log10(kmin) )/number_binning;
    weight_pow=np.zeros((number_binning,3));
    k = kmin;
    for i in range(number_binning):
        k = k * 10** (du);
        dk = k * (10** (du) -1.0);
        for j in range (np.int(np.shape(arr_1[:,0])[0])):
            if(arr_1[j,0]>k-dk/2. and arr_1[j,0]<k+dk/2. and arr_1[j,0]<(arr_1[:,0].max())/nyqvist_skip and j> 1*cosmicvar_skip ):
                weight_pow[i,1] = weight_pow[i,1] + arr_1[j,1] * arr_1[j,4];
                weight_pow[i,2] = weight_pow[i,2] + arr_1[j,4];

            if(arr_2[j,0]>k-dk/2. and arr_2[j,0]<k+dk/2. and arr_2[j,0]<(arr_2[:,0].max())/nyqvist_skip and j> 1*cosmicvar_skip):
                weight_pow[i,1] = weight_pow[i,1] + arr_2[j,1] * arr_2[j,4];
                weight_pow[i,2] = weight_pow[i,2] + arr_2[j,4];

            if(arr_3[j,0]>k-dk/2. and arr_3[j,0]<k+dk/2. and arr_3[j,0]<(arr_3[:,0].max())/nyqvist_skip and j> 1*cosmicvar_skip):
                weight_pow[i,1] = weight_pow[i,1] + arr_3[j,1] * arr_3[j,4];
                weight_pow[i,2] = weight_pow[i,2] + arr_3[j,4];

        weight_pow[i,0] =k;
        if (weight_pow[i,1] !=0):
            weight_pow[i,1] =weight_pow[i,1]/weight_pow[i,2];
        else:
            weight_pow[i,1] =weight_pow[i,1];

    return weight_pow
