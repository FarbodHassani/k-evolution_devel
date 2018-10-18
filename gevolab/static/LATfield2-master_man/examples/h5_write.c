#include "hdf5.h"
#include "iostream"

using namespace std;

int main ()
{

string H5FILE_NAME="SDS.h5",DATASETNAME="path";

int NX=5;                     /* dataset dimensions */
int NY=3;
int RANK=2;


    hid_t       file, dataset;         /* file and dataset handles */
    hid_t       datatype, dataspace,plist_id;   /* handles */
    hsize_t     dimsf[2];              /* dataset dimensions */
    herr_t      status;
    double         data[NX][NY];          /* data to write */
    int         i, j;

    /*
     * Data  and output buffer initialization.
     */
for(j = 0; j < NX; j++)
	for(i = 0; i < NY; i++)
	    data[j][i] = (double)(i + j);

    /*
     * Create a new file using H5F_ACC_TRUNC access,
     * default file creation properties, and default file
     * access properties.
     */

plist_id = H5Pcreate(H5P_FILE_ACCESS);

    file = H5Fcreate(H5FILE_NAME.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
H5Pclose(plist_id);
    /*
     * Describe the size of the array and create the data space for fixed
     * size dataset.
     */
    dimsf[0] = NX;
    dimsf[1] = NY;
    dataspace = H5Screate_simple(RANK, dimsf, NULL);

    /*
     * Define datatype for the data in the file.
     * We will store little endian INT numbers.
     */
    datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    status = H5Tset_order(datatype, H5T_ORDER_LE);

    /*
     * Create a new dataset within the file using defined dataspace and
     * datatype and default dataset creation properties.
     */
    dataset = H5Dcreate2(file, DATASETNAME.c_str(), datatype, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Write the data to the dataset using default transfer properties.
     */
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    /*
     * Close/release resources.
     */
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Fclose(file);

    return 0;
}

