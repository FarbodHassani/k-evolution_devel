#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <array>
#include "chealpix.h"
#include "metadata.hpp"
#include "parser.hpp"
#include "background.hpp"

// g++ lcmap.cpp -o lcmap -std=c++11 -O3 -DCOLORTERMINAL -I[Healpix Directory]/include -L[Healpix Directory]/lib -lchealpix -lcfitsio -lgsl -lgslcblas

using namespace std;

struct healpix_header
{
	uint32_t Nside;
	uint32_t Npix;
	uint32_t precision;
	uint32_t Ngrid;
	double direction[3];
	double distance;
	double boxsize;
	uint32_t Nside_ring;
	char fill[256 - 5 * 4 - 5 * 8]; /* fills to 256 Bytes */
	uint32_t header_blocksize;
	uint32_t data_blocksize; /* to get rid of two fread commands per header/data pair */
};

struct background_data
{
	int cycle;
	double tau;
	double a;
};

struct metric_data
{
	healpix_header hdr;
	float * pixel;
};

struct metric_container
{
	std::map<int,metric_data> healpix_data;
	background_data cinfo;
	char name[4];
	char * dir;
	char * basename;
	void init(background_data & c, char * d, char * b, const char * n)
	{
		cinfo = c;
		strcpy(name, n);
		dir = d;
		basename = b;
	}
	void truncate(double min_dist)
	{
		std::map<int,metric_data>::iterator it;

		for (it = healpix_data.begin(); it != healpix_data.end(); it++)
		{
			if (it->second.hdr.distance >= min_dist) break;
			else free(it->second.pixel);
		}

		if (it != healpix_data.begin()) healpix_data.erase(healpix_data.begin(), it);
	}
	void clear()
	{
		for (std::map<int,metric_data>::iterator it = healpix_data.begin(); it != healpix_data.end(); it++)
			free(it->second.pixel);
		healpix_data.clear();
	}
};

bool pixgrad(float * pixel, const int64_t Nside, const int64_t Npix, int64_t ipix, double * result);

int loadHealpixData(metric_container * field, double min_dist, double max_dist, char * lightconeparam = NULL);


float lin_int(float a, float b, float f)
{
	return a + (f*(b-a));
}


double distance_calc(const double a, const cosmology cosmo)
{	
	return (C_SPEED_OF_LIGHT/cosmo.h)  * (1. / sqrt(((cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) * pow(a,-3)) + (cosmo.Omega_Lambda) + (cosmo.Omega_rad * pow(a,-4)) + ((cosmo.Omega_fld * exp(3. * cosmo.wa_fld * (a - 1.)) / pow(a, 1. + 3. * (cosmo.w0_fld + cosmo.wa_fld)))/(a * a))));
}


void rungekutta4distancez(double &dist, double z, const cosmology cosmo, const double h)
{
	double k1a, k2a, k4a;

	//a=1/(1+z);

	k1a = distance_calc(1/(1+z), cosmo); 
	k2a = distance_calc(1/(1+z+(h*0.5)), cosmo);
	k4a = distance_calc(1/(1+z+h), cosmo);

	dist += h * (k1a + 4. * k2a + k4a) / 6.;
}


int main(int argc, char **argv)
{

	char * settingsfile = NULL;
	char * redshiftparam = NULL;
	char outputfile[1024];
	char * distanceparam = NULL;
	char * lightconeparam = NULL;

	parameter * params = NULL;
	metadata sim;
	cosmology cosmo;
	icsettings ic;
	int numparam = 0;
	int usedparams = 0;

	int numoutputs=0;

	char filename0[1024];
	char filename1[1024];

	char coordsys = 'G';

	if (argc < 2)
	{
		cout << COLORTEXT_WHITE << " LCARS tools: lcmap" << COLORTEXT_RESET << endl;
		cout << " extracts local and integrated potential terms from metric light-cone output" << endl << endl;
		
		cout << " List of command-line options:" << endl;
		cout << " -s <filename>       : gevolution settings file of the simulation (mandatory)" << endl;
		cout << " -d <distance>[,...] : conformal distance to source field [Mpc/h] (alternative 1)" << endl;
		cout << " -z <redshift>[,...] : redshift of source field (alternative 2)" << endl;
		cout << " -l <ID1>[,<ID2>,...]: IDs of light cones to be included (optional, must refer" << endl;
		cout << "                       to the same observation event)" << endl;
		cout << " The output will be written to HEALPix FITS files that follow the naming conventions" << endl;
		cout << " specified in the settings file." << endl;

		return 0;
	}


	for (int i = 1 ; i < argc ; i++ )
	{
		if ( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) {
			case 'z':
				redshiftparam = argv[++i]; //redshifts for calculation 
				break;
			case 's':
				settingsfile = argv[++i]; //settings file name
				break;
			case 'd':
				distanceparam= argv[++i]; //distances for calculation
				break;
			case 'l':
				lightconeparam = argv[++i];
		}
	}

	if (settingsfile == NULL)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no settings file specified!" << endl;
		return -1;
	}

	cout << COLORTEXT_WHITE << " LCARS tools: lcmap" << endl << endl << " opening settings file of simulation: " << settingsfile << endl << " parser output:" << COLORTEXT_RESET << endl << endl;

	numparam = loadParameterFile(settingsfile, params);
	usedparams = parseMetadata(params, numparam, sim, cosmo, ic);

	cout << endl << " file contains " << numparam << " parameters, " << usedparams << " of which could be parsed." << endl << endl;

	cout << " number of lightcones: " << sim.num_lightcone << endl << endl;

	for (int i = 0; i < sim.num_lightcone; i++)
	{
		cout << " lightcone " << i << " parameters:" << endl << "  vertex = " << sim.lightcone[i].vertex[0] << ", " << sim.lightcone[i].vertex[1] << ", " << sim.lightcone[i].vertex[2] << endl;
		cout << "  redshift of observation = " << sim.lightcone[i].z << endl;
		cout << "  direction = " << sim.lightcone[i].direction[0] << ", " << sim.lightcone[i].direction[1] << ", " << sim.lightcone[i].direction[2] << endl;
		cout << "  opening half-angle = " << ((sim.lightcone[i].opening > -1.) ? acos(sim.lightcone[i].opening) * 180. / M_PI : 180.) << " degrees" << endl;
		cout << "  distance interval = " << sim.lightcone[i].distance[0] << ", " << sim.lightcone[i].distance[1] << endl << endl;
	}

	char * token = NULL;

	if(redshiftparam != NULL && distanceparam != NULL)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": please specify EITHER redshift OR distance to source field, not both!" << endl;
		return 0;
	}
	else if(redshiftparam == NULL && distanceparam == NULL)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": please specify redshift (parameter -z) OR distance (parameter -d) to source field!" << endl;
		return 0;		
	}
	else if(redshiftparam == NULL && distanceparam != NULL)
	{
		sprintf(filename0, "%s", distanceparam);

		token = strtok(filename0, ",");

		while (token != NULL)
		{
			numoutputs++;
			token = strtok(NULL, ",");
		}		
	}
	else if(redshiftparam != NULL && distanceparam == NULL)
	{
		sprintf(filename0, "%s", redshiftparam);

		token = strtok(filename0, ",");

		while (token != NULL)
		{
			numoutputs++;
			token = strtok(NULL, ",");
		}		
	}


	double redshifts[numoutputs];
	double distances[numoutputs];

	int n = 0;

	double maxredshift = 190;


	if(redshiftparam == NULL && distanceparam != NULL)
	{

		sprintf(filename0, "%s", distanceparam);

		token = strtok(filename0, ",");

		while (token != NULL)
		{
			distances[n] = atof(token);
			if (distances[n] > 200000000)
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": distance of " << distances[n] << " is too large." << endl;
				return -1;
			}
			n++;
			token = strtok(NULL, ",");
		}
		
		qsort((void *) distances, numoutputs, sizeof(double), [](const void * a, const void * b)->int{ return (int) (*((double*) a) > *((double*) b)); });

		cout <<" distances (Mpc/h) chosen are: ";

		for (int i = 0; i < numoutputs; i++)
		{
			if(i!=0) cout << ", ";
			cout << distances[i];
			distances[i] /= sim.boxsize;
		}
		cout << "." << endl << endl;
	}
	else if(redshiftparam != NULL && distanceparam == NULL)
	{

		sprintf(filename0, "%s", redshiftparam);

		token = strtok(filename0, ",");

		while (token != NULL)
		{
			redshifts[n] = atof(token);
			if (redshifts[n] > maxredshift)
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": redshift of " << redshifts[n] << " is too large." << endl;
				return -1;
			}
			n++;
			token = strtok(NULL, ",");
		}
	
		
		qsort((void *) redshifts, numoutputs, sizeof(double), [](const void * a, const void * b)->int{ return (int) (*((double*) a) > *((double*) b)); });

		cout <<" redshifts chosen are:";

		for (int i = 0; i < numoutputs; i++)
		{
			distances[i] = 0;
			double h = redshifts[i]/1000;
			
			for (int j=1;  j< 1001; j++)
			{
				rungekutta4distancez(distances[i], 0+(j*h),cosmo, h);
			}
			if(i!=0) cout << ", ";
			cout << redshifts[i];
			distances[i] /= sim.boxsize;
		}
		cout << "." << endl << endl;
	}

	FILE * background_file;
	char * buffer;

	int numlines = -2;

	sprintf(filename0,"%s%s_background.dat",sim.output_path,sim.basename_generic);
	background_file = fopen(filename0, "r");

	if (background_file==NULL)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": background file " << filename0 << " cannot be opened!" << endl;
		return -1;
	}

	while(!feof(background_file))
	{
		if(fscanf(background_file, "%*[^\n]\n") != 0)
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": read error! " << endl;
			return 0;
		}
		numlines++;
	}

	rewind(background_file);

	if(fscanf(background_file, "%*[^\n]\n") != 0)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": read error! " << endl;
		return 0;
	}
	if(fscanf(background_file, "%*[^\n]\n") != 0)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": read error! " << endl;
		return 0;
	}

	background_data back[numlines];

	for (int i=0;i<numlines;i++)
	{
		if(fscanf(background_file, "%i %lf %lf %*f %*f %*f \n", &back[numlines-1-i].cycle, &back[numlines-1-i].tau, &back[numlines-1-i].a) != 3)
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": read error! " << endl;
			return 0;
		}
	}

	fclose(background_file);


	metric_container phi0;
	metric_container phi1;

	phi0.init(back[0], sim.output_path, sim.basename_lightcone, "phi");
	phi1.init(back[1], sim.output_path, sim.basename_lightcone, "phi");

	float * map_phi_final[numoutputs];
	float * map_isw_final = NULL;
	float * map_shapiro_final = NULL;
	float pot_obs = 0.;

	uint32_t Nside_final=2;
	int64_t Npix_final=48;

	int64_t q, ipix, jpix, ring;

	int step=0;
	int outcnt = 0;
	int cnt = 0;
	int thresh = 1;

	double dist = 0;
	double tauobs = particleHorizon(1, 1.5 * sim.boxsize * sim.boxsize / C_SPEED_OF_LIGHT / C_SPEED_OF_LIGHT, cosmo);
	double monopole;
	double v1[3], v2[3], grad_phi[5], grad_isw[5], grad_shapiro[5], dp[2];

	std::map<int,metric_data>::iterator it0;
	std::map<int,metric_data>::iterator it1;
 	
	cout << COLORTEXT_YELLOW << " processing light-cone output..." << COLORTEXT_RESET << endl << endl;


	while(outcnt < numoutputs)
	{
		cout << " interpolating cycle number " << COLORTEXT_CYAN << back[step].cycle << COLORTEXT_RESET << " and " << COLORTEXT_CYAN << back[step+1].cycle << COLORTEXT_RESET << "." << endl << " distance interval: " << (tauobs - back[step].tau) * sim.boxsize << " to " << (tauobs - back[step + 1].tau) * sim.boxsize << endl << endl;

		loadHealpixData(&phi0, tauobs - back[step].tau, tauobs - back[step + 1].tau, lightconeparam);
		loadHealpixData(&phi1, tauobs - back[step].tau, tauobs - back[step + 1].tau, lightconeparam);

		for (it0 = phi0.healpix_data.begin(); it0 != phi0.healpix_data.end() && it0->second.hdr.distance < tauobs - back[step].tau; it0++);
		for (it1 = phi1.healpix_data.begin(); it1 != phi1.healpix_data.end() && it1->second.hdr.distance < tauobs - back[step].tau; it1++);
		
		if (it0 == phi0.healpix_data.end() || it1 == phi0.healpix_data.end())
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": missing HEALPix data beyond " << dist * sim.boxsize << " Mpc/h!" << endl;
			return -1;
		}

		if (map_isw_final == NULL)
		{
			Nside_final = it0->second.hdr.Nside;
			Npix_final = it0->second.hdr.Npix;
			for (int m = 0; m < numoutputs; m++)
				map_phi_final[m] = (float *) malloc(Nside_final * Nside_final * 12 * sizeof(float));
			map_isw_final = (float *) malloc(Nside_final * Nside_final * 12 * sizeof(float));
			map_shapiro_final = (float *) malloc(Nside_final * Nside_final * 12 * sizeof(float));
			
			monopole = 0.;

#pragma omp parallel for reduction(+:monopole)
			for (int l = 0; l < it0->second.hdr.Npix; l++)
				monopole += lin_int(it0->second.pixel[l], it1->second.pixel[l], (it0->second.hdr.distance - (tauobs - back[step].tau))/((back[step].tau - back[step+1].tau)));
				
			pot_obs = monopole / it0->second.hdr.Npix;
			
#pragma omp parallel for
			for(int i = 0; i < (Nside_final * Nside_final * 12); i++)
			{
				for (int m = 0; m < numoutputs; m++)
					map_phi_final[m][i] = 0;
				map_isw_final[i] = 0;
				map_shapiro_final[i] = 0;
			}
		}

		for(; it0->second.hdr.distance < tauobs - back[step+1].tau && it1->second.hdr.distance < tauobs - back[step+1].tau && dist < distances[numoutputs-1]; it0++, it1++)
		{
			if (it0 == phi0.healpix_data.end() || it1 == phi0.healpix_data.end())
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": missing HEALPix data beyond " << dist * sim.boxsize << " Mpc/h!" << endl;
				return -1;
			}

			if(it0->second.hdr.Nside != Nside_final)
			{
				cout << COLORTEXT_CYAN << " map resolution change" << COLORTEXT_RESET;
				cout << ": interpolating from Nside of " << COLORTEXT_WHITE << Nside_final << COLORTEXT_RESET << " to " << COLORTEXT_WHITE <<  it0->second.hdr.Nside << COLORTEXT_RESET << "." << endl << endl;
					
				float * map_phi_temp[numoutputs-outcnt];
				float * map_isw_temp = (float *) malloc(it0->second.hdr.Nside * it0->second.hdr.Nside * 12 * sizeof(float));
				float * map_shapiro_temp = (float *) malloc(it0->second.hdr.Nside * it0->second.hdr.Nside * 12 * sizeof(float));
					
				for (int m = 0; m < numoutputs-outcnt; m++)
					 map_phi_temp[m] = (float *) malloc(it0->second.hdr.Nside * it0->second.hdr.Nside * 12 * sizeof(float));

#pragma omp parallel for private(q,v1,v2,grad_phi,grad_isw,grad_shapiro,dp,jpix,ring)					
				for(int l = 0; l <  it0->second.hdr.Nside *  it0->second.hdr.Nside * 12; l++)
				{
					pix2vec_ring64(it0->second.hdr.Nside, l, v1);
					vec2pix_ring64(Nside_final, v1, &q);
					
					if (pixgrad(map_isw_final, (int64_t) Nside_final, Npix_final, q, grad_isw) && pixgrad(map_shapiro_final, (int64_t) Nside_final, Npix_final, q, grad_shapiro))
					{
						pix2vec_ring64(Nside_final, q, v2);
							
						if (v1[2] > 2./3.)
						{
							dp[0] = sqrt((1.-v1[2]) * 3.) * Nside_final;
							ring = (1 + (int64_t) sqrt(1.5 + 2 * (double) l)) / 2;
							jpix = 2 * ((l - 2 * ring * (ring-1)) % ring) + 1;
							if (jpix > ring)
								jpix -= 2*ring;
									
							dp[1] = jpix * 0.5 / ring;
						}
						else if (v1[2] >= -2./3.)
						{
							dp[0] = (2. - 1.5 * v1[2]) * Nside_final;
							jpix = l - 2 * it0->second.hdr.Nside * (it0->second.hdr.Nside + 1);
							ring = it0->second.hdr.Nside + 1 + jpix / (4 * it0->second.hdr.Nside);
							if (jpix < 0)
							{
								ring = it0->second.hdr.Nside;
								jpix += 4 * it0->second.hdr.Nside;
							}
							jpix = 2 * (jpix % it0->second.hdr.Nside);
							if (ring%2 == 0)
								jpix++;
							if (jpix > it0->second.hdr.Nside)
								jpix -= 2*it0->second.hdr.Nside;
									
							dp[1] = jpix * 0.5 / it0->second.hdr.Nside;
						}
						else
						{
							dp[0] = (4. - sqrt((1.+v1[2]) * 3.)) * Nside_final;
							ring = (1 + (int64_t) sqrt(2 * (12l*it0->second.hdr.Nside*it0->second.hdr.Nside-l) - 0.5)) / 2; // counted from south pole
							jpix = 2 * (((12l*it0->second.hdr.Nside*it0->second.hdr.Nside - 1 - l) - 2 * ring * (ring-1)) % ring) + 1;
							if (jpix > ring)
								jpix -= 2*ring;
									
							dp[1] = -jpix * 0.5 / ring; //-(0.5 * jpix * Nside_final) / it0->second.hdr.Nside;
						}
							
						if (v2[2] > 2./3.)
						{
							dp[0] -= sqrt((1.-v2[2]) * 3.) * Nside_final;
							ring = (1 + (int64_t) sqrt(1.5 + 2 * (double) q)) / 2;
							jpix = 2 * ((q - 2 * ring * (ring-1)) % ring) + 1;
							dp[1] *= ring;
								
							if (dp[1] < 0)
								jpix -= 2*ring;
									
							dp[1] -= 0.5 * jpix;
						}
						else if (v2[2] >= -2./3.)
						{
							dp[0] -= (2. - 1.5 * v2[2]) * Nside_final;
							jpix = q - 2 * Nside_final * (Nside_final + 1);
							ring = Nside_final + 1 + jpix / (4 * Nside_final);
							if (jpix < 0)
							{
								ring = Nside_final;
								jpix += 4 * Nside_final;
							}
							jpix = 2 * (jpix % Nside_final);
							dp[1] *= Nside_final;
								
							if (ring%2 == 0)
								jpix++;
							if (jpix > Nside_final)
								jpix -= 2*Nside_final;
								
							if (jpix > 0 && dp[1] < 0)
								dp[1] += Nside_final;
							else if (jpix < 0 && dp[1] > 0)
								dp[1] -= Nside_final;
									
							dp[1] -= 0.5 * jpix;
						}
						else
						{
							dp[0] -= (4. - sqrt((1.+v2[2]) * 3.)) * Nside_final;
							ring = (1 + (int64_t) sqrt(2 * (12l*Nside_final*Nside_final-q) - 0.5)) / 2; // counted from south pole
							jpix = 2 * (((12l*Nside_final*Nside_final - 1 - q) - 2 * ring * (ring-1)) % ring) + 1;
							dp[1] *= ring;
								
							if (jpix > ring)
								jpix -= 2*ring;
								
							if (jpix < 0 && dp[1] < 0)
								dp[1] += ring;
							else if (jpix > 0 && dp[1] > 0)
								dp[1] -= ring;
									
							dp[1] += 0.5 * jpix;
						}
							
						for (int m = outcnt; m < numoutputs; m++)
						{
							if (pixgrad(map_phi_final[m], (int64_t) Nside_final, Npix_final, q, grad_phi))							
								map_phi_temp[m-outcnt][l] = map_phi_final[m][q] + grad_phi[0] * dp[0] + grad_phi[1] * dp[1] + 0.5 * grad_phi[2] * dp[0] * dp[0] + grad_phi[3] * dp[0] * dp[1] + 0.5 * grad_phi[4] * dp[1] * dp[1];
							else
								map_phi_temp[m-outcnt][l] = -1.6375e30;
						}
						map_isw_temp[l] = map_isw_final[q] + grad_isw[0] * dp[0] + grad_isw[1] * dp[1] + 0.5 * grad_isw[2] * dp[0] * dp[0] + grad_isw[3] * dp[0] * dp[1] + 0.5 * grad_isw[4] * dp[1] * dp[1];
						map_shapiro_temp[l] = map_shapiro_final[q] + grad_shapiro[0] * dp[0] + grad_shapiro[1] * dp[1] + 0.5 * grad_shapiro[2] * dp[0] * dp[0] + grad_shapiro[3] * dp[0] * dp[1] + 0.5 * grad_shapiro[4] * dp[1] * dp[1];
					}
					else
					{
						for (int m = 0; m < numoutputs-outcnt; m++)
							map_phi_temp[m][l] = -1.6375e30;
						map_isw_temp[l] = -1.6375e30;
						map_shapiro_temp[l] = -1.6375e30;
					}
				}
					
					
				free(map_isw_final);
				free(map_shapiro_final);
				for (int m = outcnt; m < numoutputs; m++)
				{
					free(map_phi_final[m]);
					map_phi_final[m] = map_phi_temp[m-outcnt];
				}
				map_isw_final = map_isw_temp;
				map_shapiro_final = map_shapiro_temp;
					
				Nside_final = it0->second.hdr.Nside;
			}

			if (it0->second.hdr.distance <= 0)
			{
				continue;
			}
			
			if (it0->second.hdr.distance != it1->second.hdr.distance || it0->second.hdr.Nside != it1->second.hdr.Nside || it0->second.hdr.Npix != it1->second.hdr.Npix)
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": map properties do not match on iteration " << cnt << "! distances: " << it0->second.hdr.distance << ", " << it1->second.hdr.distance << "; Nside: " << it0->second.hdr.Nside << ", " << it1->second.hdr.Nside << "; Npix: " << it0->second.hdr.Npix << ", " << it1->second.hdr.Npix << endl << endl;
				return -1;
			}
			
			Npix_final = it0->second.hdr.Npix;

#pragma omp parallel for
			for (int l = 0; l < Npix_final; l++)
			{
				if(map_isw_final[l] < -1.5e29 || map_shapiro_final[l] < -1.5e29)
				{
					continue;
				}
				else if(it0->second.pixel[l] < -1.5e29)
				{
					for (int m = outcnt; m < numoutputs; m++)
						map_phi_final[m][l] = it0->second.pixel[l];
					map_isw_final[l] = it0->second.pixel[l];
					map_shapiro_final[l] = it0->second.pixel[l];
				}
				else if(it1->second.pixel[l] < -1.5e29)
				{
					for (int m = outcnt; m < numoutputs; m++)
						map_phi_final[m][l] = it1->second.pixel[l];
					map_isw_final[l] = it1->second.pixel[l];
					map_shapiro_final[l] = it1->second.pixel[l];
				}
				else
				{
					for (int m = outcnt; m < numoutputs; m++)
						map_phi_final[m][l] -= 2.*(it0->second.hdr.distance - dist)*((distances[m]-it0->second.hdr.distance)/(distances[m]*it0->second.hdr.distance))*lin_int(it0->second.pixel[l],it1->second.pixel[l], (it0->second.hdr.distance - (tauobs - back[step].tau))/((back[step].tau - back[step+1].tau)));
					map_isw_final[l] -= 2.*(it0->second.hdr.distance - dist) * (it0->second.pixel[l] - it1->second.pixel[l]) / (back[step].tau - back[step+1].tau);
					map_shapiro_final[l] += 2.*sim.boxsize*(it0->second.hdr.distance - dist) * lin_int(it0->second.pixel[l],it1->second.pixel[l], (it0->second.hdr.distance - (tauobs - back[step].tau))/((back[step].tau - back[step+1].tau)));
				}
			}
			
			if (it0->second.hdr.distance >= distances[outcnt])
			{
				cout << " source distance reached of " << distances[outcnt]*sim.boxsize << " Mpc/h!" << endl << endl;

#pragma omp parallel for	
				for (long l = Npix_final; l < 12l * Nside_final * Nside_final; l++)
					map_phi_final[outcnt][l] = -1.6375e30;
	
				if (numoutputs > 1)
					sprintf(outputfile, "%s%s_lensingphi_%d.fits", sim.output_path, sim.basename_lightcone, outcnt);	
				else
					sprintf(outputfile, "%s%s_lensingphi.fits", sim.output_path, sim.basename_lightcone);	
				write_healpix_map(map_phi_final[outcnt], Nside_final, outputfile, 0, &coordsys);
	
#pragma omp parallel for
				for (long l = 0; l < Npix_final; l++)
					map_phi_final[outcnt][l] = (map_isw_final[l] < -1.5e29) ? map_isw_final[l] : map_isw_final[l] - 2.*(distances[outcnt]+0.5*dist-1.5*it0->second.hdr.distance) * (it0->second.pixel[l] - it1->second.pixel[l]) / (back[step].tau - back[step+1].tau);
	
				if (numoutputs > 1)
					sprintf(outputfile, "%s%s_isw_%d.fits", sim.output_path, sim.basename_lightcone, outcnt);
				else
					sprintf(outputfile, "%s%s_isw.fits", sim.output_path, sim.basename_lightcone);	
				write_healpix_map(map_phi_final[outcnt], Nside_final, outputfile, 0, &coordsys);

				if (it0 == phi0.healpix_data.begin() || it1 == phi1.healpix_data.begin())
				{
					cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": evaluating potential at boundary of covered range" << endl;
#pragma omp parallel for
					for (long l = 0; l < Npix_final; l++)
						map_phi_final[outcnt][l] = (map_isw_final[l] < -1.5e29) ? map_isw_final[l] : pot_obs - lin_int(it0->second.pixel[l],it1->second.pixel[l], (distances[outcnt] - (tauobs - back[step].tau))/((back[step].tau - back[step+1].tau)));
				}
				else if (std::prev(it0)->second.hdr.Nside != Nside_final)
				{
#pragma omp parallel for private(v1,q)
					for (long l = 0; l < Npix_final; l++)
					{
						pix2vec_ring64(Nside_final, l, v1);
						vec2pix_ring64(std::prev(it0)->second.hdr.Nside, v1, &q);
						map_phi_final[outcnt][l] = (map_isw_final[l] < -1.5e29) ? map_isw_final[l] : pot_obs - ((it0->second.hdr.distance-distances[outcnt]) * lin_int(it0->second.pixel[l],it1->second.pixel[l], (distances[outcnt] - (tauobs - back[step].tau))/((back[step].tau - back[step+1].tau))) + (distances[outcnt]-dist) * lin_int(std::prev(it0)->second.pixel[q],std::prev(it1)->second.pixel[q], (distances[outcnt] - (tauobs - back[step].tau))/((back[step].tau - back[step+1].tau)))) / (it0->second.hdr.distance - dist);
					}
				}
				else
				{
#pragma omp parallel for
					for (long l = 0; l < Npix_final; l++)
						map_phi_final[outcnt][l] = (map_isw_final[l] < -1.5e29) ? map_isw_final[l] : pot_obs - ((it0->second.hdr.distance-distances[outcnt]) * lin_int(it0->second.pixel[l],it1->second.pixel[l], (distances[outcnt] - (tauobs - back[step].tau))/((back[step].tau - back[step+1].tau))) + (distances[outcnt]-dist) * lin_int(std::prev(it0)->second.pixel[l],std::prev(it1)->second.pixel[l], (distances[outcnt] - (tauobs - back[step].tau))/((back[step].tau - back[step+1].tau)))) / (it0->second.hdr.distance - dist);
				}
	
				if (numoutputs > 1)
					sprintf(outputfile, "%s%s_potential_%d.fits", sim.output_path, sim.basename_lightcone, outcnt);
				else
					sprintf(outputfile, "%s%s_potential.fits", sim.output_path, sim.basename_lightcone);	
				write_healpix_map(map_phi_final[outcnt], Nside_final, outputfile, 0, &coordsys);
	
				#pragma omp parallel for
				for (long l = 0; l < Npix_final; l++)
					map_phi_final[outcnt][l] = (map_shapiro_final[l] < -1.5e29) ? map_shapiro_final[l] : map_shapiro_final[l] + 2.*sim.boxsize*(distances[outcnt]+0.5*dist-1.5*it0->second.hdr.distance) * lin_int(it0->second.pixel[l],it1->second.pixel[l], (it0->second.hdr.distance - (tauobs - back[step].tau))/((back[step].tau - back[step+1].tau)));
				
				if (numoutputs > 1)
					sprintf(outputfile, "%s%s_shapiro_%d.fits", sim.output_path, sim.basename_lightcone, outcnt);	
				else
					sprintf(outputfile, "%s%s_shapiro.fits", sim.output_path, sim.basename_lightcone);
				write_healpix_map(map_phi_final[outcnt], Nside_final, outputfile, 0, &coordsys);
				free(map_phi_final[outcnt]);
			
				outcnt++;
			}
			
			dist = it0->second.hdr.distance;
			
			if (cnt == thresh && it0->second.hdr.Npix == 12 * it0->second.hdr.Nside * it0->second.hdr.Nside)
			{
				for (int m = outcnt; m < numoutputs; m++)
				{
					monopole = 0;
				
#pragma omp parallel for reduction(+:monopole)			
					for (int l = 0; l < it0->second.hdr.Npix; l++)
					{
						if (map_phi_final[m][l] > -1.5e29)
							monopole += map_phi_final[m][l];
					}
					
					monopole /= (double) it0->second.hdr.Npix;

#pragma omp parallel for				
					for (int l = 0; l < it0->second.hdr.Npix; l++)
						map_phi_final[m][l] -= monopole;
				}
				
				thresh *= 2;
			}
			
			cnt++;
		}

		phi0.clear();
		phi1.clear();

		step++;
		
		if (step > numlines-2)
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": reached cycle " << back[step].cycle << ", no further integration steps possible!" << endl;
			return -1;
		}
		
		phi0.cinfo = back[step];
		phi1.cinfo = back[step+1];
	}
	
	free(map_isw_final);
	free(map_shapiro_final);
	
	cout << COLORTEXT_GREEN << " normal completion." << COLORTEXT_RESET << endl;

	return 0;
}


int loadHealpixData(metric_container * field, double min_dist, double max_dist, char * lightconeparam)
{
	metric_data metric;
	char filename[1024];
	char tokens[128];
	char * token = NULL;
	uint32_t blocksize[2];
	FILE * infile = NULL;
	FILE * infile2 = NULL;
	int count;
	long backtrack;
	double vec[3];
	int64_t j, q;
	int ring;
	int pixbatch_delim[3];
	int pixbatch_size[3] = {0, 0, 0};


	metric.pixel = NULL;

	if (!field->healpix_data.empty())
	{
		if (field->healpix_data.begin()->second.hdr.distance < min_dist)
		{
			if ((min_dist = field->healpix_data.rbegin()->second.hdr.distance) > max_dist)
				return field->healpix_data.rbegin()->first;
		}	
		else if (field->healpix_data.begin()->second.hdr.distance > max_dist)
			max_dist = field->healpix_data.begin()->second.hdr.distance;
	}
	
	if (lightconeparam != NULL)
	{
		strcpy(tokens, lightconeparam);
		token = strtok(tokens, ",");
		do
		{
			sprintf(filename, "%s%s%d_%04d_%s.map", field->dir, field->basename, atoi(token), field->cinfo.cycle, field->name);
			infile = fopen(filename, "rb");
			token = strtok(NULL, ",");
		}
		while (infile == NULL && token != NULL);
	}
	else
	{
		sprintf(filename, "%s%s_%04d_%s.map", field->dir, field->basename, field->cinfo.cycle, field->name);
		infile = fopen(filename, "rb");
	}

	if (infile == NULL)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": map file " << filename << " could not be opened for reading!" << endl;
		return -1;
	}

	if (fread(blocksize+1, sizeof(uint32_t), 1, infile) != 1)
	{
		fclose(infile);
		return -1;
	}

	for (count = 0; !feof(infile) && !ferror(infile); count++)
	{
		if (blocksize[1] != 256)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid header block size in map file " << filename << "!" << endl;
			return -1;
		}

		if (fread(&metric.hdr, sizeof(metric.hdr), 1, infile) != 1)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read header block in map file " << filename << "!" << endl;
			return -1;
		}

		if (metric.hdr.header_blocksize != 256)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid header block size in map file " << filename << "!" << endl;
			return -1;
		}

		if (metric.hdr.data_blocksize != metric.hdr.precision * metric.hdr.Npix)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid data block size in map file " << filename << "!" << endl;
			return -1;
		}

		if (metric.hdr.distance > min_dist) break;

		backtrack = ftell(infile) - (256 + 2 * sizeof(uint32_t));
		
		if (infile2 != NULL)
		{
			fclose(infile2);
			infile2 = NULL;
		}

		if (fseek(infile, metric.hdr.data_blocksize, SEEK_CUR))
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to skip data block in map file " << filename << "!" << endl;
			fclose(infile);
			return -1;
		}

		if (fread(blocksize, sizeof(uint32_t), 2, infile) != 2)
		{	
			if (token != NULL && feof(infile))
			{
				infile2 = infile;
				sprintf(filename, "%s%s%d_%04d_%s.map", field->dir, field->basename, atoi(token), field->cinfo.cycle, field->name);
				infile = fopen(filename, "rb");
				token = strtok(NULL, ",");
				if (infile == NULL)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": map file " << filename << " could not be opened for reading!" << endl;
					return -1;
				}
				if (fread(blocksize+1, sizeof(uint32_t), 1, infile) != 1)
				{
					fclose(infile);
					return -1;
				}
			}
			else
			{
				fclose(infile);
				return -1;
			}
		}
	}

	if (feof(infile) || ferror(infile))
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": read error occured in map file " << filename << "!" << endl;
		fclose(infile);
		return -1;
	}

	if (count > 0)
	{
		if (infile2 != NULL)
		{
			fclose(infile);
			infile = infile2;
		}
		
		if (fseek(infile, backtrack, SEEK_SET))
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to backtrack in map file " << filename << "!" << endl;
			fclose(infile);
			return -1;
		}

		if (fread(&metric.hdr, sizeof(metric.hdr), 1, infile) != 1)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read header block in map file " << filename << "!" << endl;
			return -1;
		}

		if (metric.hdr.header_blocksize != 256)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid header block size in map file " << filename << "!" << endl;
			return -1;
		}

		if (metric.hdr.data_blocksize != metric.hdr.precision * metric.hdr.Npix)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid data block size in map file " << filename << "!" << endl;
			return -1;
		}

		count--;
	}

	while (true)
	{
		if (field->healpix_data.find(count) == field->healpix_data.end()) // data not present
		{
			if (metric.hdr.Nside < 2)
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid Nside = " << metric.hdr.Nside << " in map file " << filename << "!" << endl;
				fclose(infile);
				return -1;
			}

			if (metric.hdr.distance < 0)
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid distance = " << metric.hdr.distance << " in map file " << filename << "!" << endl;
				fclose(infile);
				return -1;
			}

			metric.pixel = (float *) malloc(metric.hdr.Npix * sizeof(float));

			if (metric.hdr.Nside_ring > 0 && metric.hdr.Nside_ring < metric.hdr.Nside)
			{	
				if ((long) metric.hdr.Npix <= 2 * (long) metric.hdr.Nside * (metric.hdr.Nside + 1))
					ring = (int) floor((sqrt(2. * metric.hdr.Npix + 1.01) - 1.) / 2.);
				else if ((long) metric.hdr.Npix <= 2 * (long) metric.hdr.Nside * (metric.hdr.Nside + 1) + 4 * (2 * metric.hdr.Nside - 1) * (long) metric.hdr.Nside)
					ring = ((metric.hdr.Npix - 2 * metric.hdr.Nside * (metric.hdr.Nside + 1)) / 4 / metric.hdr.Nside) + metric.hdr.Nside;
				else if ((long) metric.hdr.Npix < 12 * (long) metric.hdr.Nside * metric.hdr.Nside)
				{
					ring = 12 * (long) metric.hdr.Nside * metric.hdr.Nside - (long) metric.hdr.Npix;
					ring = (int) floor((sqrt(2. * ring + 1.01) - 1.) / 2.);
					ring = 4 * metric.hdr.Nside - 1 - ring;
				}
				else
					ring = 4 * metric.hdr.Nside - 1;
					
				pixbatch_size[0] = (metric.hdr.Nside / metric.hdr.Nside_ring);
						
				pixbatch_delim[1] = ring / pixbatch_size[0];
				pixbatch_delim[0] = (pixbatch_delim[1] > 0) ? pixbatch_delim[1]-1 : 0;
				pixbatch_delim[2] = pixbatch_delim[1]+1;
				pixbatch_size[1] = (pixbatch_size[0] * (pixbatch_size[0]+1) + (2*pixbatch_size[0] - 1 - ring%pixbatch_size[0]) * (ring%pixbatch_size[0])) / 2;
				pixbatch_size[2] = ((ring%pixbatch_size[0] + 1) * (ring%pixbatch_size[0])) / 2;
				pixbatch_size[0] *= pixbatch_size[0];
				for (int p = 0; p < 3; p++)
				{
					if (pixbatch_delim[p] <= metric.hdr.Nside_ring)
						pixbatch_delim[p] = 2 * pixbatch_delim[p] * (pixbatch_delim[p]+1);
					else if (pixbatch_delim[p] <= 3 * metric.hdr.Nside_ring)
						pixbatch_delim[p] = 2 * metric.hdr.Nside_ring * (metric.hdr.Nside_ring+1) + (pixbatch_delim[p]-metric.hdr.Nside_ring) * 4 * metric.hdr.Nside_ring;
					else if (pixbatch_delim[p] < 4 * metric.hdr.Nside_ring)
						pixbatch_delim[p] = 12 * metric.hdr.Nside_ring * metric.hdr.Nside_ring - 2 * (4 * metric.hdr.Nside_ring - 1 - pixbatch_delim[p]) * (4 * metric.hdr.Nside_ring - pixbatch_delim[p]);
					else
						pixbatch_delim[p] = 12 * metric.hdr.Nside_ring * metric.hdr.Nside_ring;
				}
					
				if (metric.hdr.precision == 4)
				{
					float * fpix = (float *) malloc (metric.hdr.Npix * sizeof(float));
					if (fread(fpix, sizeof(float), metric.hdr.Npix, infile) != metric.hdr.Npix)
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read data block in map file " << filename << "!" << endl;
						fclose(infile);
						free(metric.pixel);
						return -1;
					}

#pragma omp parallel for private(j) collapse(2)
					for (int p = 0; p < pixbatch_delim[0]; p++)
					{
						for (int i = 0; i < pixbatch_size[0]; i++)
						{
							ring2nest64(metric.hdr.Nside_ring, p, &j);
							j = j*pixbatch_size[0] + i;
							nest2ring64(metric.hdr.Nside, j, &j);
							metric.pixel[j] = fpix[pixbatch_size[0]*p + i];
						}
					}
#pragma omp parallel for private(q,j)
					for (int p = pixbatch_delim[0]; p < pixbatch_delim[1]; p++)
					{
						q = 0;
						for (int i = 0; i < pixbatch_size[0]; i++)
						{
							ring2nest64(metric.hdr.Nside_ring, p, &j);
							j = j*pixbatch_size[0] + i;
							nest2ring64(metric.hdr.Nside, j, &j);
							if (j < metric.hdr.Npix)
							{
								metric.pixel[j] = fpix[pixbatch_size[0]*pixbatch_delim[0] + pixbatch_size[1]*(p-pixbatch_delim[0]) + q];
								q++;
							}
						}
					}
#pragma omp parallel for private(q,j)
					for (int p = pixbatch_delim[1]; p < pixbatch_delim[2]; p++)
					{
						q = 0;
						for (int i = 0; i < pixbatch_size[0]; i++)
						{
							ring2nest64(metric.hdr.Nside_ring, p, &j);
							j = j*pixbatch_size[0] + i;
							nest2ring64(metric.hdr.Nside, j, &j);
							if (j < metric.hdr.Npix)
							{
								metric.pixel[j] = fpix[pixbatch_size[0]*pixbatch_delim[0] + pixbatch_size[1]*(pixbatch_delim[1]-pixbatch_delim[0]) + pixbatch_size[2]*(p-pixbatch_delim[1]) + q];
								q++;
							}
						}
					}
	
					free(fpix);
				}
				else if (metric.hdr.precision == 8)
				{
					double * dpix = (double *) malloc (metric.hdr.Npix * sizeof(double));
					if (fread(dpix, sizeof(double), metric.hdr.Npix, infile) != metric.hdr.Npix)
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read data block in map file " << filename << "!" << endl;
						fclose(infile);
						free(metric.pixel);
						return -1;
					}
	
#pragma omp parallel for private(j) collapse(2)
					for (int p = 0; p < pixbatch_delim[0]; p++)
					{
						for (int i = 0; i < pixbatch_size[0]; i++)
						{
							ring2nest64(metric.hdr.Nside_ring, p, &j);
							j = j*pixbatch_size[0] + i;
							nest2ring64(metric.hdr.Nside, j, &j);
							metric.pixel[j] = dpix[pixbatch_size[0]*p + i];
						}
					}
#pragma omp parallel for private(q,j)
					for (int p = pixbatch_delim[0]; p < pixbatch_delim[1]; p++)
					{
						q = 0;
						for (int i = 0; i < pixbatch_size[0]; i++)
						{
							ring2nest64(metric.hdr.Nside_ring, p, &j);
							j = j*pixbatch_size[0] + i;
							nest2ring64(metric.hdr.Nside, j, &j);
							if (j < metric.hdr.Npix)
							{
								metric.pixel[j] = dpix[pixbatch_size[0]*pixbatch_delim[0] + pixbatch_size[1]*(p-pixbatch_delim[0]) + q];
								q++;
							}
						}
					}
#pragma omp parallel for private(q,j)
					for (int p = pixbatch_delim[1]; p < pixbatch_delim[2]; p++)
					{
						q = 0;
						for (int i = 0; i < pixbatch_size[0]; i++)
						{
							ring2nest64(metric.hdr.Nside_ring, p, &j);
							j = j*pixbatch_size[0] + i;
							nest2ring64(metric.hdr.Nside, j, &j);
							if (j < metric.hdr.Npix)
							{
								metric.pixel[j] = dpix[pixbatch_size[0]*pixbatch_delim[0] + pixbatch_size[1]*(pixbatch_delim[1]-pixbatch_delim[0]) + pixbatch_size[2]*(p-pixbatch_delim[1]) + q];
								q++;
							}
						}
					}
	
					free(dpix);
				}
				else
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": precision " << metric.hdr.precision << " bytes not supported for map files!" << endl;
					free(metric.pixel);
				}
			}
			else
			{
				if (metric.hdr.precision == 4)
				{
					if (fread(metric.pixel, sizeof(float), metric.hdr.Npix, infile) != metric.hdr.Npix)
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read data block in map file " << filename << "!" << endl;
						fclose(infile);
						free(metric.pixel);
						return -1;
					}
				}
				else if (metric.hdr.precision == 8)
				{
					double * dpix = (double *) malloc (metric.hdr.Npix * sizeof(double));
					if (fread(dpix, sizeof(double), metric.hdr.Npix, infile) != metric.hdr.Npix)
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read data block in map file " << filename << "!" << endl;
						fclose(infile);
						free(metric.pixel);
						return -1;
					}
					for (int p = 0; p < metric.hdr.Npix; p++)
						metric.pixel[p] = dpix[p];
					free(dpix);
				}
				else
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": precision " << metric.hdr.precision << " bytes not supported for map files!" << endl;
					free(metric.pixel);
				}
			}	

			field->healpix_data.insert(std::pair<int,metric_data>(count, metric));
		}
		else
		{
			if (fseek(infile, metric.hdr.data_blocksize, SEEK_CUR))
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to skip data block in map file " << filename << "!" << endl;
				fclose(infile);
				return -1;
			}
		}

		if (metric.hdr.distance > max_dist) break;

		if (fread(blocksize, sizeof(uint32_t), 2, infile) != 2)
		{
			if (infile2 != NULL)
			{
				fclose(infile2);
				infile2 = NULL;
				infile = fopen(filename, "rb");
				if (infile == NULL)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": map file " << filename << " could not be opened for reading!" << endl;
					return -1;
				}
				if (fread(blocksize+1, sizeof(uint32_t), 1, infile) != 1)
				{
					fclose(infile);
					return -1;
				}
			}
			else if (token != NULL && feof(infile))
			{
				fclose(infile);
				sprintf(filename, "%s%s%d_%04d_%s.map", field->dir, field->basename, atoi(token), field->cinfo.cycle, field->name);
				infile = fopen(filename, "rb");
				token = strtok(NULL, ",");
				if (infile == NULL)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": map file " << filename << " could not be opened for reading!" << endl;
					return -1;
				}
				if (fread(blocksize+1, sizeof(uint32_t), 1, infile) != 1)
				{
					fclose(infile);
					return -1;
				}
			}
			else
			{
				fclose(infile);
				return -2;
			}
		}

		if (blocksize[1] != 256)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid header block size in map file " << filename << "!" << endl;
			return -1;
		}

		if (fread(&metric.hdr, sizeof(metric.hdr), 1, infile) != 1)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read header block in map file " << filename << "!" << endl;
			return -1;
		}

		if (metric.hdr.header_blocksize != 256)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid header block size in map file " << filename << "!" << endl;
			return -1;
		}

		if (metric.hdr.data_blocksize != metric.hdr.precision * metric.hdr.Npix)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid data block size in map file " << filename << "!" << endl;
			return -1;
		}

		count++;
	}

	fclose(infile);

	return count;
}

bool pixgrad(float * pixel, const int64_t Nside, const int64_t Npix, int64_t ipix, double * result)
{
	int64_t j, k, l, m, q, ring;
	float w1, w2, dring = 2.;

	if (ipix >= Npix || pixel[ipix] < -1e30) return false;
	
	if (ipix < Nside * (Nside + 1) * 2l) // in north polar cap
	{
		ring = (1 + (int64_t) sqrt(1.5 + 2 * (double) ipix)) / 2;
		j = ipix - 2 * ring * (ring-1);
		q = j / ring;
		j %= ring;
		
		// j-derivative
		k = ipix+1;
		l = ipix-1;
		if (q == 3 && j == ring-1)
			k -= 4*ring;
		if (q == 0 && j == 0)
			l += 4*ring;
			
		if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
		
		result[1] = (pixel[k] - pixel[l]) / 2.;
		
		result[4] = pixel[k] + pixel[l] - 2. * pixel[ipix];
		
		// ring derivative
		if (ring == Nside)
		{
			k = Nside * (Nside+1) * 2l + q * Nside + j;
			l = k+1;
			if (q == 3 && j == Nside-1)
				l -= 4*Nside;
				
			if (k >= Npix || l >= Npix)
			{
				result[0] = 0;//pixel[ipix];
				dring = 1.;
				k = (j == 0 && q == 0) ? ipix-1+4*Nside : ipix-1;
				l = (q == 3 && j == Nside-1) ? ipix-1+4*Nside : ipix+1;
				if (pixel[k] < -1e30 || pixel[l] < -1e30)
					return false;
				result[3] = (pixel[l] - pixel[k]) / 2.;
			}
			else if (pixel[k] < -1e30 || pixel[l] < -1e30)
				return false;
			else
			{
				result[0] = (0.5 * (pixel[k] + pixel[l]) - pixel[ipix]);// * 8. / 3.;
				result[2] = 0.5 * (pixel[k] + pixel[l]);
				result[3] = pixel[l] - pixel[k];
				w1 = 0.125;
			}
		}
		else
		{
			k = ring * (ring+1) * 2l + q * (ring+1) + j;
			l = k+1;
			if (j < ring/2)
			{
				m = k-1;
				if (q == 0 && j == 0)
					m += 4*(ring+1);
			}
			else
			{
				m = l+1;
				if (q == 3 && j == ring-1)
					m -= 4*(ring+1);
			}
			
			if (k >= Npix || l >= Npix || m >= Npix)
			{
				result[0] = pixel[ipix];
				dring = 1.;
				k = (j == 0 && q == 0) ? ipix-1+4*ring : ipix-1;
				l = (q == 3 && j == ring-1) ? ipix-1+4*ring : ipix+1;
				if (pixel[k] < -1e30 || pixel[l] < -1e30)
					return false;
				result[3] = (pixel[l] - pixel[k]) / 2.;
				w1 = 0;
			}
			else if (pixel[k] < -1e30 || pixel[l] < -1e30 || pixel[m] < -1e30)
				return false;
			else
			{
				result[0] = (1. - (j+0.5)/ring) * pixel[k] + ((j+0.5)/ring) * pixel[l];
				if (j < ring/2)
					result[3] = ((0.5 + (j+0.5)/ring) * (pixel[l] - pixel[k]) + (0.5 - (j+0.5)/ring) * (pixel[k] - pixel[m])) * (ring+1) * (ring+1) / ring / ring;
				else
					result[3] = ((1.5 - (j+0.5)/ring) * (pixel[l] - pixel[k]) + ((j+0.5)/ring - 0.5) * (pixel[m] - pixel[l])) * (ring+1) * (ring+1) / ring / ring;
					
				w1 = (ring-j-0.5) * (j+0.5) * 0.5 / (ring+1) / (ring+1);
				
				/*if (ring == 1)
				{
					result[1] *= 1. - sqrt(2.) * sin(M_PI / 8.);
					result[1] += 0.5 * (pixel[l] - pixel[k]) / sqrt(2.);
					result[1] *= M_PI / 2.;
				}*/
			}
			
			result[2] = result[0];
		}
		
		if (ring == 1)
		{
			if (pixel[0] < -1e30 || pixel[1] < -1e30 || pixel[2] < -1e30 || pixel[3] < -1e30) return false;
			
			/*result[0] -= 0.25 * (pixel[0] + pixel[1] + pixel[2] + pixel[3]);
			result[2] += 0.25 * (pixel[0] + pixel[1] + pixel[2] + pixel[3]);
			w2 = 0;*/
			result[0] = result[0]/0.75 - pixel[ipix] - pixel[(ipix+2)%4]/3.;
			result[2] = result[2]/1.5 + pixel[ipix] + pixel[(ipix+2)%4]/3.;
			w2 = -1./96.;
			
			/*w1 = 0;
			w2 = 0;
			result[0] /= 0.75;
			result[0] -= (2. - sqrt(2.)) * (pixel[0] + pixel[1] + pixel[2] + pixel[3]) / 1.5 + ((sqrt(2.) + cos(M_PI/8.)) / 0.75 - 3.) * pixel[ipix] + ((sqrt(2.) - cos(M_PI/8.)) / 0.75 - 1.) * pixel[(ipix+2)%4];
			result[2] /= 1.5;
			result[2] += (sqrt(2.) - 2.) * (pixel[0] + pixel[1] + pixel[2] + pixel[3]) / 3. - ((sqrt(2.) + cos(M_PI/8.)) / 1.5 - 3.) * pixel[ipix] - ((sqrt(2.) - cos(M_PI/8.)) / 1.5 - 1.) * pixel[(ipix+2)%4];
			result[4] -= result[2] - 2. * pixel[ipix] - result[0];
			result[4] *= M_PI * M_PI / 4.;*/
			//result[2] -= result[0] * 0.5 / (6l * Nside * Nside - 1);
		}
		else
		{
			k = (ring-2) * (ring-1) * 2l + q * (ring-1) + j;
			l = k-1;
			if (j < ring/2)
				m = (k != 3) ? k+1 : 0;
			else if (q == 0 && j < 2)
				m = l-1+4*(ring-1);
			else
				m = l-1;
			if (q == 0 && j == 0)
				l += 4*(ring-1);
			if (q == 3 && j == ring-1)
				k -= 4*(ring-1);
				
			if (pixel[k] < -1e30 || pixel[l] < -1e30  || pixel[m] < -1e30) return false;
			
			w2 = (ring-j-0.5) * (j+0.5) * 0.5 / (ring-1) / (ring-1);
			
			if (ring == Nside)
			{
				//result[2] += (result[0] - w1 * result[4]) / (2.5 * Nside);
				result[0] += pixel[ipix] - ((1. - (j+0.5)/ring) * pixel[k] + ((j+0.5)/ring) * pixel[l]);
				//result[2] -= (w2 * result[4] + pixel[ipix] - ((1. - (j+0.5)/ring) * pixel[k] + ((j+0.5)/ring) * pixel[l])) / (10 * Nside);
			}
			else
			{
				result[0] -= (1. - (j+0.5)/ring) * pixel[k] + ((j+0.5)/ring) * pixel[l];
				//result[2] -= (result[0] - (w1 - w2) * result[4]) * ring * 0.5 / (6l * Nside * Nside - ring * ring);
			}
			result[2] += (1. - (j+0.5)/ring) * pixel[k] + ((j+0.5)/ring) * pixel[l];
			if (j < ring/2)
				result[3] -= ((0.5 + (j+0.5)/ring) * (pixel[k] - pixel[l]) + (0.5 - (j+0.5)/ring) * (pixel[m] - pixel[k])) * (ring-1) * (ring-1) / ring / ring;
			else
				result[3] -= ((1.5 - (j+0.5)/ring) * (pixel[k] - pixel[l]) + ((j+0.5)/ring - 0.5) * (pixel[l] - pixel[m])) * (ring-1) * (ring-1) / ring / ring;
				
			//result[3] -= ((dring-1.) * 0.5 * ring / (ring == Nside ? Nside : ring+1) - 0.5 * ring / (ring-1)) * result[4];
		}
		
		if (dring > 1.)
		{
			result[0] -= (w1 - w2) * result[4];
			result[0] /= dring;
			result[2] -= 2. * pixel[ipix] + (w1 + w2) * result[4];
			result[3] -= 2. * result[1] / ring;
			result[3] /= dring;
		}
		else
		{
			result[2] = 0.;
			result[3] = (ring > 1) ? (ring * result[3] - result[1]) / (ring-1) : 0.;
		}
			
		//result[4] += result[0] * ring * (1. - (0.5 - ring * ring / 18. / Nside / Nside) * ring * ring / Nside / Nside);
		//result[3] -= result[1] * (6l * Nside * Nside - 2l * ring * ring) / ring / (6l * Nside * Nside - ring * ring);
	}
	else if (ipix < Nside * (10 * Nside - 2)) // in equatorial region
	{
		j = ipix - 2 * Nside * (Nside + 1);
		ring = Nside + 1 + j / (4 * Nside);
		j %= 4 * Nside;
		
		// j-derivative
		k = (j == 4*Nside-1) ? ipix+1-4*Nside : ipix+1;
		l = (j == 0) ? ipix-1+4*Nside : ipix-1;
		
		if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
		
		result[1] = (pixel[k] - pixel[l]) / 2.;
		
		result[4] = pixel[k] + pixel[l] - 2. * pixel[ipix];
		
		// ring derivative
		k = ipix + 4*Nside;
		l = (ring%2) ? k-1 : k+1;
		if (j == 0 && ring%2)
			l += 4*Nside;
		else if (j == 4*Nside-1 && ring%2 == 0)
			l -= 4*Nside;
			
		if (k >= Npix || l >= Npix)
		{
			result[0] = pixel[ipix];
			dring = 1.;
			
			k = (j == 4*Nside-1) ? ipix+1-4*Nside : ipix+1;
			l = (j == 0) ? ipix-1+4*Nside : ipix-1;
		
			if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
			
			result[2] = 0;
			result[3] = (pixel[k] - pixel[l]) / 2.;
		}
		else if (pixel[k] < -1e30 || pixel[l] < -1e30)
			return false;
		else
		{		
			result[0] = 0.25 * (pixel[k] + pixel[l]);
			result[2] = 2. * result[0];
			result[3] = (ring%2) ? (pixel[k] - pixel[l]) : (pixel[l] - pixel[k]);
		}
		
		k = ipix - 4*Nside;
		l = (ring%2) ? k-1 : k+1;
		if (j == 0 && ring%2)
			l += 4*Nside;
		else if (j == 4*Nside-1 && ring%2 == 0)
			l -= 4*Nside;
		
		if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
		
		result[0] -= 0.5 * (pixel[k] + pixel[l]) / dring;
		result[3] -= (ring%2) ? (pixel[k] - pixel[l]) : (pixel[l] - pixel[k]);
		if (dring > 1.)
		{
			result[2] += 0.5 * (pixel[k] + pixel[l]) - 2. * pixel[ipix] - 0.25 * result[4];
			//result[2] -= result[0] * (8l * Nside - 4l * ring) / (7l * Nside * Nside - 16l * Nside * ring + 4l * ring * ring);
			result[3] /= dring;
		}
		//result[3] += result[1] * (8l * Nside - 4l * ring) / (7l * Nside * Nside - 16l * Nside * ring + 4l * ring * ring);
		//result[4] -= result[0] * (2l * Nside - ring) * (7l * Nside * Nside - 16l * Nside * ring + 4l * ring * ring) / 9. / Nside / Nside;
	}
	else
	{
		ring = (1 + (int64_t) sqrt(2 * (12l*Nside*Nside-ipix) - 0.5)) / 2; // counted from south pole
		j = (12l*Nside*Nside - 1 - ipix) - 2 * ring * (ring-1);
		q = j / ring;
		j %= ring;
		
		// j-derivative
		k = ipix+1;
		l = ipix-1;
		if (q == 0 && j == 0)
			k -= 4*ring;
		if (q == 3 && j == ring-1)
			l += 4*ring;
			
		if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
		
		result[1] = (pixel[k] - pixel[l]) / 2.;
		
		result[4] = pixel[k] + pixel[l] - 2. * pixel[ipix];
		
		// ring derivative
		if (ring == Nside)
		{
			k = Nside * (Nside+1) * 2l + q * Nside + j;
			l = k-1;
			if (q == 0 && j == 0)
				l += 4*Nside;
			
			k = 12l*Nside*Nside - 1 - k;
			l = 12l*Nside*Nside - 1 - l;
				
			if (k >= Npix || l >= Npix)
				return false;
			else if (pixel[k] < -1e30 || pixel[l] < -1e30)
				return false;
			else
			{
				result[0] = (pixel[ipix] - 0.5 * (pixel[k] + pixel[l]));// * 8. / 3.;
				result[2] = 0.5 * (pixel[k] + pixel[l]);
				result[3] = pixel[k] - pixel[l];
				w1 = 0.125;
			}
		}
		else
		{
			k = ring * (ring+1) * 2l + q * (ring+1) + j;
			l = k+1;
			if (j < ring/2)
			{
				m = k-1;
				if (q == 0 && j == 0)
					m += 4*(ring+1);
			}
			else
			{
				m = l+1;
				if (q == 3 && j == ring-1)
					m -= 4*(ring+1);
			}
			
			k = 12l*Nside*Nside - 1 - k;
			l = 12l*Nside*Nside - 1 - l;
			m = 12l*Nside*Nside - 1 - m;
			
			if (k >= Npix || l >= Npix || m >= Npix)
				return false;
			else if (pixel[k] < -1e30 || pixel[l] < -1e30 || pixel[m] < -1e30)
				return false;
			else
			{
				result[0] = -(1. - (j+0.5)/ring) * pixel[k] - ((j+0.5)/ring) * pixel[l];
				if (j < ring/2)
					result[3] = ((0.5 + (j+0.5)/ring) * (pixel[l] - pixel[k]) + (0.5 - (j+0.5)/ring) * (pixel[k] - pixel[m])) * (ring+1) / ring;
				else
					result[3] = ((1.5 - (j+0.5)/ring) * (pixel[l] - pixel[k]) + ((j+0.5)/ring - 0.5) * (pixel[m] - pixel[l])) * (ring+1) / ring;
				w1 = (ring-j-0.5) * (j+0.5) * 0.5 / (ring+1) / (ring+1);
			}
			
			result[2] = -result[0];
		}
		
		if (ring == 1)
		{
			if (pixel[12l*Nside*Nside-1] < -1e30 || pixel[12l*Nside*Nside-2] < -1e30 || pixel[12l*Nside*Nside-3] < -1e30 || pixel[12l*Nside*Nside-4] < -1e30) return false;
			
			result[0] += 0.25 * (pixel[12l*Nside*Nside-1] + pixel[12l*Nside*Nside-2] + pixel[12l*Nside*Nside-3] + pixel[12l*Nside*Nside-4]);
			result[2] += 0.25 * (pixel[12l*Nside*Nside-1] + pixel[12l*Nside*Nside-2] + pixel[12l*Nside*Nside-3] + pixel[12l*Nside*Nside-4]);
			w2 = 0;
			//result[2] += result[0] * 0.5 / (6l * Nside * Nside - 1);
		}
		else
		{
			k = (ring-2) * (ring-1) * 2l + q * (ring-1) + j;
			l = k-1;
			if (j < ring/2)
				m = (k != 3) ? k+1 : 0;
			else if (q == 0 && j < 2)
				m = l-1+4*(ring-1);
			else
				m = l-1;
			if (q == 0 && j == 0)
				l += 4*(ring-1);
			if (q == 3 && j == ring-1)
				k -= 4*(ring-1);
				
			k = 12l*Nside*Nside - 1 - k;
			l = 12l*Nside*Nside - 1 - l;
			m = 12l*Nside*Nside - 1 - m;
				
			if (k >= Npix || l >= Npix || m >= Npix)
				return false;
			else if (pixel[k] < -1e30 || pixel[l] < -1e30|| pixel[m] < -1e30)
				return false;
				
			w2 = (ring-j-0.5) * (j+0.5) * 0.5 / (ring-1) / (ring-1);
				
			if (ring == Nside)
			{
				//result[2] -= (result[0] + w1 * result[4]) / (2.5 * Nside);
				result[0] -= pixel[ipix] - ((1. - (j+0.5)/ring) * pixel[k] + ((j+0.5)/ring) * pixel[l]);
				//result[2] -= (w2 * result[4] + pixel[ipix] - ((1. - (j+0.5)/ring) * pixel[k] + ((j+0.5)/ring) * pixel[l])) / (10 * Nside);
			}
			else
			{
				result[0] += (1. - (j+0.5)/ring) * pixel[k] + ((j+0.5)/ring) * pixel[l];
				//result[2] += (result[0] + (w1 - w2) * result[4]) * ring * 0.5 / (6l * Nside * Nside - ring * ring);
			}
			result[2] += (1. - (j+0.5)/ring) * pixel[k] + ((j+0.5)/ring) * pixel[l];
			if (j < ring/2)
				result[3] -= ((0.5 + (j+0.5)/ring) * (pixel[k] - pixel[l]) + (0.5 - (j+0.5)/ring) * (pixel[m] - pixel[k])) * (ring-1) / ring;
			else
				result[3] -= ((1.5 - (j+0.5)/ring) * (pixel[k] - pixel[l]) + ((j+0.5)/ring - 0.5) * (pixel[l] - pixel[m])) * (ring-1) / ring;
			//result[3] += ((dring-1.) * 0.5 * ring / (ring == Nside ? Nside : ring+1) - 0.5 * ring / (ring-1)) * result[4];
		}
		
		if (dring > 1.)
		{
			result[0] += (w1 - w2) * result[4];
			result[0] /= dring;
			result[2] -= 2. * pixel[ipix] + (w1 + w2) * result[4];
			result[3] /= dring;
		}
		else
			result[2] = 0.;
			
		//result[4] -= result[0] * ring * (1. - (0.5 - ring * ring / 18. / Nside / Nside) * ring * ring / Nside / Nside);
		//result[3] -= result[1] * (6l * Nside * Nside - 2l * ring * ring) / ring / (6l * Nside * Nside - ring * ring);
	}
	
	return true;
}



