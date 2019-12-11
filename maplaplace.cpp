#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "chealpix.h"

using namespace std;

bool laplacian(float * pixel, const int64_t Nside, int64_t ipix, float & result);

int main(int argc, char **argv)
{
	float * map = NULL;
	float * kappa = NULL;
	long Nside;
	long Npix = 0;
	char coordsys;
	char ordering[10];
	
	if (argc < 3) return -1;
	
	map = read_healpix_map(argv[1], &Nside, &coordsys, ordering);
	
	cout << " Nside = " << Nside << ", coord sys = " << coordsys << ", ordering = " << ordering << endl;
	
	kappa = (float *) malloc(sizeof(float) * 12 * Nside * Nside);
	
	/*long j = 12l * Nside * Nside - 4l;
	long ring = 1;
	long q = 0;
	
	for (long i = 0; i < 12l * Nside * Nside; i++)
	{
		map[i] = kappa[j];
		j++;
		q++;
		
		if (i < 2l * Nside * (Nside-1l))
		{
			if (q == 4 * ring)
			{
				j = j - 4 * ring - 4 * (ring+1);
				ring++;
				q = 0;
			}
		}
		else if (i < 2l * Nside * (5l * Nside - 1l))
		{
			if (q == 4 * Nside)
			{
				j = j - 8 * Nside;
				q = 0;
			}
		}
		else if (q == 4 * ring)
		{
			j = j - 4 * ring - 4 * (ring-1);
			ring--;
			q = 0;
		}
	}*/
	
	if (argc > 3) Npix = atoi(argv[3]);
	else Npix = 12 * Nside * Nside;
		
	cout << " computing kappa..." << endl;
		
	for (long i = 0; i < Npix; i++)
	{
		if(!laplacian(map, Nside, i, kappa[i]))
			kappa[i] = -1.6375e30;
	}
	
	free(map);
	
	for (long i = Npix; i < 12 * Nside * Nside; i++)
		kappa[i] = -1.6375e30;
	
	cout << " writing data to " << argv[2] << endl;
	
	write_healpix_map(kappa, Nside, argv[2], 0, &coordsys);
	
	free(kappa);
	
	return 0;
}


bool laplacian(float * pixel, const int64_t Nside, int64_t ipix, float & result)
{
	int64_t j, k, l, q, ring;
	float temp, temp2, w1, w2; //, temp3;

	if (pixel[ipix] < -1e30) return false;

	if (ipix < Nside * (Nside + 1) * 2l) // north polar cap
	{
		ring = (1 + (int64_t) sqrt(1.5 + 2 * (double) ipix)) / 2;
		j = ipix - 2 * ring * (ring-1);
		q = j / ring;
		j %= ring;
		
		// phi-derivative
		k = ipix+1;
		l = ipix-1;
		if (q == 3 && j == ring-1)
			k -= 4*ring;
		if (q == 0 && j == 0)
			l += 4*ring;
		
		if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
	
		//temp = (4. * (pixel[k] + pixel[l]) / 3. - 2.5 * pixel[ipix]) * 4. / M_PI / M_PI;
		temp2 = (pixel[k] + pixel[l] - 2. * pixel[ipix]);
		
		// ring derivative
		if (ring == Nside)
		{
			k = ring * (ring+1) * 2l + q * ring + j;
			l = k+1;
			
			if (q == 3 && j == Nside-1)
				l -= 4*Nside;
		
			if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
		
			result = (0.5 * (pixel[k] + pixel[l]) - pixel[ipix]) * 8. / 3.;
			temp = (0.5 * (pixel[k] + pixel[l]));
			w1 = 0.125;
		}
		else
		{
			k = ring * (ring+1) * 2l + q * (ring+1) + j;
			l = k+1;
		
			if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
		
			result = (1. - (j+0.5)/ring) * pixel[k] + ((j+0.5)/ring) * pixel[l];
			temp = result;
			w1 = (ring-j-0.5) * (j+0.5) * 0.5 / (ring+1) / (ring+1);
		}
	
		/*k = (ring+1) * (ring+2) * 2l + ((q * ring + j) * (ring+2) + 1) / ring;
		l = k+1;
	
		if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
	
		temp2 = (1. + ((2 * (q * ring + j) + 1) / ring) - (2. * (q * ring + j) - 1.) / ring) * pixel[k] - (((2 * (q * ring + j) + 1) / ring) - (2. * (q * ring + j) - 1.) / ring) * pixel[l];*/
	
		if (ring == 1)
		{
			if (pixel[0] < -1e30 || pixel[1] < -1e30 || pixel[2] < -1e30 || pixel[3] < -1e30) return false;
				
			result -= 0.25 * (pixel[0] + pixel[1] + pixel[2] + pixel[3]);
			temp += 0.25 * (pixel[0] + pixel[1] + pixel[2] + pixel[3]);
			//temp3 = pixel[(ipix+2)%4];
			w2 = 0;
		}
		else
		{
			k = (ring-2) * (ring-1) * 2l + q * (ring-1) + j;
			l = k-1;
			if (q == 0 && j == 0)
				l += 4*(ring-1);
			if (q == 3 && j == ring-1)
				k -= 4*(ring-1);
				
			if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
			
			if (ring == Nside)
				result += pixel[ipix] - ((1. - (j+0.5)/ring) * pixel[k] + ((j+0.5)/ring) * pixel[l]);
			else
				result -= (1. - (j+0.5)/ring) * pixel[k] + ((j+0.5)/ring) * pixel[l];
			temp += (1. - (j+0.5)/ring) * pixel[k] + ((j+0.5)/ring) * pixel[l];
			
			w2 = (ring-j-0.5) * (j+0.5) * 0.5 / (ring-1) / (ring-1);
		
			/*if (ring == 2)
			{
				if (pixel[0] < -1e30 || pixel[1] < -1e30 || pixel[2] < -1e30 || pixel[3] < -1e30) return false;
			
				temp3 = 0.25 * (pixel[0] + pixel[1] + pixel[2] + pixel[3]);
			}
			else
			{
				l = ((q * ring + j) * (ring-2) + ring-1) / ring;
				k = (l == 0) ? 2l * (ring-1) * (ring-2) - 1 : 2l * (ring-2) * (ring-3) + l - 1;
				l = (l == 4 * (ring-2)) ? 2l * (ring-2) * (ring-3) : 2l * (ring-2) * (ring-3) + l;
			
				if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
			
				temp3 = ((2. * (q * ring + j) + 1.) / ring - (2 * (q * ring + j) + 1) / ring) * pixel[k] + (1 - (2. * (q * ring + j) + 1.) / ring + (2 * (q * ring + j) + 1) / ring) * pixel[l];
			}*/
		}
		/*result -= (temp2-temp3) / 8;
		result /= 6 * ring;
	
		result += (4. * temp / 3. - (temp2+temp3) / 12. - 2.5 * pixel[ipix]) / 4.;*/
		
		result -= (w1 - w2) * temp2;
	
		// FIXME
		result *= (6 * Nside * Nside - 3 * ring * ring) / 8. / ring;
		result += (6 * Nside * Nside - ring * ring) * (temp - 2. * pixel[ipix] - (w1 + w2) * temp2) / 4.;
		/*
		result /= 8 * ring;
		result += (temp - 2. * pixel[ipix]) / 4.;
	
		// keep this:
		result *= (6 * Nside * Nside - ring * ring);*/
	
		/*if (ring == 1)
		{
			k = (ipix+2) % 4;
			l = k;
		}
		else
		{
			k = (q == 3 && j == ring-2) ? k+1-4*ring : k+1;
			l = (q == 0 && j == 1) ? l-1+4*ring : l-1;
		}
	
		if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
	
		temp += (pixel[k] + pixel[l]) / 3. / M_PI / M_PI;*/
	
		result += 36. * Nside * Nside * Nside * Nside * temp2 / M_PI / M_PI / (6 * Nside * Nside - ring * ring);
	}
	else if (ipix < 2l * Nside * (5l * Nside - 1l)) // equatorial region
	{
		ring = (ipix - 2l * Nside * (Nside-1)) / (4l * Nside); // + Nside
		j = ipix - 2l * Nside * (Nside-1) - 4l * Nside * ring;
		
		k = (j == 4l*Nside-1) ? ipix+1-4l*Nside : ipix+1;
		l = (j == 0) ? ipix+4l*Nside-1 : ipix-1;
		
		if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
		
		temp2 = (pixel[k] + pixel[l] - 2.*pixel[ipix]);
		
		k = ipix + 4l * Nside;
		
		if (ring % 2)
		{
			l = (j == 0) ? k+(4l*Nside-1) : k-1;
			
			if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
			
			result = 0.5 * (pixel[k] + pixel[l]);
			temp = result;
			
			k = ipix - 4l * Nside;
			l = (j == 0) ? ipix-1 : k-1;
		}
		else
		{
			l = (j == 4l*Nside-1) ? ipix+1 : k+1;
			
			if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
			
			result = 0.5 * (pixel[k] + pixel[l]);
			temp = result;
			
			k = ipix - 4l * Nside;
			l = (j == 4l*Nside-1) ? k+1-4l*Nside : k+1;
		}
		
		if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
			
		result -= 0.5 * (pixel[k] + pixel[l]);
		temp += 0.5 * (pixel[k] + pixel[l]);
		
		result *= Nside-ring;
		result += (temp - 2.*pixel[ipix] - 0.25 * temp2) * (2.25*Nside*Nside - (Nside-ring)*(Nside-ring));
		
		result += temp2 * 4. * Nside * Nside / M_PI / M_PI / (1. - (Nside-ring)*(Nside-ring)/2.25/Nside/Nside);
	}
	else // south polar cap
	{
		ring = (1 + (int64_t) sqrt(1.5 + 2 * (double) (12l*Nside*Nside-1 - ipix))) / 2;
		j = 12l*Nside*Nside-1 - ipix - 2 * ring * (ring-1);
		q = j / ring;
		j %= ring;
		
		// phi-derivative
		k = ipix+1;
		l = ipix-1;
		if (q == 3 && j == ring-1)
			l += 4*ring;
		if (q == 0 && j == 0)
			k -= 4*ring;
		
		if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
	
		temp2 = (pixel[k] + pixel[l] - 2. * pixel[ipix]);
		
		// ring derivative
		if (ring == Nside)
		{
			k = ring * (ring+1) * 2l + q * ring + j;
			l = k-1;
			if (q == 0 && j == 0)
				l += 4*Nside;
		
			if (pixel[12l*Nside*Nside-1-k] < -1e30 || pixel[12l*Nside*Nside-1-l] < -1e30) return false;
		
			result = (pixel[ipix] - 0.5 * (pixel[12l*Nside*Nside-1-k] + pixel[12l*Nside*Nside-1-l])) * 8. / 3.;
			temp = (0.5 * (pixel[12l*Nside*Nside-1-k] + pixel[12l*Nside*Nside-1-l]));
			w1 = 0.125;
		}
		else
		{
			k = ring * (ring+1) * 2l + q * (ring+1) + j;
			l = k+1;
		
			if (pixel[12l*Nside*Nside-1-k] < -1e30 || pixel[12l*Nside*Nside-1-l] < -1e30) return false;
		
			result = (1. - (j+0.5)/ring) * pixel[12l*Nside*Nside-1-k] + ((j+0.5)/ring) * pixel[12l*Nside*Nside-1-l];
			temp = result;
			w1 = (ring-j-0.5) * (j+0.5) * 0.5 / (ring+1) / (ring+1);
		}
	
		if (ring == 1)
		{
			if (pixel[12l*Nside*Nside-1] < -1e30 || pixel[12l*Nside*Nside-2] < -1e30 || pixel[12l*Nside*Nside-3] < -1e30 || pixel[12l*Nside*Nside-4] < -1e30) return false;
				
			result -= 0.25 * (pixel[12l*Nside*Nside-1] + pixel[12l*Nside*Nside-2] + pixel[12l*Nside*Nside-3] + pixel[12l*Nside*Nside-4]);
			temp += 0.25 * (pixel[12l*Nside*Nside-1] + pixel[12l*Nside*Nside-2] + pixel[12l*Nside*Nside-3] + pixel[12l*Nside*Nside-4]);
			w2 = 0;
		}
		else
		{
			k = (ring-2) * (ring-1) * 2l + q * (ring-1) + j;
			l = k-1;
			if (q == 0 && j == 0)
				l += 4*(ring-1);
			if (q == 3 && j == ring-1)
				k -= 4*(ring-1);
				
			if (pixel[12l*Nside*Nside-1-k] < -1e30 || pixel[12l*Nside*Nside-1-l] < -1e30) return false;
			
			if (ring == Nside)
				result -= pixel[ipix] - ((1. - (j+0.5)/ring) * pixel[12l*Nside*Nside-1-k] + ((j+0.5)/ring) * pixel[12l*Nside*Nside-1-l]);
			else
				result -= (1. - (j+0.5)/ring) * pixel[12l*Nside*Nside-1-k] + ((j+0.5)/ring) * pixel[12l*Nside*Nside-1-l];
			temp += (1. - (j+0.5)/ring) * pixel[12l*Nside*Nside-1-k] + ((j+0.5)/ring) * pixel[12l*Nside*Nside-1-l];
			w2 = (ring-j-0.5) * (j+0.5) * 0.5 / (ring-1) / (ring-1);
		}

		result += (w1 - w2) * temp2;
		result *= (6 * Nside * Nside - 3 * ring * ring) / 8. / ring;
		result += (6 * Nside * Nside - ring * ring) * (temp - 2. * pixel[ipix] - (w1 + w2) * temp2) / 4.;
		
	
		result += 36. * Nside * Nside * Nside * Nside * temp2  / M_PI / M_PI / (6 * Nside * Nside - ring * ring);
	}
	
	return true;
}

