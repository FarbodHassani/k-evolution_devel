#include <stdlib.h>
#include <iostream>
#include <math.h>

using namespace std;

int main(int argc, char **argv)
{
	double direction[3];// = {0.735, 0.485, 1.};
	int pos[3];
	double x[3];
	int wrap[3];// = {0, 0, 0};
	double vertex[3];// = {0., 0., 0.};
	double temp;
	
	if (argc < 5) return -1;
	
	sscanf(argv[1], " %lf ", direction);
	sscanf(argv[2], " %lf ", direction+1);
	sscanf(argv[3], " %lf ", direction+2);
	sscanf(argv[4], " %lf ", vertex);
	
	vertex[1] = vertex[0];
	vertex[2] = vertex[0];
	wrap[0] = (int) vertex[0];
	wrap[1] = (int) vertex[1];
	wrap[2] = (int) vertex[2];
	
	temp = sqrt(direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);
	
	direction[0] /= temp;
	direction[1] /= temp;
	direction[2] /= temp;
	
	for (int i = 1; i < 20000; i++)
	{
		pos[0] = (int) floor(i * direction[0]);
		pos[1] = (int) floor(i * direction[1]);
		pos[2] = (int) floor(i * direction[2]);
		
		if (fabs(pos[0]/direction[0] - pos[1]/direction[1]) < 1e-8 && fabs(pos[1]/direction[1] - pos[2]/direction[2]) < 1e-8)
		{
			cout << " trying " << pos[0] << ", " << pos[1] << ", " << pos[2] << ": ";
			x[0] = (double) pos[0] / 7680.;
			x[1] = (double) pos[1] / 7680.;
			x[2] = (double) pos[2] / 7680.;
			
			temp = sqrt((x[0]-(vertex[0]-wrap[0]))*(x[0]-(vertex[0]-wrap[0])) + (x[1]-(vertex[1]-wrap[1]))*(x[1]-(vertex[1]-wrap[1])) + (x[2]-(vertex[2]-wrap[2]))*(x[2]-(vertex[2]-wrap[2])));
			
			temp = acos(((x[0]-(vertex[0]-wrap[0]))*direction[0] + (x[1]-(vertex[1]-wrap[1]))*direction[1] + (x[2]-(vertex[2]-wrap[2]))*direction[2])/temp);
			
			cout << temp << endl;
		}
	}

	return 1;
}

