#ifndef LATFIELD2_MOVE_FUNCTION_HPP
#define LATFIELD2_MOVE_FUNCTION_HPP




/**
 * \addtogroup prartClass
 * @{
 */

using namespace LATfield2;


void move_particles_simple(double dtau,
                               double lat_resolution,
                               part_simple * part,
                               double * ref_dist,
                               part_simple_info partInfo,
                               Field<Real> ** fields,
                               Site * sites,
                               int nfield,
                               double * params,
                               double * outputs,
                               int noutputs){


    //double a;
    //a = 1 + 23;
    for (int l=0;l<3;l++) (*part).pos[l] += dtau*(*part).vel[l];
   
}

void move_particles_round(double dtau,
                               double lat_resolution,
                               part_simple * part,
                               double * ref_dist,
                               part_simple_info partInfo,
                               Field<Real> ** fields,
                               Site * sites,
                               int nfield,
                               double * params,
                               double * outputs,
                               int noutputs){


    //double a;
    //a = 1 + 23;
		double r,r2;
    double v2 = 0;
    double x,y,omg,omgt;

		x = (*part).pos[0]-0.5;
		y = (*part).pos[1]-0.5;
		r2 = x*x+y*y;
		r = sqrt(r2);

    for(int i=0;i<2;i++)
    { v2 += (*part).vel[i] * (*part).vel[i];}

		omg = sqrt(v2/r2);
		omg = 2;
		omgt = acos( x / r );
		if (y<0) omgt = -omgt;
		
		omgt += omg*dtau;
		
    (*part).pos[0] = r*cos(omgt)+0.5; 
    (*part).pos[1] = r*sin(omgt)+0.5; 
}

/**@}*/

#endif
