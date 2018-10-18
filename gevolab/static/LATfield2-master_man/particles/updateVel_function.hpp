#ifndef LATFIELD2_UPDATEVEL_FUNCTION_HPP
#define LATFIELD2_UPDATEVEL_FUNCTION_HPP



using namespace LATfield2;

/**
 * \addtogroup prartClass
 * @{
 */

Real updateVel_simple(double dtau,
                          double lat_resolution,
                          part_simple * part,
                          double * ref_dist,
                          part_simple_info partInfo,
                          Field<Real> ** fields,
                          Site * sites,
                          int nfield,
                          double * params,
                          double * outputs,
                          int noutputs)
{

    double v2;

    
    for(int i=0;i<3;i++)
    {
        (*part).vel[i] = (*part).vel[i];
        v2 += (*part).vel[i] * (*part).vel[i];
    }
    
    
    return v2;

}

Real updateVel_round(double dtau,
                          double lat_resolution,
                          part_simple * part,
                          double * ref_dist,
                          part_simple_info partInfo,
                          Field<Real> ** fields,
                          Site * sites,
                          int nfield,
                          double * params,
                          double * outputs,
                          int noutputs)
{

    double v2 = 0;

    for(int i=0;i<2;i++)
    { v2 += (*part).vel[i] * (*part).vel[i];}

    double x = (*part).pos[0]-0.5;
    double y = (*part).pos[1]-0.5;
    double r = sqrt(x*x+y*y);

    (*part).vel[0] = -sqrt(v2)*(y/r);
    (*part).vel[1] = sqrt(v2)*(x/r);
    (*part).vel[2] = 0;
    
    return v2;

}

/**@}*/

#endif
