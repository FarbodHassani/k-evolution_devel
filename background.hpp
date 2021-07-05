//////////////////////////
// background.hpp
//////////////////////////
//
// code components related to background evolution
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: September 2018
//
//////////////////////////

#ifndef BACKGROUND_HEADER
#define BACKGROUND_HEADER

#include <gsl/gsl_integration.h>


//TODO: get it from Class
double FermiDiracIntegrand(double q, void * w)
{
	return q * q * sqrt(q * q + *(double *)w) / (exp(q) + 1.0l);
}

//////////////////////////
// FermiDiracIntegral
//////////////////////////
// Description:
//   computes the integral of the relativistic Fermi-Dirac distribution
//
// Arguments:
//   w          parameter in the F-D distribution, "(m a / kB T)^2"
//
// Returns: value for the integral
//
//////////////////////////

//TODO: get it from Class
double FermiDiracIntegral(double &w)
{
	double result;
	gsl_function f;
	double err;
	size_t n;

	f.function = &FermiDiracIntegrand;
	f.params = &w;

	gsl_integration_qng(&f, 0.0l, 24.0l, 5.0e-7, 1.0e-7, &result, &err, &n);

	return result;
}


//////////////////////////
// bg_ncdm (1)
//////////////////////////
// Description:
//   computes the background model for one ncdm species by integrating the relativistic
//   Fermi-Dirac distribution
//
// Arguments:
//   a          scale factor at which to compute the background model
//   cosmo      structure containing the cosmological parameters
//   p          index of the ncdm species
//
// Returns: value for the background model
//
//////////////////////////

//TODO: get it from Class
double bg_ncdm(const double a, const cosmology cosmo, const int p)
{
	if (p < 0 || p >= cosmo.num_ncdm)
		return 0;
	else
	{
		double w = a * cosmo.m_ncdm[p] / (pow(cosmo.Omega_g * cosmo.h * cosmo.h / C_PLANCK_LAW, 0.25) * cosmo.T_ncdm[p] * C_BOLTZMANN_CST);
		w *= w;

		return FermiDiracIntegral(w) * cosmo.Omega_ncdm[p] * pow(cosmo.Omega_g * cosmo.h * cosmo.h / C_PLANCK_LAW, 0.25) * cosmo.T_ncdm[p] * C_BOLTZMANN_CST / cosmo.m_ncdm[p] / C_FD_NORM / a;
	}
}


//////////////////////////
// bg_ncdm (2)
//////////////////////////
// Description:
//   computes the background model for all ncdm species by integrating the relativistic
//   Fermi-Dirac distribution
//
// Arguments:
//   a          scale factor at which to compute the background model
//   cosmo      structure containing the cosmological parameters
//
// Note:
//   For optimization, the last value of a is stored in a static variable such that
//   multiple calls at the same value of a will not result in multiple integrations
//   being carried out. This assumes that the cosmological model should not change!
//
// Returns: value for the background model
//
//////////////////////////

//TODO: get it from Class
double bg_ncdm(const double a, const cosmology cosmo)
{
	double w;
	static double result = -1.0;
	static double a_prev = -1.0;

	if (a != a_prev)
	{
		result = 0.0;
		a_prev = a;

		for (int p = 0; p < cosmo.num_ncdm; p++)
			result += bg_ncdm(a, cosmo, p);
	}

	return result;
}


//////////////////////////
// Hconf
//////////////////////////
// Description:
//   computes the conformal Hubble rate at given scale factor
//
// Arguments:
//   a          scale factor
//   fourpiG    "4 pi G"
//   H_spline   Class spline with physical H
//   acc        Gsl acc parameter
//   cosmo      structure containing the cosmological parameters
//
// Returns: conformal Hubble rate
//
//////////////////////////

// Hconf normalized to critial density so we have H0^2= 8piG/3
double Hconf(const double a, const double fourpiG,
	#ifdef HAVE_CLASS_BG
	gsl_spline * H_spline, gsl_interp_accel * acc
	#else
	const cosmology cosmo
	#endif
)
{
	#ifdef HAVE_CLASS_BG
	double norm = sqrt(2./3.*fourpiG)/gsl_spline_eval(H_spline, 1., acc);
	double Ha = gsl_spline_eval(H_spline, a, acc);
	// The extra scale factor comes from converting physical to conformal H
	return norm*a*Ha;
	#else
	return sqrt((2. * fourpiG / 3.) * (((cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) / a) + (cosmo.Omega_Lambda * a * a)
	+ (cosmo.Omega_rad / a / a)+ (cosmo.Omega_kessence * pow(a,-3.-3. * cosmo.w_kessence)* a * a)));
	// cout<<"Omega_rad: "<<cosmo.Omega_rad<<" cosmo.Omega_Lambda"<<cosmo.Omega_Lambda<<endl;
	#endif
}

// Here the normalization factor is not \rho_crit=1, it is what it should be in th enormal unit.
// So Omega_m is the matter density at arbitrary redshift and is not normalized, since we did not use Hconf in the fomrula
// While Hconf is normalized to critical density 1 so H^2/H_0^2= H^2/(8piG/3) which is used in the last formula.
//TODO: check where this is used and take it from Class
double Omega_m(const double a, const cosmology cosmo) { return cosmo.Omega_m / (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo) + cosmo.Omega_kessence * pow(a,-3.-3. * cosmo.w_kessence)* a * a * a + cosmo.Omega_Lambda * a * a * a + cosmo.Omega_rad / a); }
//
//TODO: check where this is used and take it from Class
double Omega_rad(const double a, const cosmology cosmo) { return (cosmo.Omega_rad + (bg_ncdm(a, cosmo) + cosmo.Omega_cdm + cosmo.Omega_b - cosmo.Omega_m) * a) / ((cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) * a + cosmo.Omega_kessence * pow(a,-3.-3. * cosmo.w_kessence)* a * a * a * a + cosmo.Omega_Lambda * a * a * a * a + cosmo.Omega_rad); }

//Here Omega_Lambda is just Lambda
//TODO: check where this is used and take it from Class
double Omega_Lambda(const double a, const cosmology cosmo) { return cosmo.Omega_Lambda / ((cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) / a / a / a + cosmo.Omega_Lambda + cosmo.Omega_kessence * pow(a,-3.-3. * cosmo.w_kessence) + cosmo.Omega_rad / a / a / a / a);}

//TODO: remove this
double Hconf_class(const double a, const cosmology cosmo)
{
  double H0_class=100*cosmo.h/(C_SPEED_OF_LIGHT*100.);
 //0.00022593979933110373; // H0 in unit of 1/Mpc H0=100h/c;
	return H0_class * sqrt( ((cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) / a) + (cosmo.Omega_Lambda * a * a)
	+ (cosmo.Omega_rad / a / a)+ (cosmo.Omega_kessence * pow(a,-3.-3. * cosmo.w_kessence)* a * a) );
}

//////////////////////////
// Hconf_prime
//////////////////////////
// Description:
//   computes the conformal Hubble rate derivative at given scale factor
//
// Arguments:
//   a          scale factor
//   fourpiG    "4 pi G"
//   H_spline   Class spline with physical H
//   acc        Gsl acc parameter
//   cosmo      structure containing the cosmological parameters
//
// Returns: conformal Hubble rate
//
//////////////////////////
// Hconf normalized to critial density so we have H0^2= 8piG/3
double Hconf_prime(const double a, const double fourpiG,
	#ifdef HAVE_CLASS_BG
	gsl_spline * H_spline, gsl_interp_accel * acc
	#else
	const cosmology cosmo
	#endif
)
{
	#ifdef HAVE_CLASS_BG
	double norm = sqrt(2./3.*fourpiG)/gsl_spline_eval(H_spline, 1., acc);
	double Hc = Hconf(a, fourpiG,
		#ifdef HAVE_CLASS_BG
		H_spline, acc
		#else
		cosmo
		#endif
		);
	// dHc/da = H + a*d(H)/da
	double dHcda = Hc/a + a*norm*gsl_spline_eval_deriv(H_spline, a, acc);
	// dHconf/dtau = a*Hc*dHc/da
	return a*Hc*dHcda;
	#else
	return (2. * fourpiG / (6. * a * a)) * (  -(cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) * a + 2. * cosmo.Omega_Lambda * a *  a * a * a - 2. * cosmo.Omega_rad - (1. + 3. * cosmo.w_kessence) * cosmo.Omega_kessence * pow( a, 1.-3.* cosmo.w_kessence));
	#endif
}

//////////////////////////
// rungekutta4bg
//////////////////////////
// Description:
//   integrates the Friedmann equation for the background model using a fourth-order
//   Runge-Kutta method
//
// Arguments:
//   a          scale factor (will be advanced by dtau)
//   fourpiG    "4 pi G"
//   H_spline   Class spline with physical H
//   acc        Gsl acc parameter
//   cosmo      structure containing the cosmological parameters
//   dtau       time step by which the scale factor should be advanced
//
// Returns:
//
//////////////////////////

void rungekutta4bg(double &a, const double fourpiG,
	#ifdef HAVE_CLASS_BG
	gsl_spline * H_spline, gsl_interp_accel * acc,
	#else
	const cosmology cosmo,
	#endif
	const double dtau)
{
	double k1a, k2a, k3a, k4a;

	k1a = a * Hconf(a, fourpiG,//TODO_EB
		#ifdef HAVE_CLASS_BG
			H_spline, acc
		#else
			cosmo
		#endif
		);
	k2a = (a + k1a * dtau / 2.) * Hconf(a + k1a * dtau / 2., fourpiG,//TODO_EB
		#ifdef HAVE_CLASS_BG
			H_spline, acc
		#else
			cosmo
		#endif
		);
	k3a = (a + k2a * dtau / 2.) * Hconf(a + k2a * dtau / 2., fourpiG,//TODO_EB
		#ifdef HAVE_CLASS_BG
			H_spline, acc
		#else
			cosmo
		#endif
		);
	k4a = (a + k3a * dtau) * Hconf(a + k3a * dtau, fourpiG,//TODO_EB
		#ifdef HAVE_CLASS_BG
			H_spline, acc
		#else
			cosmo
		#endif
		);

	a += dtau * (k1a + 2. * k2a + 2. * k3a + k4a) / 6.;
}


#ifndef HAVE_CLASS_BG
double particleHorizonIntegrand(double sqrta, void * cosmo)
{
	double Hc = Hconf(sqrta*sqrta, 1., *(cosmology *)cosmo);
	return 2. / (sqrta * Hc);
}
#endif

//////////////////////////
// particleHorizon
//////////////////////////
// Description:
//   computes the particle horizon (tau) at given scale factor
//
// Arguments:
//   a          scale factor
//   fourpiG    "4 pi G"
//   cosmo      structure containing the cosmological parameters
//
// Returns: particle horizon (tau)
//
//////////////////////////

double particleHorizon(const double a, const double fourpiG,
	#ifdef HAVE_CLASS_BG
	const double H_spline_0, background & class_background
	#else
	cosmology & cosmo
	#endif
)
{
	#ifdef HAVE_CLASS_BG
	double tau;
	background_tau_of_z(&class_background, 1./a - 1., &tau);
	return tau*H_spline_0/sqrt(2./3.*fourpiG);
	#else
	double result;
	gsl_function f;
	double err;
	size_t n;

	f.function = &particleHorizonIntegrand;
	f.params = &cosmo;

	gsl_integration_qng(&f, sqrt(a) * 1.0e-7, sqrt(a), 5.0e-7, 1.0e-7, &result, &err, &n);

	return result / sqrt(fourpiG);
	#endif
}

#endif
