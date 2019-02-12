/**************************************************************/
/*      Finite difference code for Twist Bend Nematics        */
/*         created November 2018, Gareth Alexander            */
/**************************************************************/

#include "TriCubicInterpolator.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <complex>
#include <string>
// for parallelisation
#include <omp.h>
// for the minimiser
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

using namespace std;

#ifndef twistbend_H
#define twistbend_H

// ============================================================

const int Nmax = 300000;       // Number of timesteps
const int stepskip = 5000;   // print pitch and bend every stepskip timesteps
const int Lx = 100;          // System size
const int Ly = 100;          // System size
const int Lz = 50;          // System size


#define BC 0 // periodic (0) or fixed (1) boundary conditions -- currently only fixed along z

// user defined settings
const int    LL=Lx*Ly*Lz;     // total system size
const double    K = 0.04;       // elastic constant
const double dt = 0.65;        // integration timestep
const double thetah = M_PI/11.0;  // heliconical angle
const double qh = 2.0*(2.0*M_PI/Lz);  // heliconical pitch
// lambda = (Kq/2) tan(2thetah)
const double lambda = 0.5*K*qh*tan(2.0*thetah);
// C = K sin^4(thetah)/cos(2thetah)
const double C = K*pow(sin(thetah),4)/cos(2.0*thetah);
const double U = C/9.0; // say

// do we want to read in a file? If so, whats its name?
enum InitialisationType {FROM_FUNCTION,FROM_FILE};
const InitialisationType InitialisationMethod = FROM_FUNCTION;
const string director_filename="vtk_director_100000.vtk";
const string polarisation_filename="vtk_polarisation_100000.vtk";
const int starttime=0;
const char prefix[] = ""; // CHANGE THIS TO A FILE ON YOUR COMPUTER

// ============================================================
// structures which are used in the code
struct parameters
{
	gsl_vector *v,*f,*b;
    likely::TriCubicInterpolator* ucvmag;
};

struct knotpoint
{
    double xcoord;   //position vector x coord
    double ycoord;   //position vector y coord
    double zcoord;   //position vector z coord

    double tx;   //position vector x coord
    double ty;   //position vector y coord
    double tz;   //position vector z coord

    double gradmagbx;   //position vector x coord
    double gradmagby;   //position vector x coord
    double gradmagbz;   //position vector x coord

    double gradmagbperpx;   //position vector x coord
    double gradmagbperpy;   //position vector x coord
    double gradmagbperpz;   //position vector x coord

};
struct knotcurve
{
    std::vector<knotpoint> knotcurve; // the actual data of the curve
    // global data for the knot component
};

/* functions */
void startconfig(int& n ,double* nx, double* ny,double* nz,double* px, double* py,double* pz);
void update(double* nx, double* ny,double* nz,double* px, double* py,double* pz, double* hx, double* hy,double* hz,double* hpx, double* hpy,double* hpz);
void computeBendAndCurlofCirculation(const int n, const double* nx,const double* ny,const double* nz, double* bx, double* by, double* bz, double* bmag, const double* px, const double* py, const double* pz, double* pmag, double* tx, double* ty, double* tz);
void FindBendZeros(double *magb,double* tx,double* ty,double* tz, vector<knotcurve>& knotcurves,double n, gsl_multimin_fminimizer* minimizerstate);
double my_minimisation_function(const gsl_vector* minimum, void* params);
int pt(const int k,const  int l,const  int m);       //convert i,j,k to single index

#endif //twistbend_H
