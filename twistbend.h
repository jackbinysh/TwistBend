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

const int Nmax = 50000;       // Number of timesteps
const int vtkstepskip = 1000;   // print pitch and bend every stepskip timesteps
const int curvestepskip = 1000;   // print pitch and bend every stepskip timesteps
const int curvestarttime = 2000;   // print pitch and bend every stepskip timesteps
const int Lx = 130;          // System size
const int Ly = 130;          // System size
const int Lz = 130;          // System size

#define BC 1 // periodic (0) or fixed (1) boundary conditions -- currently only fixed along z

// user defined settings
const int    LL=Lx*Ly*Lz;     // total system size
const double    K = 0.04;       // elastic constant
const double dt = 0.65;        // integration timestep
const double thetah = M_PI/11.0;  // heliconical angle
const double qh = 1.1*(2.0*M_PI/Lz);  // heliconical pitch
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
	gsl_vector *mypt,*t,*e1,*e2;
    likely::TriCubicInterpolator* ucvmag;
    double sphereradius;
};

struct knotpoint
{
    // point
    double xcoord;   //position vector x coord
    double ycoord;   //position vector y coord
    double zcoord;   //position vector z coord

    // tangent vector
    double tx;   
    double ty;   
    double tz;   

    // the solid angle framing
    double omegax;
    double omegay;
    double omegaz;

    // the bend vector (mainly for the pushoff)
    double bx;
    double by;
    double bz;

    // some other framing
    double e1x;   //position vector x coord
    double e1y;   //position vector x coord
    double e1z;   //position vector x coord

    double e2x;   //position vector x coord
    double e2y;   //position vector x coord
    double e2z;   //position vector x coord

    double length;

};
struct knotcurve
{
    std::vector<knotpoint> knotcurve; // the actual data of the curve
    // global data for the knot component
    // is this an open curve, i.e one which terminates at the boundaries, or not?
    bool closed;
};

struct Link
{
    std::vector<knotcurve> Components;
}; 

/* functions */
void startconfig(int& n ,double* nx, double* ny,double* nz,double* px, double* py,double* pz);
void update(double* nx, double* ny,double* nz,double* px, double* py,double* pz, double* hx, double* hy,double* hz,double* hpx, double* hpy,double* hpz);
void computeBendTwistEnergyOrientation(const int n, const double* nx,const double* ny,const double* nz, double* bx, double* by, double* bz,const double* px,const double* py, const double* pz, double* bmag, double* twist,double* FreeEnergy, double* tx, double* ty, double* tz);
void FindBendZeros(Link& Curve, Link& PushOffCurve, double* bx,double* by,double* bz, double *magb,double* tx,double* ty,double* tz, bool* mask);
double my_minimisation_function(const gsl_vector* minimum, void* params);
int pt(const int k,const  int l,const  int m);       //convert i,j,k to single index
int incp(int i, int p, int N);    //increment i with p for periodic boundary
int mod(int i, int N);   //my own mod fn
void CurveSmoothing(Link& Curve, int filterlength);
bool KnotpointInMask(const knotpoint knotpoint, const bool* mask);
void setupmask(bool* mask);

#endif //twistbend_H
