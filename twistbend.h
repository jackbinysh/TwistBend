/**************************************************************/
/*      Finite difference code for Twist Bend Nematics        */
/*         created November 2018, Gareth Alexander            */
/**************************************************************/

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <complex>
#include <string>
#include <omp.h>

using namespace std;

#ifndef twistbend_H
#define twistbend_H

// ============================================================

const int Nmax = 300000;       // Number of timesteps
const int stepskip = 1000;   // print pitch and bend every stepskip timesteps
const int Lx = 100;          // System size
const int Ly = 100;          // System size
const int Lz = 50;          // System size


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

#define BC 0 // periodic (0) or fixed (1) boundary conditions -- currently only fixed along z

/* functions */
void startconfig(int& n ,double* nx, double* ny,double* nz,double* px, double* py,double* pz);
void update(double* nx, double* ny,double* nz,double* px, double* py,double* pz, double* hx, double* hy,double* hz,double* hpx, double* hpy,double* hpz);
void computeBendAndCurlofCirculation(const int n, const double* nx,const double* ny,const double* nz, double* bx, double* by, double* bz, double* bmag, double* tx, double* ty, double* tz);
int pt(const int k,const  int l,const  int m);       //convert i,j,k to single index

#endif //twistbend_H
