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

using namespace std;

// ============================================================

const int Nmax = 8000;       // Number of timesteps
const int stepskip = 1000;   // print pitch and bend every stepskip timesteps
const int FePrintInt = 100;  // print free energy every FePrintInt timesteps
const int Lx = 80;          // System size
const int Ly = 80;          // System size
const int Lz = 40;          // System size

// ============================================================

#define BC 0 // periodic (0) or fixed (1) boundary conditions -- currently only fixed along z
#define CYLINDER 0 // default (0); cylinder (1) 

/* functions */
void initialise(void);
void startconfig(void);
void update(void);
void writeVTKfiles(void);
void writeBENDfiles(void);
void cholestericPitch(void);
void writeTwist(void);
// modified numerical recipes routine
#define n 3     
void jacobi(double (*a)[n], double d[], double (*v)[n], int *nrot);
#undef n

/* global variables */
int n,LL;
double K,q0,C,lambda,theta,qh,U,Gamma,dt;
double *nx,*ny,*nz,*px,*py,*pz;
double *hx,*hy,*hz,*hpx,*hpy,*hpz;
ofstream output,output2,outputFE;

char free_energy[200],vtk_director[200],vtk_polarisation[200],vtk_TWIST[200],vtk_BEND[200],pitch_defect[200],pitch_p[200];
char prefix[] = ""; // CHANGE THIS TO A FILE ON YOUR COMPUTER

// ============================================================

int main(int argc, char** argv) 
{
  LL=Lx*Ly*Lz;     // total system size
  K = 0.04;       // elastic constant
  q0 = 0.0*M_PI/Lx;    // chirality
  Gamma = 0.65;     // relaxation constant
  dt = 1.0;        // integration timestep
  theta = M_PI/6.0;  // heliconical angle
  qh = 1.0*(2.0*M_PI/Lz);  // heliconical pitch
  // lambda = (Kq/2) tan(2theta)
  lambda = 0.5*K*qh*tan(2.0*theta);
  // C = K sin^4(theta)/cos(2theta)
  C = K*pow(sin(theta),4)/cos(2.0*theta);
  U = C/9.0; // say

  int step=0;

  cout << "initialising" << endl;
  initialise();
  startconfig();
  cout << "starting simulation" << endl;

  //#if CYLINDER // open file for the free energy
  //  sprintf(free_energy,"%sfree_energy_%dx%dx%d.dat",prefix,Lx,Ly,Lz); 
  //  outputFE.open(free_energy);
  //  outputFE.precision(12);
  //#endif
  
  for (n=0; n<=Nmax; n++) {

    if (n%1000==0) {
      cout << "timestep " << n << endl;
    }

    if (step==stepskip || step==0) {
      cout << "writing VTK files at timestep " << n << endl;
      writeVTKfiles();       // output VTK files for use with ParaView
      //      cholestericPitch();    // output VTK files for use with ParaView
      writeBENDfiles();      // output VTK files for use with ParaView
      //      writeTwist();      // output VTK files for use with ParaView
      step=0;
    }

    step++;
    update();
  }
  //#if CYLINDER  // close file for free energy
  //  outputFE.close();
  //#endif
} // end main


/**********************************************************************/
void initialise(void)
{
  nx=new double[LL];
  ny=new double[LL];
  nz=new double[LL];

  px=new double[LL];
  py=new double[LL];
  pz=new double[LL];

  hx=new double[LL];
  hy=new double[LL];
  hz=new double[LL];

  hpx=new double[LL];
  hpy=new double[LL];
  hpz=new double[LL];
}

/**********************************************************************/
void startconfig(void)
{
  int j,k,l,m;
  int HELICONICAL,CYL,SKYRMION,HOPF,SQUARE;
  double k0,l0,m0,rr,RR,norm;     // variables for a droplet

  // define texture to simulate -- 0 or 1
  HELICONICAL = 1; // standard heliconical texture
  CYL = 0; // escape in a cylinder
  SKYRMION = 0; // hexagonal skyrmion lattice
  HOPF = 0; // Hopf texture
  SQUARE = 0;

  // initial configuration
  k=l=m=0;
  // variables for a droplet
  k0 = Lx/2.0 - 0.5; l0 = Ly/2.0 - 0.5; m0 = Lz/2.0 - 0.5; // offset from the lattice points
  RR = 0.98*k0; //q0 = 2.0*M_PI/(1.8*RR); // scale
  
  // standard heliconical texture
  if (HELICONICAL == 1) {
    for (j=0; j<LL; j++) {
      // director field
      nx[j] = sin(theta)*cos(qh*m);
      ny[j] = sin(theta)*sin(qh*m);
      nz[j] = cos(theta);
      // normalise 
      //      norm = sqrt(nx[j]*nx[j]+ny[j]*ny[j]+nz[j]*nz[j]); 
      //      nx[j] /= norm;
      //      ny[j] /= norm;
      //      nz[j] /= norm;
      // polarisation
      px[j] = -sin(qh*m);
      py[j] = cos(qh*m);
      pz[j] = 0.0;
      // normalise 
      //      norm = sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]); 
      //      px[j] /= norm;
      //      py[j] /= norm;
      //      pz[j] /= norm;
      // perversion -- try this
      //      if ((m>0.25*Lz)&&(m<0.75*Lz)) {
//      if (m>0.5*Lz) {
//	ny[j] *= -1.0;
//	py[j] *= -1.0;
	// test
	//	nx[j] *= -1.0;
	//	nz[j] *= -1.0;
     // }
      // deal with the periodic boundaries            
      k++;
      if (k==Lx) {l++; k=0;}
      if (l==Ly) {m++; l=0;}
    }
  }

  // escape in a cylinder
  if (CYL == 1) {
    for (j=0; j<LL; j++) {
      // director field
      nx[j] = sin(theta)*cos(qh*m);
      ny[j] = sin(theta)*sin(qh*m);
      nz[j] = cos(theta);
      // polarisation
      px[j] = -sin(qh*m);
      py[j] = cos(qh*m);
      pz[j] = 0.0;
      // normalise -- not strictly needed
      norm = sqrt(nx[j]*nx[j]+ny[j]*ny[j]+nz[j]*nz[j]); 
      nx[j] /= norm;
      ny[j] /= norm;
      nz[j] /= norm;
      //      norm = sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]); 
      //      px[j] /= norm;
      //      py[j] /= norm;
      //      pz[j] /= norm;
      // deal with the periodic boundaries            
      k++;
      if (k==Lx) {l++; k=0;}
      if (l==Ly) {m++; l=0;}
    }
  }

  // hexagonal Skyrmion lattice
  if (SKYRMION == 1) {
    qh = 4.0*M_PI/97.0; // size of an hexagonal cell -- 56 x 97
    lambda = 0.5*K*qh*tan(2.0*theta);
    for (j=0; j<LL; j++) {
      // director field
      nx[j] = -(sin(qh*(1.0*l-l0))-0.5*sin(qh*(sqrt(3.0)*(1.0*k-k0)-(1.0*l-l0))/2.0)+0.5*sin(qh*(sqrt(3.0)*(1.0*k-k0)+(1.0*l-l0))/2.0));
      ny[j] = -(-sqrt(3.0)*0.5*sin(qh*(sqrt(3.0)*(1.0*k-k0)-(1.0*l-l0))/2.0)-sqrt(3.0)*0.5*sin(qh*(sqrt(3.0)*(1.0*k-k0)+(1.0*l-l0))/2.0));
      nz[j] = -(cos(qh*(1.0*l-l0))+cos(qh*(sqrt(3.0)*(1.0*k-k0)-(1.0*l-l0))/2.0)+cos(qh*(sqrt(3.0)*(1.0*k-k0)+(1.0*l-l0))/2.0));
      // normalise 
      norm = sqrt(nx[j]*nx[j]+ny[j]*ny[j]+nz[j]*nz[j]); 
      nx[j] /= norm;
      ny[j] /= norm;
      nz[j] /= norm;
      // polarisation -- try this ??
      px[j] = -(qh*ny[j]*cos(qh*(1.0*l-l0)) + qh*0.25*(sqrt(3.0)*nx[j]+ny[j])*cos(qh*(sqrt(3.0)*(1.0*k-k0)+(1.0*l-l0))/2.0) - 0.25*qh*(sqrt(3.0)*nx[j]-ny[j])*cos(qh*(sqrt(3.0)*(1.0*k-k0)-(1.0*l-l0))/2.0));
      py[j] = -(-qh*sqrt(3.0)*0.25*(sqrt(3.0)*nx[j]+ny[j])*cos(qh*(sqrt(3.0)*(1.0*k-k0)+(1.0*l-l0))/2.0) - qh*0.25*sqrt(3.0)*(sqrt(3.0)*nx[j]-ny[j])*cos(qh*(sqrt(3.0)*(1.0*k-k0)-(1.0*l-l0))/2.0));
      pz[j] = -(-qh*ny[j]*sin(qh*(1.0*l-l0)) - qh*0.5*(sqrt(3.0)*nx[j]+ny[j])*sin(qh*(sqrt(3.0)*(1.0*k-k0)+(1.0*l-l0))/2.0) - qh*0.5*(sqrt(3.0)*nx[j]-ny[j])*sin(qh*(sqrt(3.0)*(1.0*k-k0)-(1.0*l-l0))/2.0));
      // normalise
      norm = sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]); 
      px[j] /= norm;
      py[j] /= norm;
      pz[j] /= norm;
      // deal with the periodic boundaries            
      k++;
      if (k==Lx) {l++; k=0;}
      if (l==Ly) {m++; l=0;}
    }
    qh *= 0.72; // simple hack
  }

  // hopf texture      
  if (HOPF == 1) { 
    int xup,xdwn,yup,ydwn,zup,zdwn;
    double rho,theta,phi;    // variables for hopf texture 
    RR = 16.0; qh = M_PI/RR;  // scale 

    for (j=0; j<LL; j++) {
      nx[j] = 0.0; ny[j] = 0.0; nz[j] = 1.0;
      // cheap trick
      px[j] = 0.0;
      py[j] = 0.0;
      pz[j] = 0.0;
      // back to Hopf
      rr = sqrt((k-k0)*(k-k0)+(l-l0)*(l-l0));
      rho = sqrt((rr-RR)*(rr-RR)+(m-m0)*(m-m0));
      
      if (rho < RR) {
        theta = atan2(m-m0,rr-RR);
        if (theta < 0.0) {theta+=2.0*M_PI;} // put in the range [0,2pi)
        phi = atan2(l-l0,k-k0);
        if (phi < 0.0) {phi+=2.0*M_PI;} // put in the range [0,2pi)
       
	// cross-section looks like double twist and an escaped -1
	// for -theta these are ordered -1 inside and double twist outside; for +theta it is the other way around  
	// Hopf invariant -1
	//	nx[j] = -sin(M_PI*rho/RR)*sin(phi-theta); // could also use -theta; they seem to be equivalent in energy
	//	ny[j] = sin(M_PI*rho/RR)*cos(phi-theta); // could also use -theta
	// Hopf invariant +1
	nx[j] = -sin(M_PI*rho/RR)*sin(phi+theta); // could also use -theta; they seem to be equivalent in energy
	ny[j] = sin(M_PI*rho/RR)*cos(phi+theta); // could also use -theta
	nz[j] = -cos(M_PI*rho/RR);
      }

#if BC // normal anchoring along z
      if (m==0) {nx[j] = 0.0; ny[j] = 0.0; nz[j] = 1.0;}
      if (m==Lz-1) {nx[j] = 0.0; ny[j] = 0.0; nz[j] = 1.0;}
#endif
           
      // deal with the periodic boundaries            
      k++;
      if (k==Lx) {l++; k=0;}
      if (l==Ly) {m++; l=0;}
    }
    // set the polarisation to the bend of the director 
    k=l=m=0;
    for (j=0; j<LL; j++) {
      // hack here ??
      px[j] = 0.0; py[j] = 0.0; pz[j] = 1.0;
      rr = sqrt((k-k0)*(k-k0)+(l-l0)*(l-l0));
      rho = sqrt((rr-RR)*(rr-RR)+(m-m0)*(m-m0));
      
      // define neighbouring nodes
      xup=j+1; xdwn=j-1; yup=j+Lx; ydwn=j-Lx; zup=j+Lx*Ly; zdwn=j-Lx*Ly;
      // correct for periodic boundaries
      if (k==0) {xdwn+=Lx;}
      if (k==Lx-1) {xup-=Lx;}
      if (l==0) {ydwn+=Lx*Ly;}
      if (l==Ly-1) {yup-=Lx*Ly;}
      if (m==0) {zdwn+=LL;}
      if (m==Lz-1) {zup-=LL;}

      if (rho < RR) {
	px[j] = nx[j]*(nx[xup]-nx[xdwn])+ny[j]*(nx[yup]-nx[ydwn])+nz[j]*(nx[zup]-nx[zdwn]);
	py[j] = nx[j]*(ny[xup]-ny[xdwn])+ny[j]*(ny[yup]-ny[ydwn])+nz[j]*(ny[zup]-ny[zdwn]);
	pz[j] = nx[j]*(nz[xup]-nz[xdwn])+ny[j]*(nz[yup]-nz[ydwn])+nz[j]*(nz[zup]-nz[zdwn]);
	// normalise
	norm = sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]);
	px[j] /= norm; py[j] /= norm; pz[j] /= norm; 
      }
      // deal with the periodic boundaries            
      k++;
      if (k==Lx) {l++; k=0;}
      if (l==Ly) {m++; l=0;}
    }
  } // end hopf

  // square
  if (SQUARE == 1) {
    double norm,r,R,phi;
    R = 0.8*Lx/2.0;

    for (j=0; j<LL; j++) {
      // default is heliconical texture
      // director field
      nx[j] = -sin(theta)*cos(qh*m);
      ny[j] = -sin(theta)*sin(qh*m);
      nz[j] = -cos(theta);
      // polarisation
      px[j] = -sin(qh*m);
      py[j] = cos(qh*m);
      pz[j] = 0.0;
      // replace a cylinder with a twisted Skyrmion -- testing!
      r = sqrt((k-k0)*(k-k0)+(l-l0)*(l-l0));
      if (r<R) {
	phi = atan2(l-l0,k-k0);
	// degree +1
	//	nx[j] = sin(M_PI*r/R)*sin(phi);
	//	ny[j] = -sin(M_PI*r/R)*cos(phi);
	// twisted 
	nx[j] = sin(M_PI*r/R)*sin(phi+2.0*M_PI*m/Lz);
	ny[j] = -sin(M_PI*r/R)*cos(phi+2.0*M_PI*m/Lz);
	// degree -1
	//	nx[j] = -sin(M_PI*r/R)*sin(phi);
	//	ny[j] = -sin(M_PI*r/R)*cos(phi);
	nz[j] = cos(M_PI*r/R);
	// polarisation
	px[j] = cos(M_PI*r/R)*sin(2.0*M_PI*m/Lz)*sin(phi+2.0*M_PI*m/Lz)+((2.0*R/Lz)*cos(M_PI*r/R)-(R/(M_PI*r))*sin(M_PI*r/R)*cos(2.0*M_PI*m/Lz))*cos(phi+2.0*M_PI*m/Lz);
	py[j] = -cos(M_PI*r/R)*sin(2.0*M_PI*m/Lz)*cos(phi+2.0*M_PI*m/Lz)+((2.0*R/Lz)*cos(M_PI*r/R)-(R/(M_PI*r))*sin(M_PI*r/R)*cos(2.0*M_PI*m/Lz))*sin(phi+2.0*M_PI*m/Lz);
	pz[j] = -sin(M_PI*r/R)*sin(2.0*M_PI*m/Lz);
      }
      /*
      // normalise 
      norm = sqrt(nx[j]*nx[j]+ny[j]*ny[j]+nz[j]*nz[j]);
      nx[j] /= norm;
      ny[j] /= norm;
      nz[j] /= norm;
      */
      // deal with the periodic boundaries            
      k++;
      if (k==Lx) {l++; k=0;}
      if (l==Ly) {m++; l=0;}
    }
  } // end square

} // end startconfig

/**********************************************************************/
void update(void)
{
  int j,k,l,m,jj;
  int xup,xdwn,yup,ydwn,zup,zdwn;
  double Dxnx,Dynx,Dznx,Dxxnx,Dyynx,Dzznx;
  double Dxny,Dyny,Dzny,Dxxny,Dyyny,Dzzny;
  double Dxnz,Dynz,Dznz,Dxxnz,Dyynz,Dzznz;
  double hdotn,sqrtndotn;
  double Dxpx,Dypx,Dzpx,Dxxpx,Dyypx,Dzzpx;
  double Dxpy,Dypy,Dzpy,Dxxpy,Dyypy,Dzzpy;
  double Dxpz,Dypz,Dzpz,Dxxpz,Dyypz,Dzzpz;
  double hpdotp;

  k=l=m=0;

#if CYLINDER
  double k0,l0,m0,rr,RR,energy;
  k0 = Lx/2.0 - 0.5; l0 = Ly/2.0 - 0.5; m0 = Lz/2.0 - 0.5; // offset from the lattice points
  RR = 0.98*k0; // NEEDS TO MATCH WITH STARTCONFIG
  energy = 0.0;
#endif

  /* Calculate derivatives, molecular field and the energy */
  for (j=0; j<LL; j++) {
    // define neighbouring nodes
    xup=j+1; xdwn=j-1; yup=j+Lx; ydwn=j-Lx; zup=j+Lx*Ly; zdwn=j-Lx*Ly;
    // correct for periodic boundaries
    if (k==0) {xdwn+=Lx;}
    if (k==Lx-1) {xup-=Lx;}
    if (l==0) {ydwn+=Lx*Ly;}
    if (l==Ly-1) {yup-=Lx*Ly;}
    if (m==0) {zdwn+=LL;}
    if (m==Lz-1) {zup-=LL;}

#if BC // cheap fix for the Dirichlet boundary conditions
    if (m==0) {xup=j; xdwn=j; yup=j; ydwn=j; zup=j; zdwn=j;}
    if (m==Lz-1) {xup=j; xdwn=j; yup=j; ydwn=j; zup=j; zdwn=j;}
#endif

#if CYLINDER // cheap fix for droplet Dirichlet boundary conditions
    rr = sqrt((k-k0)*(k-k0)+(l-l0)*(l-l0));
    if (rr>=RR) {xup=j; xdwn=j; yup=j; ydwn=j; zup=j; zdwn=j;}
#endif
        
    // calculate first order derivatives
    Dxnx = (nx[xup]-nx[xdwn])/2.0;
    Dynx = (nx[yup]-nx[ydwn])/2.0;
    Dznx = (nx[zup]-nx[zdwn])/2.0;

    Dxny = (ny[xup]-ny[xdwn])/2.0;
    Dyny = (ny[yup]-ny[ydwn])/2.0;
    Dzny = (ny[zup]-ny[zdwn])/2.0;

    Dxnz = (nz[xup]-nz[xdwn])/2.0;
    Dynz = (nz[yup]-nz[ydwn])/2.0;
    Dznz = (nz[zup]-nz[zdwn])/2.0;

    Dxpx = (px[xup]-px[xdwn])/2.0;
    Dypx = (px[yup]-px[ydwn])/2.0;
    Dzpx = (px[zup]-px[zdwn])/2.0;

    Dxpy = (py[xup]-py[xdwn])/2.0;
    Dypy = (py[yup]-py[ydwn])/2.0;
    Dzpy = (py[zup]-py[zdwn])/2.0;

    Dxpz = (pz[xup]-pz[xdwn])/2.0;
    Dypz = (pz[yup]-pz[ydwn])/2.0;
    Dzpz = (pz[zup]-pz[zdwn])/2.0;

    // calculate second order derivatives
    Dxxnx = nx[xup]-2.0*nx[j]+nx[xdwn];
    Dyynx = nx[yup]-2.0*nx[j]+nx[ydwn];
    Dzznx = nx[zup]-2.0*nx[j]+nx[zdwn];

    Dxxny = ny[xup]-2.0*ny[j]+ny[xdwn];
    Dyyny = ny[yup]-2.0*ny[j]+ny[ydwn];
    Dzzny = ny[zup]-2.0*ny[j]+ny[zdwn];

    Dxxnz = nz[xup]-2.0*nz[j]+nz[xdwn];
    Dyynz = nz[yup]-2.0*nz[j]+nz[ydwn];
    Dzznz = nz[zup]-2.0*nz[j]+nz[zdwn];
     
    Dxxpx = px[xup]-2.0*px[j]+px[xdwn];
    Dyypx = px[yup]-2.0*px[j]+px[ydwn];
    Dzzpx = px[zup]-2.0*px[j]+px[zdwn];

    Dxxpy = py[xup]-2.0*py[j]+py[xdwn];
    Dyypy = py[yup]-2.0*py[j]+py[ydwn];
    Dzzpy = py[zup]-2.0*py[j]+py[zdwn];

    Dxxpz = pz[xup]-2.0*pz[j]+pz[xdwn];
    Dyypz = pz[yup]-2.0*pz[j]+pz[ydwn];
    Dzzpz = pz[zup]-2.0*pz[j]+pz[zdwn];
     
    // calculate molecular field
    hx[j] = K*(Dxxnx+Dyynx+Dzznx) - 2.0*K*q0*(Dynz-Dzny) + lambda*(px[j]*Dxnx+py[j]*Dxny+pz[j]*Dxnz - px[j]*(Dxnx+Dyny+Dznz) - (nx[j]*Dxpx+ny[j]*Dypx+nz[j]*Dzpx));
    hy[j] = K*(Dxxny+Dyyny+Dzzny) - 2.0*K*q0*(Dznx-Dxnz) + lambda*(px[j]*Dynx+py[j]*Dyny+pz[j]*Dynz - py[j]*(Dxnx+Dyny+Dznz) - (nx[j]*Dxpy+ny[j]*Dypy+nz[j]*Dzpy));
    hz[j] = K*(Dxxnz+Dyynz+Dzznz) - 2.0*K*q0*(Dxny-Dynx) + lambda*(px[j]*Dznx+py[j]*Dzny+pz[j]*Dznz - pz[j]*(Dxnx+Dyny+Dznz) - (nx[j]*Dxpz+ny[j]*Dypz+nz[j]*Dzpz));

    hdotn = nx[j]*hx[j] + ny[j]*hy[j] + nz[j]*hz[j];
    hx[j] -= nx[j]*hdotn;
    hy[j] -= ny[j]*hdotn;
    hz[j] -= nz[j]*hdotn;

    // molecular field for polarisation
    hpx[j] = C*(Dxxpx+Dyypx+Dzzpx) + U*(1.0-px[j]*px[j]-py[j]*py[j]-pz[j]*pz[j])*px[j] + lambda*(nx[j]*Dxnx+ny[j]*Dynx+nz[j]*Dznx); 
    hpy[j] = C*(Dxxpy+Dyypy+Dzzpy) + U*(1.0-px[j]*px[j]-py[j]*py[j]-pz[j]*pz[j])*py[j] + lambda*(nx[j]*Dxny+ny[j]*Dyny+nz[j]*Dzny); 
    hpz[j] = C*(Dxxpz+Dyypz+Dzzpz) + U*(1.0-px[j]*px[j]-py[j]*py[j]-pz[j]*pz[j])*pz[j] + lambda*(nx[j]*Dxnz+ny[j]*Dynz+nz[j]*Dznz); 

#if BC // Dirichlet boundary conditions along z
    if (m==0) {hx[j] = 0.0; hy[j] = 0.0; hz[j] = 0.0;}
    if (m==Lz-1) {hx[j] = 0.0; hy[j] = 0.0; hz[j] = 0.0;}
    if (m==0) {hpx[j] = 0.0; hpy[j] = 0.0; hpz[j] = 0.0;}
    if (m==Lz-1) {hpx[j] = 0.0; hpy[j] = 0.0; hpz[j] = 0.0;}
#endif

#if CYLINDER 
    rr = sqrt((k-k0)*(k-k0)+(l-l0)*(l-l0));
    // Dirichlet boundary conditions for droplet
    if (rr>=RR) {hx[j] = 0.0; hy[j] = 0.0; hz[j] = 0.0;}
    // energy calculation -- one elastic constant -- only for the droplet at present
    if (rr<RR) {
      energy += 0.5*K*(Dxnx*Dxnx+Dxny*Dxny+Dxnz*Dxnz+Dynx*Dynx+Dyny*Dyny+Dynz*Dynz+Dznx*Dznx+Dzny*Dzny+Dznz*Dznz+q0*q0) 
	+ K*q0*(nx[j]*(Dynz-Dzny)+ny[j]*(Dznx-Dxnz)+nz[j]*(Dxny-Dynx));
    }
#endif

    // keep track of boundaries
    k++;
    if (k==Lx) {l++; k=0;}
    if (l==Ly) {m++; l=0;}
  }

#if CYLINDER 
  if ((n % FePrintInt == 0) && (n>0)) {
    outputFE << n << " " << energy << endl;
  }
#endif 

  // do the update -- first order Euler
  for (j=0; j<LL; j++) {
    // director
    nx[j] += Gamma*hx[j]*dt;
    ny[j] += Gamma*hy[j]*dt;
    nz[j] += Gamma*hz[j]*dt;
    //normalise
    sqrtndotn = sqrt(nx[j]*nx[j] + ny[j]*ny[j] + nz[j]*nz[j]);
    nx[j] /= sqrtndotn;
    ny[j] /= sqrtndotn;
    nz[j] /= sqrtndotn;
    // polarisation
    px[j] += Gamma*hpx[j]*dt;
    py[j] += Gamma*hpy[j]*dt;
    pz[j] += Gamma*hpz[j]*dt;
  }
    
} // end update

/**********************************************************************/
void writeVTKfiles(void)
{
  int j;

  sprintf(vtk_director,"%svtk_director_%d.vtk",prefix,n);
  output.open(vtk_director);
  output.precision(12);

  // header data for the VTK file
  output << "# vtk DataFile Version 3.0" << endl;
  output << "Director Field" << endl;
  output << "ASCII" << endl;
  output << "DATASET STRUCTURED_POINTS" << endl;
  output << "DIMENSIONS " << Lx << " " << Ly << " " << Lz << endl;
  output << "ASPECT_RATIO 1 1 1" << endl;
  output << "ORIGIN 0 0 0" << endl;
  output << "POINT_DATA " << LL << endl;
  output << "VECTORS n double" << endl; // output the director field 
  // currently using NORMALS could also use VECTORS
  //  output << "LOOKUP_TABLE default" << endl;

  sprintf(vtk_polarisation,"%svtk_polarisation_%d.vtk",prefix,n);
  output2.open(vtk_polarisation);
  output2.precision(12);

  // header data for the VTK file
  output2 << "# vtk DataFile Version 3.0" << endl;
  output2 << "Polarisation Field" << endl;
  output2 << "ASCII" << endl;
  output2 << "DATASET STRUCTURED_POINTS" << endl;
  output2 << "DIMENSIONS " << Lx << " " << Ly << " " << Lz << endl;
  output2 << "ASPECT_RATIO 1 1 1" << endl;
  output2 << "ORIGIN 0 0 0" << endl;
  output2 << "POINT_DATA " << LL << endl;
  output2 << "VECTORS p double" << endl; // output the polarisation field 
  // currently using NORMALS could also use VECTORS
  //  output2 << "LOOKUP_TABLE default" << endl;

  for (j=0; j<LL; j++) {
    output << nx[j] << " " << ny[j] << " " << nz[j] << endl;
    output2 << px[j] << " " << py[j] << " " << pz[j] << endl;
  }

  output.close();
  output2.close();
} // end writeVTKfiles

/**********************************************************************/
void writeTwist(void)
{
  int j,k,l,m;
  int xup,xdwn,yup,ydwn,zup,zdwn;
  double Dxnx,Dxny,Dxnz,Dynx,Dyny,Dynz,Dznx,Dzny,Dznz;
  double twist;

  sprintf(vtk_TWIST,"%svtk_TWIST_%d.vtk",prefix,n);
  output.open(vtk_TWIST);
  output.precision(12);

  // header data for the VTK file
  output << "# vtk DataFile Version 3.0" << endl;
  output << "Twist Field" << endl;
  output << "ASCII" << endl;
  output << "DATASET STRUCTURED_POINTS" << endl;
  output << "DIMENSIONS " << Lx << " " << Ly << " " << Lz << endl;
  output << "ASPECT_RATIO 1 1 1" << endl;
  output << "ORIGIN 0 0 0" << endl;
  output << "POINT_DATA " << LL << endl;
  output << "SCALARS twist double 1" << endl; 
  output << "LOOKUP_TABLE default" << endl;

  k=l=m=0;

  for (j=0;j<LL;j++) {
    // define neighbouring nodes 
    xup=j+1; xdwn=j-1; yup=j+Lx; ydwn=j-Lx; zup=j+Lx*Ly; zdwn=j-Lx*Ly;
  
    // correct for periodic boundaries
    if (k==0) {xdwn+=Lx;}
    if (k==Lx-1) {xup-=Lx;}
    if (l==0) {ydwn+=Lx*Ly;}
    if (l==Ly-1) {yup-=Lx*Ly;}
    if (m==0) {zdwn+=LL;}
    if (m==Lz-1) {zup-=LL;}
  
    // cheap fix
    if (Lz==1) {zup=j; zdwn=j;}

#if BC // use one-sided derivatives at boundaries
    if (m==0) {zdwn=j+2*Lx*Ly;}
    if (m==Lz-1) {zup=j-2*Lx*Ly;}
#endif
  
    // calculate first order derivatives
    //    Dxnx = (nx[xup]-nx[xdwn])/2.0;
    Dynx = (nx[yup]-nx[ydwn])/2.0;
    Dznx = (nx[zup]-nx[zdwn])/2.0;

    Dxny = (ny[xup]-ny[xdwn])/2.0;
    //    Dyny = (ny[yup]-ny[ydwn])/2.0;
    Dzny = (ny[zup]-ny[zdwn])/2.0;

    Dxnz = (nz[xup]-nz[xdwn])/2.0;
    Dynz = (nz[yup]-nz[ydwn])/2.0;
    //    Dznz = (nz[zup]-nz[zdwn])/2.0;

#if BC // one-sided derivates at boundaries
    if (m==0) {
      Dznx = (-3.0*nx[j]+4.0*nx[zup]-nx[zdwn])/2.0; 
      Dzny = (-3.0*ny[j]+4.0*ny[zup]-ny[zdwn])/2.0; 
      //      Dznz = (-3.0*nz[j]+4.0*nz[zup]-nz[zdwn])/2.0;
    }
    if (m==Lz-1) {
      Dznx = (3.0*nx[j]-4.0*nx[zdwn]+nx[zup])/2.0; 
      Dzny = (3.0*ny[j]-4.0*ny[zdwn]+ny[zup])/2.0; 
      //      Dznz = (3.0*nz[j]-4.0*nz[zdwn]+nz[zup])/2.0;
    }
#endif

    // calculate twist
    twist = nx[j]*(Dynz-Dzny) + ny[j]*(Dznx-Dxnz) + nz[j]*(Dxny-Dynx);

    output << twist << endl;

    // keep track of boundaries
    k++;
    if (k==Lx) {l++; k=0;}
    if (l==Ly) {m++; l=0;}
  }
  output.close();
} // end writeTwist

/**********************************************************************/
void cholestericPitch(void)
{
  int j,k,l,m;
  int xup,xdwn,yup,ydwn,zup,zdwn;
  double Dxnx,Dxny,Dxnz,Dynx,Dyny,Dynz,Dznx,Dzny,Dznz;
  double Pixx,Pixy,Pixz,Piyy,Piyz,Pizz;
  double trPi,trPisq,normPisq;

  int nrots,emax,enxt;
  double mm[3][3],d[3],v[3][3];

  sprintf(pitch_defect,"%spitch_defect_%d.vtk",prefix,n);
  output.open(pitch_defect);
  output.precision(12);
  // header data for the VTK file
  output << "# vtk DataFile Version 3.0" << endl;
  output << "Pitch Defects" << endl;
  output << "ASCII" << endl;
  output << "DATASET STRUCTURED_POINTS" << endl;
  output << "DIMENSIONS " << Lx << " " << Ly << " " << Lz << endl;
  output << "ASPECT_RATIO 1 1 1" << endl;
  output << "ORIGIN 0 0 0" << endl;
  output << "POINT_DATA " << LL << endl;
  output << "SCALARS DP double 1" << endl; // output the value of -log(normPisq)
  output << "LOOKUP_TABLE default" << endl;
  
  sprintf(pitch_p,"%spitch_p_%d.vtk",prefix,n);
  output2.open(pitch_p);
  output2.precision(12);
  // header data for the VTK file
  output2 << "# vtk DataFile Version 3.0" << endl;
  output2 << "Pitch Axis" << endl;
  output2 << "ASCII" << endl;
  output2 << "DATASET STRUCTURED_POINTS" << endl;
  output2 << "DIMENSIONS " << Lx << " " << Ly << " " << Lz << endl;
  output2 << "ASPECT_RATIO 1 1 1" << endl;
  output2 << "ORIGIN 0 0 0" << endl;
  output2 << "POINT_DATA " << LL << endl;
  output2 << "VECTORS p double" << endl; // output the pitch axis
  // currently using NORMALS could also use VECTORS
  //  output2 << "LOOKUP_TABLE default" << endl;
 
  k=l=m=0;

  for (j=0; j<LL; j++) {
    // define neighbouring nodes 
    xup=j+1; xdwn=j-1; yup=j+Lx; ydwn=j-Lx; zup=j+Lx*Ly; zdwn=j-Lx*Ly;
  
    // correct for periodic boundaries
    if (k==0) {xdwn+=Lx;}
    if (k==Lx-1) {xup-=Lx;}
    if (l==0) {ydwn+=Lx*Ly;}
    if (l==Ly-1) {yup-=Lx*Ly;}
    if (m==0) {zdwn+=LL;}
    if (m==Lz-1) {zup-=LL;}
    
    // cheap fix
    if (Lz==1) {zup=j; zdwn=j;}

#if BC // use one-sided derivatives at boundaries
    if (m==0) {zdwn=j+2*Lx*Ly;}
    if (m==Lz-1) {zup=j-2*Lx*Ly;}
#endif
    
    // calculate first order derivatives
    Dxnx = (nx[xup]-nx[xdwn])/2.0;
    Dynx = (nx[yup]-nx[ydwn])/2.0;
    Dznx = (nx[zup]-nx[zdwn])/2.0;

    Dxny = (ny[xup]-ny[xdwn])/2.0;
    Dyny = (ny[yup]-ny[ydwn])/2.0;
    Dzny = (ny[zup]-ny[zdwn])/2.0;

    Dxnz = (nz[xup]-nz[xdwn])/2.0;
    Dynz = (nz[yup]-nz[ydwn])/2.0;
    Dznz = (nz[zup]-nz[zdwn])/2.0;
    
#if BC // one-sided derivates at boundaries
    if (m==0) {
      Dznx = (-3.0*nx[j]+4.0*nx[zup]-nx[zdwn])/2.0; 
      Dzny = (-3.0*ny[j]+4.0*ny[zup]-ny[zdwn])/2.0; 
      Dznz = (-3.0*nz[j]+4.0*nz[zup]-nz[zdwn])/2.0;
    }
    if (m==Lz-1) {
      Dznx = (3.0*nx[j]-4.0*nx[zdwn]+nx[zup])/2.0; 
      Dzny = (3.0*ny[j]-4.0*ny[zdwn]+ny[zup])/2.0; 
      Dznz = (3.0*nz[j]-4.0*nz[zdwn]+nz[zup])/2.0;
    }
#endif
    
    // Pi is -- but i do not include the factor of 1/4
    // Pi_{AB} = (1/4) e_{Alk} ( nl[j]*DknB + nl[j]*DBnk - nB[j]*nl[j]*(nx[j]*Dxnk + ny[j]*Dynk + nz[j]*Dznk) ) + (1/4) e_{Blk} ( nl[j]*DknA + nl[j]*DAnk - nA[j]*nl[j]*(nx[j]*Dxnk + ny[j]*Dynk + nz[j]*Dznk) )
    Pixx = 2.0*(ny[j]*Dznx + ny[j]*Dxnz - nx[j]*ny[j]*(nx[j]*Dxnz + ny[j]*Dynz + nz[j]*Dznz)) - 2.0*(nz[j]*Dynx + nz[j]*Dxny - nx[j]*nz[j]*(nx[j]*Dxny + ny[j]*Dyny + nz[j]*Dzny)); 
    Pixy = (ny[j]*Dzny + ny[j]*Dynz - ny[j]*ny[j]*(nx[j]*Dxnz + ny[j]*Dynz + nz[j]*Dznz)) - (nz[j]*Dyny + nz[j]*Dyny - ny[j]*nz[j]*(nx[j]*Dxny + ny[j]*Dyny + nz[j]*Dzny)) + (nz[j]*Dxnx + nz[j]*Dxnx - nx[j]*nz[j]*(nx[j]*Dxnx + ny[j]*Dynx + nz[j]*Dznx)) - (nx[j]*Dznx + nx[j]*Dxnz - nx[j]*nx[j]*(nx[j]*Dxnz + ny[j]*Dynz + nz[j]*Dznz));
    Pixz = (ny[j]*Dznz + ny[j]*Dznz - nz[j]*ny[j]*(nx[j]*Dxnz + ny[j]*Dynz + nz[j]*Dznz)) - (nz[j]*Dynz + nz[j]*Dzny - nz[j]*nz[j]*(nx[j]*Dxny + ny[j]*Dyny + nz[j]*Dzny)) + (nx[j]*Dynx + nx[j]*Dxny - nx[j]*nx[j]*(nx[j]*Dxny + ny[j]*Dyny + nz[j]*Dzny)) - (ny[j]*Dxnx + ny[j]*Dxnx - nx[j]*ny[j]*(nx[j]*Dxnx + ny[j]*Dynx + nz[j]*Dznx));
    Piyy = 2.0*(nz[j]*Dxny + nz[j]*Dynx - ny[j]*nz[j]*(nx[j]*Dxnx + ny[j]*Dynx + nz[j]*Dznx)) - 2.0*(nx[j]*Dzny + nx[j]*Dynz - ny[j]*nx[j]*(nx[j]*Dxnz + ny[j]*Dynz + nz[j]*Dznz));
    Piyz = (nz[j]*Dxnz + nz[j]*Dznx - nz[j]*nz[j]*(nx[j]*Dxnx + ny[j]*Dynx + nz[j]*Dznx)) - (nx[j]*Dznz + nx[j]*Dznz - nz[j]*nx[j]*(nx[j]*Dxnz + ny[j]*Dynz + nz[j]*Dznz)) + (nx[j]*Dyny + nx[j]*Dyny - ny[j]*nx[j]*(nx[j]*Dxny + ny[j]*Dyny + nz[j]*Dzny)) - (ny[j]*Dxny + ny[j]*Dynx - ny[j]*ny[j]*(nx[j]*Dxnx + ny[j]*Dynx + nz[j]*Dznx));
    Pizz = 2.0*(nx[j]*Dynz + nx[j]*Dzny - nz[j]*nx[j]*(nx[j]*Dxny + ny[j]*Dyny + nz[j]*Dzny)) - 2.0*(ny[j]*Dxnz + ny[j]*Dznx - nz[j]*ny[j]*(nx[j]*Dxnx + ny[j]*Dynx + nz[j]*Dznx));

    // the 2-norm of Pi vanishes on umbilics
    normPisq = Pixx*Pixx + 2.0*Pixy*Pixy + 2.0*Pixz*Pixz + Piyy*Piyy + 2.0*Piyz*Piyz + Pizz*Pizz;

    output << -log(normPisq) << endl;

    // find eigenvalues and eigenvectors of Pi
    mm[0][0]=Pixx;
    mm[0][1]=Pixy;
    mm[0][2]=Pixz;
    mm[1][0]=Pixy;
    mm[1][1]=Piyy;
    mm[1][2]=Piyz;
    mm[2][0]=Pixz;
    mm[2][1]=Piyz;
    mm[2][2]=Pizz;
    jacobi(mm,d,v,&nrots);

    if (d[0] > d[1]) {
      emax=0;
      enxt=1;
    }
    else {
      emax=1;
      enxt=0;
    }
    if (d[2] > d[emax]) {
      emax=2;
    }
    else if (d[2] > d[enxt]) {
      enxt=2;
    }

    output2 << v[0][emax] << " " << v[1][emax] << " " << v[2][emax] << endl;
    //    output2 << v[0][enxt] << " " << v[1][enxt] << " " << v[2][enxt] << endl;

    // keep track of boundaries
    k++;
    if (k==Lx) {l++; k=0;}
    if (l==Ly) {m++; l=0;}
  }
  output.close();
  output2.close();
} // end cholestericPitch

/**********************************************************************/
void writeBENDfiles(void)
{
  int j,k,l,m;
  int xup,xdwn,yup,ydwn,zup,zdwn;
  double Dxnx,Dxny,Dxnz,Dynx,Dyny,Dynz,Dznx,Dzny,Dznz;
  double bx,by,bz;

  sprintf(vtk_BEND,"%svtk_BEND_%d.vtk",prefix,n);
  output.open(vtk_BEND);
  output.precision(12);

  // header data for the VTK file
  output << "# vtk DataFile Version 3.0" << endl;
  output << "Bend Field" << endl;
  output << "ASCII" << endl;
  output << "DATASET STRUCTURED_POINTS" << endl;
  output << "DIMENSIONS " << Lx << " " << Ly << " " << Lz << endl;
  output << "ASPECT_RATIO 1 1 1" << endl;
  output << "ORIGIN 0 0 0" << endl;
  output << "POINT_DATA " << LL << endl;
  output << "VECTORS b double" << endl; // output the bend vector field 
  // currently using NORMALS could also use VECTORS
  //  output << "LOOKUP_TABLE default" << endl;

  k=l=m=0;

  for (j=0;j<LL;j++) {
    // define neighbouring nodes 
    xup=j+1; xdwn=j-1; yup=j+Lx; ydwn=j-Lx; zup=j+Lx*Ly; zdwn=j-Lx*Ly;
  
    // correct for periodic boundaries
    if (k==0) {xdwn+=Lx;}
    if (k==Lx-1) {xup-=Lx;}
    if (l==0) {ydwn+=Lx*Ly;}
    if (l==Ly-1) {yup-=Lx*Ly;}
    if (m==0) {zdwn+=LL;}
    if (m==Lz-1) {zup-=LL;}
  
    // cheap fix
    if (Lz==1) {zup=j; zdwn=j;}

#if BC // use one-sided derivatives at boundaries
    if (m==0) {zdwn=j+2*Lx*Ly;}
    if (m==Lz-1) {zup=j-2*Lx*Ly;}
#endif
  
    // calculate first order derivatives
    Dxnx = (nx[xup]-nx[xdwn])/2.0;
    Dynx = (nx[yup]-nx[ydwn])/2.0;
    Dznx = (nx[zup]-nx[zdwn])/2.0;

    Dxny = (ny[xup]-ny[xdwn])/2.0;
    Dyny = (ny[yup]-ny[ydwn])/2.0;
    Dzny = (ny[zup]-ny[zdwn])/2.0;

    Dxnz = (nz[xup]-nz[xdwn])/2.0;
    Dynz = (nz[yup]-nz[ydwn])/2.0;
    Dznz = (nz[zup]-nz[zdwn])/2.0;

#if BC // one-sided derivates at boundaries
    if (m==0) {
      Dznx = (-3.0*nx[j]+4.0*nx[zup]-nx[zdwn])/2.0; 
      Dzny = (-3.0*ny[j]+4.0*ny[zup]-ny[zdwn])/2.0; 
      Dznz = (-3.0*nz[j]+4.0*nz[zup]-nz[zdwn])/2.0;
    }
    if (m==Lz-1) {
      Dznx = (3.0*nx[j]-4.0*nx[zdwn]+nx[zup])/2.0; 
      Dzny = (3.0*ny[j]-4.0*ny[zdwn]+ny[zup])/2.0; 
      Dznz = (3.0*nz[j]-4.0*nz[zdwn]+nz[zup])/2.0;
    }
#endif

    // calculate bend
    bx = nx[j]*Dxnx + ny[j]*Dynx + nz[j]*Dznx;
    by = nx[j]*Dxny + ny[j]*Dyny + nz[j]*Dzny;
    bz = nx[j]*Dxnz + ny[j]*Dynz + nz[j]*Dznz;

    output << bx << " " << by << " " << bz << endl;

    // keep track of boundaries
    k++;
    if (k==Lx) {l++; k=0;}
    if (l==Ly) {m++; l=0;}
  }
  output.close();
} // end writeBENDfiles


/**************************************************************************/
/*  The following routines are based on those given in Numerical Recipes  */
/*   and are copied here from Alex and Davide's lattice Boltzmann code    */
/**************************************************************************/

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);
#define n 3
void jacobi(double (*a)[n], double d[], double (*v)[n], int *nrot)
{
  int j,iq,ip,i;
  double tresh,theta,tau,t,sm,s,h,g,c;
  double b[n],z[n];

  for (ip=0;ip<n;ip++) {
    for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }
  for (ip=0;ip<n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=0;ip< n-1;ip++) {
      for (iq=ip+1;iq<n;iq++)
	sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0) {
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=0;ip<n-1;ip++) {
      for (iq=ip+1;iq<n;iq++) {
	g=100.0*fabs(a[ip][iq]);
	if (i > 4 && (fabs(d[ip])+g) == fabs(d[ip])
	    && (fabs(d[iq])+g) == fabs(d[iq]))
	  a[ip][iq]=0.0;
	else if (fabs(a[ip][iq]) > tresh) {
	  h=d[iq]-d[ip];
	  if ((fabs(h)+g) == fabs(h))
	    t=(a[ip][iq])/h;
	  else {
	    theta=0.5*h/(a[ip][iq]);
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq]=0.0;
	  for (j=0;j<=ip-1;j++) {
	    ROTATE(a,j,ip,j,iq)
	      }
	  for (j=ip+1;j<=iq-1;j++) {
	    ROTATE(a,ip,j,j,iq)
	      }
	  for (j=iq+1;j<n;j++) {
	    ROTATE(a,ip,j,iq,j)
	      }
	  for (j=0;j<n;j++) {
	    ROTATE(v,j,ip,j,iq)
	      }
	  ++(*nrot);
	}
      }
    }
    for (ip=0;ip<n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  cout << "Too many iterations in routine jacobi" << endl;
  exit(0);
}
#undef n
#undef ROTATE


