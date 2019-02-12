/**************************************************************/
/*      Finite difference code for Twist Bend Nematics        */
/*         created November 2018, Gareth Alexander            */
/**************************************************************/

#include "twistbend.h"
#include "ReadingWriting.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <complex>
#include <string>
#include <complex>
#include <omp.h>

using namespace std;

int main(int argc, char** argv) 
{

    double* nx=new double[LL];
    double* ny=new double[LL];
    double* nz=new double[LL];

    double* px=new double[LL];
    double* py=new double[LL];
    double* pz=new double[LL];

    double* pmag=new double[LL];

    double* bx=new double[LL];
    double* by=new double[LL];
    double* bz=new double[LL];

    double* bmag=new double[LL];

    double* tx=new double[LL];
    double* ty=new double[LL];
    double* tz=new double[LL];

    double* hx=new double[LL];
    double* hy=new double[LL];
    double* hz=new double[LL];

    double* hpx=new double[LL];
    double* hpy=new double[LL];
    double* hpz=new double[LL];


    // the structure we store the bend zeros in
    vector<knotcurve> knotcurves;
    // GSL initialization
    const gsl_multimin_fminimizer_type *Type;
    Type = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *minimizerstate;
    minimizerstate = gsl_multimin_fminimizer_alloc (Type,2);

    int n =0;
    int step=0;
    startconfig(n, nx,ny,nz,px,py,pz);
    cout << "starting simulation" << endl;

#pragma omp parallel default(none) shared(nx,ny,nz,px,py,pz,pmag, hx,hy,hz,hpx,hpy,hpz, bx,by,bz,bmag, tx,ty,tz, knotcurves,minimizerstate, Type, n,step,cout)
    {
        while(n<=Nmax)
        {
#pragma omp single
            {
                if (n%1000==0)
                {
                    cout << "timestep " << n << endl;
                }

                if (step==stepskip)
                {
                    cout << "writing VTK files at timestep " << n << endl;
                    computeBendAndCurlofCirculation(n,nx,ny,nz,bx,by,bz,bmag,px,py,pz,pmag, tx,ty,tz);
                    FindBendZeros(pmag,tx, ty,tz,knotcurves,n,minimizerstate);
                    writeVTKfiles(n,nx,ny,nz,px,py,pz,bx,by,bz,tx,ty,tz);       // output VTK files for use with ParaView
                    print_knot(n,knotcurves);

                    step=0;
                }
                step++;
                n++;
            }
            update(nx,ny,nz,px, py,pz, hx, hy,hz,hpx, hpy,hpz);
        }
    }
} // end main

/**********************************************************************/
void startconfig(int & n, double* nx, double* ny,double* nz,double* px, double* py,double* pz)
{
    switch(InitialisationMethod)
    {
        case FROM_FUNCTION:
            {
                // k,l,m  are x,y,z indices, j is 1-d loop index.
                int j,k,l,m;
                int HELICONICAL,HOPF,SQUARE,HOPF2;
                double k0,l0,m0;

                // define texture to simulate -- 0 or 1
                HELICONICAL = 0; // standard heliconical texture
                HOPF = 1; // Hopf texture
                HOPF2 = 0; // Hopf texture
                SQUARE = 0;

                // initial configuration
                k=l=m=0;
                k0 = Lx/2.0 - 0.5; l0 = Ly/2.0 - 0.5; m0 = Lz/2.0 - 0.5; // offset from the lattice points

                // standard heliconical texture
                if (HELICONICAL == 1)
                {
                    for (j=0; j<LL; j++)
                    {
                        // director field
                        nx[j] = sin(thetah)*cos(qh*m);
                        ny[j] = sin(thetah)*sin(qh*m);
                        nz[j] = cos(thetah);
                        // polarisation
                        px[j] = -sin(qh*m);
                        py[j] = cos(qh*m);
                        pz[j] = 0.0;
                        // perversion -- try this
                        //      if ((m>0.25*Lz)&&(m<0.75*Lz)) {
                        //if (m>0.5*Lz) {
                        //	ny[j] *= -1.0;
                        //py[j] *= -1.0;
                        // test
                        //	nx[j] *= -1.0;
                        //	nz[j] *= -1.0;
                        // }}
                        // deal with the periodic boundaries
                        k++;
                        if (k==Lx) {l++; k=0;}
                        if (l==Ly) {m++; l=0;}
                    }
                }
                // hopf texture
                // THIS ISNT CORRECT YET
                if (HOPF == 1)
                {
                    int xup,xdwn,yup,ydwn,zup,zdwn;
                    double rho,theta,phi, rr, norm;    // variables for hopf texture
                    double RR = 0.45*Lz;  // scale
                    double RR2 =0.25*Lz;  // scale
                    for (j=0; j<LL; j++)
                    {
                        // heliconical
                        // director field
                        nx[j] = sin(thetah)*cos(qh*m);
                        ny[j] = sin(thetah)*sin(qh*m);
                        nz[j] = cos(thetah);
                        // polarisation
                        px[j] = -sin(qh*m);
                        py[j] = cos(qh*m);
                        pz[j] = 0.0;
                        // back to Hopf
                        rr = sqrt((k-k0)*(k-k0)+(l-l0)*(l-l0));
                        rho = sqrt((rr-RR)*(rr-RR)+(m-m0)*(m-m0));
                        if (rho < RR2)
                        {
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
                            nx[j] = -sin(M_PI*rho/RR2)*sin(phi+theta); // could also use -theta; they seem to be equivalent in energy
                            ny[j] = sin(M_PI*rho/RR2)*cos(phi+theta); // could also use -theta
                            nz[j] = -cos(M_PI*rho/RR2);

                            px[j] = 0;
                            py[j] = 0;
                            pz[j] = 0.0;
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
                    for (j=0; j<LL; j++)
                    {
                        // hack here ??
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
                        if (rho < RR2) {
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

                if (HOPF2 == 1)
                {
                    for (j=0; j<LL; j++)
                    {
                        complex<double> z1 =(double)(2)*(complex<double>(k,0)+complex<double>(0,l))/(double)(k*k+l*l+m*m+1);
                        complex<double> z2= (  (double)(2.0*m)+complex<double>(0,k*k+l*l+m*m-1) )/(double)(k*k+l*l+m*m+1);
                        complex<double> nxplusiny = (double)(2)*z1*conj(z2);
                        nx[j] = real(nxplusiny);
                        ny[j] = imag(nxplusiny);
                        nz[j] = abs(z1)*abs(z1)-abs(z2)*abs(z2);
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
                    int xup,xdwn,yup,ydwn,zup,zdwn,norm;
                    for (j=0; j<LL; j++)
                    {
                        // define neighbouring nodes
                        xup=j+1; xdwn=j-1; yup=j+Lx; ydwn=j-Lx; zup=j+Lx*Ly; zdwn=j-Lx*Ly;
                        // correct for periodic boundaries
                        if (k==0) {xdwn+=Lx;}
                        if (k==Lx-1) {xup-=Lx;}
                        if (l==0) {ydwn+=Lx*Ly;}
                        if (l==Ly-1) {yup-=Lx*Ly;}
                        if (m==0) {zdwn+=LL;}
                        if (m==Lz-1) {zup-=LL;}

                        px[j] = nx[j]*(nx[xup]-nx[xdwn])+ny[j]*(nx[yup]-nx[ydwn])+nz[j]*(nx[zup]-nx[zdwn]);
                        py[j] = nx[j]*(ny[xup]-ny[xdwn])+ny[j]*(ny[yup]-ny[ydwn])+nz[j]*(ny[zup]-ny[zdwn]);
                        pz[j] = nx[j]*(nz[xup]-nz[xdwn])+ny[j]*(nz[yup]-nz[ydwn])+nz[j]*(nz[zup]-nz[zdwn]);
                        // normalise
                        norm = sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]);
                        px[j] /= norm; py[j] /= norm; pz[j] /= norm;

#if BC // normal anchoring along z, so make the polarization zero here
                        if (m==0) {px[j] = 0.0; py[j] = 0.0; pz[j] = 0.0;}
                        if (m==Lz-1) {px[j] = 0.0; py[j] = 0.0; pz[j] = 0.0;}
#endif
                        // deal with the periodic boundaries
                        k++;
                        if (k==Lx) {l++; k=0;}
                        if (l==Ly) {m++; l=0;}
                    }
                } // end hopf

                // square
                if (SQUARE == 1)
                {
                    double r,R,phi;
                    R = 0.8*Lx/2.0;

                    for (j=0; j<LL; j++) {
                        // default is heliconical texture
                        // director field
                        nx[j] = -sin(thetah)*cos(qh*m);
                        ny[j] = -sin(thetah)*sin(qh*m);
                        nz[j] = -cos(thetah);
                        // polarisation
                        px[j] = -sin(qh*m);
                        py[j] = cos(qh*m);
                        pz[j] = 0.0;
                        // replace a cylinder with a twisted Skyrmion -- testing!
                        r = sqrt((k-k0)*(k-k0)+(l-l0)*(l-l0));
                        if (r<R) {
                            phi = atan2(l-l0,k-k0);

                            // nz is common to the twisted and untwisted cases
                            nz[j] = cos(M_PI*r/R);

                            // degree +1 , untwisted
                            nx[j] = sin(M_PI*r/R)*sin(phi);
                            ny[j] = -sin(M_PI*r/R)*cos(phi);
                            // polarisation for +1 untwisted
                            px[j] = -sin(M_PI*r/R)*sin(M_PI*r/R)*cos(phi)/r;
                            py[j] = -sin(M_PI*r/R)*sin(M_PI*r/R)*sin(phi)/r;
                            pz[j] = 0.0;

                            // twisted
                            //nx[j] = sin(M_PI*r/R)*sin(phi+2.0*M_PI*m/Lz);
                            //ny[j] = -sin(M_PI*r/R)*cos(phi+2.0*M_PI*m/Lz);
                            // polarisation for twisted
                            //px[j] = cos(M_PI*r/R)*sin(2.0*M_PI*m/Lz)*sin(phi+2.0*M_PI*m/Lz)+((2.0*R/Lz)*cos(M_PI*r/R)-(R/(M_PI*r))*sin(M_PI*r/R)*cos(2.0*M_PI*m/Lz))*cos(phi+2.0*M_PI*m/Lz);
                            //py[j] = -cos(M_PI*r/R)*sin(2.0*M_PI*m/Lz)*cos(phi+2.0*M_PI*m/Lz)+((2.0*R/Lz)*cos(M_PI*r/R)-(R/(M_PI*r))*sin(M_PI*r/R)*cos(2.0*M_PI*m/Lz))*sin(phi+2.0*M_PI*m/Lz);
                            //pz[j] = -sin(M_PI*r/R)*sin(2.0*M_PI*m/Lz);

                            // degree -1
                            //	nx[j] = -sin(M_PI*r/R)*sin(phi);
                            //	ny[j] = -sin(M_PI*r/R)*cos(phi);
                        }
                        // deal with the periodic boundaries
                        k++;
                        if (k==Lx) {l++; k=0;}
                        if (l==Ly) {m++; l=0;}
                    }
                } // end square
                break;
            }
        case FROM_FILE:
            {
                cout << "Reading input file...\n";
                if(file_read(nx,ny,nz,px,py,pz)){ cout << "somethings not right with the read...\n";}
                // get the start time -  we hack this together as so:
                n = starttime;
                break;
            }
    }
}

/**********************************************************************/
void update(double* nx, double* ny,double* nz,double* px, double* py,double* pz, double* hx, double* hy,double* hz,double* hpx, double* hpy,double* hpz)
{
    int k,l,m;
    int xup,xdwn,yup,ydwn,zup,zdwn;
    double Dxnx,Dynx,Dznx,Dxxnx,Dyynx,Dzznx;
    double Dxny,Dyny,Dzny,Dxxny,Dyyny,Dzzny;
    double Dxnz,Dynz,Dznz,Dxxnz,Dyynz,Dzznz;
    double hdotn,sqrtndotn;
    double Dxpx,Dypx,Dzpx,Dxxpx,Dyypx,Dzzpx;
    double Dxpy,Dypy,Dzpy,Dxxpy,Dyypy,Dzzpy;
    double Dxpz,Dypz,Dzpz,Dxxpz,Dyypz,Dzzpz;
    /* Calculate derivatives, molecular field and the energy */

#pragma omp for
    for ( m=0; m<Lz; m++) {
        for (l=0; l<Ly; l++) {
            for (k=0; k<Lx; k++) {
                // where am I in the 1-D array?
                int j = pt(k,l,m);
                // define neighbouring nodes in the 1-D array
                xup=pt(k+1,l,m);
                xdwn=pt(k-1,l,m);
                yup=pt(k,l+1,m);
                ydwn=pt(k,l-1,m);
                zup=pt(k,l,m+1);
                zdwn=pt(k,l,m-1);
                // correct for periodic boundaries
                if (k==0) {xdwn=pt(Lx-1,l,m);}
                if (k==Lx-1) {xup=pt(0,l,m);}
                if (l==0) {ydwn=pt(k,Ly-1,m);}
                if (l==Ly-1) {yup=pt(k,0,m);}
                if (m==0) {zdwn=pt(k,l,Lz-1);}
                if (m==Lz-1) {zup=pt(k,l,0);}
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
                hx[j] = K*(Dxxnx+Dyynx+Dzznx) + lambda*(px[j]*Dxnx+py[j]*Dxny+pz[j]*Dxnz - px[j]*(Dxnx+Dyny+Dznz) - (nx[j]*Dxpx+ny[j]*Dypx+nz[j]*Dzpx));
                hy[j] = K*(Dxxny+Dyyny+Dzzny) + lambda*(px[j]*Dynx+py[j]*Dyny+pz[j]*Dynz - py[j]*(Dxnx+Dyny+Dznz) - (nx[j]*Dxpy+ny[j]*Dypy+nz[j]*Dzpy));
                hz[j] = K*(Dxxnz+Dyynz+Dzznz) + lambda*(px[j]*Dznx+py[j]*Dzny+pz[j]*Dznz - pz[j]*(Dxnx+Dyny+Dznz) - (nx[j]*Dxpz+ny[j]*Dypz+nz[j]*Dzpz));

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
            }
        }
    }
    // now do the update itself
#pragma omp for
    for (int j=0; j<Lx*Ly*Lz; j++)
    {
        // director
        nx[j] += hx[j]*dt;
        ny[j] += hy[j]*dt;
        nz[j] += hz[j]*dt;
        //normalise
        sqrtndotn = sqrt(nx[j]*nx[j] + ny[j]*ny[j] + nz[j]*nz[j]);
        nx[j] /= sqrtndotn;
        ny[j] /= sqrtndotn;
        nz[j] /= sqrtndotn;
        // polarisation
        px[j] += hpx[j]*dt;
        py[j] += hpy[j]*dt;
        pz[j] += hpz[j]*dt;

    }

} // end update

void computeBendAndCurlofCirculation(const int n, const double* nx,const double* ny,const double* nz, double* bx, double* by, double* bz, double* bmag, const double* px, const double* py, const double* pz, double* pmag, double* tx, double* ty, double* tz)
{

    // first up, compute magp
    int j;
    for (j=0;j<LL;j++) pmag[j] = sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]);

    // arrays storing the circulaiton, which we will take the curl of
    double* cx=new double[LL];
    double* cy=new double[LL];
    double* cz=new double[LL];

    int k,l,m;
    int xup,xdwn,yup,ydwn,zup,zdwn;
    double Dxnx,Dxny,Dxnz,Dynx,Dyny,Dynz,Dznx,Dzny,Dznz;
    double Dxbx,Dxby,Dxbz,Dybx,Dyby,Dybz,Dzbx,Dzby,Dzbz;
    double Dxcx,Dxcy,Dxcz,Dycx,Dycy,Dycz,Dzcx,Dzcy,Dzcz;

    k=l=m=0;
    for (j=0;j<LL;j++) 
    {
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
        bx[j] = nx[j]*Dxnx + ny[j]*Dynx + nz[j]*Dznx;
        by[j] = nx[j]*Dxny + ny[j]*Dyny + nz[j]*Dzny;
        bz[j] = nx[j]*Dxnz + ny[j]*Dynz + nz[j]*Dznz;
        bmag[j] = sqrt(bx[j]*bx[j]+by[j]*by[j]+bz[j]*bz[j]);

        // keep track of boundaries
        k++;
        if (k==Lx) {l++; k=0;}
        if (l==Ly) {m++; l=0;}
    }

    k=l=m=0;
    // now compute the ciruclation, B^i = Jb . d_i b
    for (j=0;j<LL;j++) 
    {
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
        Dxbx = (bx[xup]-bx[xdwn])/2.0;
        Dybx = (bx[yup]-bx[ydwn])/2.0;
        Dzbx = (bx[zup]-bx[zdwn])/2.0;

        Dxby = (by[xup]-by[xdwn])/2.0;
        Dyby = (by[yup]-by[ydwn])/2.0;
        Dzby = (by[zup]-by[zdwn])/2.0;

        Dxbz = (bz[xup]-bz[xdwn])/2.0;
        Dybz = (bz[yup]-bz[ydwn])/2.0;
        Dzbz = (bz[zup]-bz[zdwn])/2.0;

#if BC // one-sided derivates at boundaries
        if (m==0) {
            Dzbx = (-3.0*bx[j]+4.0*bx[zup]-bx[zdwn])/2.0;
            Dzby = (-3.0*by[j]+4.0*by[zup]-by[zdwn])/2.0;
            Dzbz = (-3.0*bz[j]+4.0*bz[zup]-bz[zdwn])/2.0;
        }
        if (m==Lz-1) {
            Dzbx = (3.0*bx[j]-4.0*bx[zdwn]+bx[zup])/2.0;
            Dzby = (3.0*by[j]-4.0*by[zdwn]+by[zup])/2.0;
            Dzbz = (3.0*bz[j]-4.0*bz[zdwn]+bz[zup])/2.0;
        }
#endif

        // get b_perp at this point via a cross product:
        double bperpx = ny[j]*bz[j] - nz[j]*by[j];
        double bperpy = nz[j]*bx[j] - nx[j]*bz[j];
        double bperpz = nx[j]*by[j] - ny[j]*bx[j];

        cx[j] = bperpx* Dxbx + bperpy* Dxby+ bperpz* Dxbz;
        cy[j] = bperpx* Dybx + bperpy* Dyby+ bperpz* Dybz;
        cz[j] = bperpx* Dzbx + bperpy* Dzby+ bperpz* Dzbz;

        // keep track of boundaries
        k++;
        if (k==Lx) {l++; k=0;}
        if (l==Ly) {m++; l=0;}
    }

    // and finally compute the curl of the circulation
    k=l=m=0;
    for (j=0;j<LL;j++) 
    {
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
        Dxcx = (cx[xup]-cx[xdwn])/2.0;
        Dycx = (cx[yup]-cx[ydwn])/2.0;
        Dzcx = (cx[zup]-cx[zdwn])/2.0;

        Dxcy = (cy[xup]-cy[xdwn])/2.0;
        Dycy = (cy[yup]-cy[ydwn])/2.0;
        Dzcy = (cy[zup]-cy[zdwn])/2.0;

        Dxcz = (cz[xup]-cz[xdwn])/2.0;
        Dycz = (cz[yup]-cz[ydwn])/2.0;
        Dzcz = (cz[zup]-cz[zdwn])/2.0;

#if BC // one-sided derivates at boundaries
        if (m==0) {
            Dzcx = (-3.0*cx[j]+4.0*cx[zup]-cx[zdwn])/2.0;
            Dzcy = (-3.0*cy[j]+4.0*cy[zup]-cy[zdwn])/2.0;
            Dzcz = (-3.0*cz[j]+4.0*cz[zup]-cz[zdwn])/2.0;
        }
        if (m==Lz-1) {
            Dzcx = (3.0*cx[j]-4.0*cx[zdwn]+cx[zup])/2.0;
            Dzcy = (3.0*cy[j]-4.0*cy[zdwn]+cy[zup])/2.0;
            Dzcz = (3.0*cz[j]-4.0*cz[zdwn]+cz[zup])/2.0;
        }
#endif
        tx[j] = Dycz - Dzcy;
        ty[j] = Dzcx - Dxcz;
        tz[j] = Dxcy - Dycx;

        // keep track of boundaries
        k++;
        if (k==Lx) {l++; k=0;}
        if (l==Ly) {m++; l=0;}
    }
}

void FindBendZeros(double *magb,double* tx,double* ty,double* tz, vector<knotcurve>& knotcurves,double n, gsl_multimin_fminimizer* minimizerstate)
{

    // first thing, clear the knotcurve object before we begin writing a new one
    knotcurves.clear(); //empty vector with knot curve points

    // initialise the tricubic interpolator for ucvmag
    likely::TriCubicInterpolator interpolatedmagb(magb, 1, Lx,Ly,Lz);

    // an array holding points we have already visited, so we dont keep tracing out the same bend zero 
    vector<int> marked(Lx*Ly*Lz,0);


    // settings
    double threshold =1.2;


    int c =0;
    bool knotexists = true;
    while(knotexists)
    {
        double   bmin = 10000000;// I dunno some big number 
        int kmin=-1,lmin=-1,mmin=-1;
        for(int k=0;k<Lx;k++)
        {
            for(int l=0; l<Ly; l++)
            {
                for(int m=0; m<Lz; m++) 
                {
                    int n = pt(k,l,m);
                    if( magb[n] < bmin && marked[n]==0)
                    {
                        bmin = magb[n];
                        kmin = k;
                        lmin = l;
                        mmin=m;

                    }
                }
            }
        }

        if(bmin>threshold)knotexists = false;

        if(knotexists)
        {
            knotcurves.push_back(knotcurve() );
            knotcurves[c].knotcurve.push_back(knotpoint());
            knotcurves[c].knotcurve[0].xcoord=kmin;
            knotcurves[c].knotcurve[0].ycoord=lmin;
            knotcurves[c].knotcurve[0].zcoord=mmin;

            int s=1;
            bool finish=false;
            /*calculate local direction of grad u x grad v (the tangent to the knot curve) at point s-1, then move to point s by moving along tangent + unit confinement force*/
            while (finish==false)
            {   
                //pop a new point on the stack
                knotcurves[c].knotcurve.push_back(knotpoint());

                // the tangent vector cirlcirculation 
                double txs=0;
                double tys=0;
                double tzs=0;

                /**Find nearest gridpoint**/
                int idwn = (int) ((knotcurves[c].knotcurve[s-1].xcoord));
                int jdwn = (int) ((knotcurves[c].knotcurve[s-1].ycoord));
                int kdwn = (int) ((knotcurves[c].knotcurve[s-1].zcoord));
                /*curve to gridpoint down distance*/
                double xd = knotcurves[c].knotcurve[s-1].xcoord - idwn;
                double yd = knotcurves[c].knotcurve[s-1].ycoord - jdwn;
                double zd = knotcurves[c].knotcurve[s-1].zcoord - kdwn;
                for(int m=0;m<8;m++)  //linear interpolation from 8 nearest neighbours
                {
                    /* Work out increments*/
                    int iinc = m%2;
                    int jinc = (m/2)%2;
                    int kinc = (m/4)%2;
                    /*Loop over nearest points*/
                    int i = idwn +iinc;
                    int j = jdwn+jinc;
                    int k = kdwn+kinc;
                    double prefactor = (1-iinc + pow(-1,1+iinc)*xd)*(1-jinc + pow(-1,1+jinc)*yd)*(1-kinc + pow(-1,1+kinc)*zd);
                    /*interpolate curlcirc over nearest points*/
                    txs += prefactor*tx[pt(i,j,k)];
                    tys += prefactor*ty[pt(i,j,k)];
                    tzs += prefactor*tz[pt(i,j,k)];
                }
                // normalise it
                double norm = sqrt(txs*txs + tys*tys + tzs*tzs);
                txs = txs/norm; 
                tys = tys/norm;
                tzs = tzs/norm; 
                // walk along the tangent to the bend zero some distance
                double walkstepsize = 1;
                double testx = knotcurves[c].knotcurve[s-1].xcoord + walkstepsize*txs;
                double testy = knotcurves[c].knotcurve[s-1].ycoord + walkstepsize*tys;
                double testz = knotcurves[c].knotcurve[s-1].zcoord + walkstepsize*tzs;

                // we are flowing along an integral curve of the curlcirculation here; if everything was perfect this would pick out the bend zero. However, discretization error etc. will 
                // cause us to flow off the bend zero over time; we need to correct ourselves. We do this by minimizing magb in a plane transverse to the curlcirculation vector. 
                // We construct a triple spanning this plane by taking grad b (if needed projecting onto the plane perp to curlc) and the cross product.

                // recompute at our test point
                idwn = (int) (testx) ;
                jdwn = (int) (testy) ;
                kdwn = (int) (testz) ;
                txs=0;
                tys=0;
                tzs=0;
                for(int m=0;m<8;m++)  //linear interpolation from 8 nearest neighbours
                {
                    /* Work out increments*/
                    int iinc = m%2;
                    int jinc = (m/2)%2;
                    int kinc = (m/4)%2;
                    /*Loop over nearest points*/
                    int i = idwn +iinc;
                    int j = jdwn+jinc;
                    int k = kdwn+kinc;
                    double prefactor = (1-iinc + pow(-1,1+iinc)*xd)*(1-jinc + pow(-1,1+jinc)*yd)*(1-kinc + pow(-1,1+kinc)*zd);
                    /*interpolate curlcirc over nearest points*/
                    txs += prefactor*tx[pt(i,j,k)];
                    tys += prefactor*ty[pt(i,j,k)];
                    tzs += prefactor*tz[pt(i,j,k)];
                }
                // now we get some frame perpendicular to this tangent vector. Ill get this by using the cartesian frame and projecting. 
                double e1x,e1y,e1z;
                double e2x,e2y,e2z;
                // the cartesian basis
                vector<double> cartesian(3);
                double maxnorm= -1;
                for(int i=0; i<2; i++)
                {
                    std::fill(cartesian.begin(), cartesian.end(), 0);
                    cartesian[i]=1;
                    // the projection of each cartesian basis vector onto the plane perp to the tangent vector
                    double projx = (cartesian[0] - (cartesian[0]*txs + cartesian[1]*tys + cartesian[2]*tzs)*txs);
                    double projy = (cartesian[1] - (cartesian[0]*txs + cartesian[1]*tys + cartesian[2]*tzs)*txs);
                    double projz = (cartesian[2] - (cartesian[0]*txs + cartesian[1]*tys + cartesian[2]*tzs)*txs);
                    double tempnorm = sqrt(projx*projx+projy*projy+projz*projz);
                    // Ill actually use the perp with the biggest norm - one of them could be colinear with the tangent vector, this ferrets that one out
                    if(tempnorm>maxnorm)
                    {
                        e1x = projx/tempnorm;
                        e1y = projy/tempnorm;
                        e1z = projz/tempnorm;
                    }
                }
                // and take a cross product to get the third in our triple
                e2x = tys*e1z - tzs*e1y;
                e2y = tzs*e1x - txs*e1z;
                e2z = txs*e1y - tys*e1x;
                // okay we have our directions in the plane we want to perfrom the line minimisation in. Time to do it
                // the point
                gsl_vector* v = gsl_vector_alloc (3);
                gsl_vector_set (v, 0, testx);
                gsl_vector_set (v, 1, testy);
                gsl_vector_set (v, 2, testz);
                // one vector in the plane we with to minimize in
                gsl_vector* f = gsl_vector_alloc (3);
                gsl_vector_set (f, 0, e1x);
                gsl_vector_set (f, 1, e1y);
                gsl_vector_set (f, 2, e1z);
                // and the other
                gsl_vector* b = gsl_vector_alloc (3);
                gsl_vector_set (b, 0, e2x);
                gsl_vector_set (b, 1, e2y);
                gsl_vector_set (b, 2, e2z);
                gsl_vector* minimum = gsl_vector_alloc (2);
                gsl_vector_set (minimum, 0, 0);
                gsl_vector_set (minimum, 1, 0);

                // set up the structure to pass into the GSL minimizer. The naming is all from old code
                struct parameters params; struct parameters* pparams = &params;
                pparams->ucvmag=&interpolatedmagb;
                pparams->v = v; pparams->f = f;pparams->b=b;
                // some initial values
                gsl_multimin_function F;
                F.n=2;
                F.f = &my_minimisation_function;
                F.params = (void*) pparams;
                gsl_vector* stepsize = gsl_vector_alloc (2);
                gsl_vector_set (stepsize, 0, 0.5);
                gsl_vector_set (stepsize, 1, 0.5);
                gsl_multimin_fminimizer_set (minimizerstate, &F, minimum, stepsize);
                int iter=0;
                int status =0;
                double minimizersize=0;

                do
                {
                    iter++;
                    status = gsl_multimin_fminimizer_iterate(minimizerstate);

                    if (status)
                        break;

                    minimizersize = gsl_multimin_fminimizer_size (minimizerstate);
                    status = gsl_multimin_test_size (minimizersize, 1e-10);

                }

                while (status == GSL_CONTINUE && iter < 1000);

                gsl_vector_scale(f,gsl_vector_get(minimizerstate->x, 0));
                gsl_vector_scale(b,gsl_vector_get(minimizerstate->x, 1));
                gsl_vector_add(f,b);
                gsl_vector_add(v,f);
                knotcurves[c].knotcurve[s].xcoord = gsl_vector_get(v, 0);
                knotcurves[c].knotcurve[s].ycoord= gsl_vector_get(v, 1);
                knotcurves[c].knotcurve[s].zcoord= gsl_vector_get(v, 2);

                knotcurves[c].knotcurve[s].tx = txs;
                knotcurves[c].knotcurve[s].ty = tys;
                knotcurves[c].knotcurve[s].tz = tzs;

                knotcurves[c].knotcurve[s].gradmagbx =e1x;
                knotcurves[c].knotcurve[s].gradmagby = e1y;
                knotcurves[c].knotcurve[s].gradmagbz = e1z;

                knotcurves[c].knotcurve[s].gradmagbperpx = e2x;
                knotcurves[c].knotcurve[s].gradmagbperpy = e2y;
                knotcurves[c].knotcurve[s].gradmagbperpz = e2z;

                gsl_vector_free(v);
                gsl_vector_free(f);
                gsl_vector_free(b);
                gsl_vector_free(stepsize);

                // mark the marked array with this curve
                marked[pt((int) knotcurves[c].knotcurve[s].xcoord,(int) knotcurves[c].knotcurve[s].ycoord,(int) knotcurves[c].knotcurve[s].zcoord)]=-1;

                // okay we got our new point. Now the question is, when should we terminate this process? Do it by whether we are close to the start point after some nontrivial number of steps
                double xdiff = knotcurves[c].knotcurve[0].xcoord - knotcurves[c].knotcurve[s].xcoord;     //distance from start/end point
                double ydiff = knotcurves[c].knotcurve[0].ycoord - knotcurves[c].knotcurve[s].ycoord;
                double zdiff = knotcurves[c].knotcurve[0].zcoord - knotcurves[c].knotcurve[s].zcoord;
                if( (sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff) <3 && s > 30 )||s>2000) finish = true;
                s++;
            }
            // for each curve, I will extract the subset of our data such that magb < threshold which contains the curve. We we 
            // only look in this connected component once
            // construct a tube around the knot, to use as an excluded region if we searching for multiple components.

            bool stillboundaryleft = true;
            while(stillboundaryleft)
            {
                // the marked array has the following values
                // 0 - not evaluated 
                // -1 - a boundary, to be grown 
                // -2 - the interrior, already grown
                // -3 - a temporary state, marked as a boundary during the update
                // positive numbers - layers of shells already marked
                for(int m=1;m<Lz-1;m++)
                {
                    for(int l=1; l<Ly-1; l++)
                    {
                        for(int k=1; k<Lx-1; k++)   //Central difference
                        {
                            int n = pt(k,l,m);
                            if(marked[n] ==-1)
                            {

                                for(int minc=-1;minc<=1;minc++)
                                {
                                    for(int linc=-1; linc<=1; linc++)
                                    {
                                        for(int kinc=-1; kinc<=1; kinc++)   //Central difference
                                        {
                                            int neighboringn = pt(k+kinc,l+linc,m+minc);
                                            if(marked[neighboringn] == 0 && magb[neighboringn] < threshold) marked[neighboringn] = -3;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                for(int n = 0; n<LL;n++)
                {
                    if(marked[n]==-1){marked[n] =-2;}
                    if(marked[n]==-3){marked[n] =-1;}
                }
                stillboundaryleft = false;
                for(int n = 0; n<LL;n++)
                {
                    if(marked[n]==-1) stillboundaryleft =true;
                }
            }
        }
    // temp
    char vtk_data[200];
    sprintf(vtk_data,"marked_%d.vtk",c);
    ofstream Bout (vtk_data);
    Bout << "# vtk DataFile Version 3.0\nKnot\nASCII\nDATASET STRUCTURED_POINTS\n";
    Bout << "DIMENSIONS " << Lx << ' ' << Ly << ' ' << Lz << '\n';
    Bout << "ORIGIN 0 0 0" << '\n';
    Bout << "SPACING " << 1 << ' ' << 1 << ' ' << 1 << '\n';
    Bout << "POINT_DATA " << LL << '\n';
    Bout << "SCALARS marked float\nLOOKUP_TABLE default\n";
    for(int k=0; k<Lz; k++)
    {
        for(int j=0; j<Ly; j++)
        {
            for(int i=0; i<Lx; i++)
            {
                int n = pt(i,j,k);
                Bout << marked[n] << '\n';
            }
        }
    }
    Bout.close();
    //
        c++;
    }
}

// some little functions which go back and forth from point to index
int pt(const int k,const  int l,const  int m)       //convert i,j,k to single index
{
    return (m*Ly*Lx+l*Lx+k);
}

double my_minimisation_function(const gsl_vector* minimum, void* params)
{
    struct parameters* myparameters = (struct parameters *) params;
    likely::TriCubicInterpolator* interpolateducvmag = myparameters->ucvmag;
    gsl_vector* tempf = gsl_vector_alloc (3);
    gsl_vector* tempv = gsl_vector_alloc (3);
    gsl_vector* tempb = gsl_vector_alloc (3);
    gsl_vector_memcpy (tempf,myparameters->f);
    gsl_vector_memcpy (tempv,myparameters->v);
    gsl_vector_memcpy (tempb,myparameters->b);

    // s gives us how much of f to add to p
    gsl_vector_scale(tempf,gsl_vector_get (minimum, 0));
    gsl_vector_scale(tempb,gsl_vector_get (minimum, 1));
    gsl_vector_add(tempf,tempb);
    gsl_vector_add(tempv,tempf);
    double px = gsl_vector_get(tempv, 0);
    double py = gsl_vector_get(tempv, 1);
    double pz = gsl_vector_get(tempv, 2);
    gsl_vector_free(tempf);
    gsl_vector_free(tempv);
    gsl_vector_free(tempb);

    double value = ((*interpolateducvmag)(px,py,pz));
    return value;
}
