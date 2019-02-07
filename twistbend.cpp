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

    int n =0;
    int step=0;
    startconfig(n, nx,ny,nz,px,py,pz);
    cout << "starting simulation" << endl;

#pragma omp parallel default(none) shared(nx,ny,nz,px,py,pz,hx,hy,hz,hpx,hpy,hpz, bx,by,bz,bmag, tx,ty,tz,n,step,cout)
    {
        while(n<=Nmax)
        {
#pragma omp single
            {
                if (n%1000==0)
                {
                    cout << "timestep " << n << endl;
                }

                if (step==stepskip || step==0)
                {
                    cout << "writing VTK files at timestep " << n << endl;
                    computeBendAndCurlofCirculation(n,nx,ny,nz,bx,by,bz,bmag,tx,ty,tz);
                    writeVTKfiles(n,nx,ny,nz,px,py,pz,bx,by,bz,tx,ty,tz);       // output VTK files for use with ParaView
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
            double RR2 = 0.25*Lz;  // scale
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

void computeBendAndCurlofCirculation(const int n, const double* nx,const double* ny,const double* nz, double* bx, double* by, double* bz, double* bmag, double* tx, double* ty, double* tz)
{
    // arrays storing the circulaiton, which we will take the curl of
    double* cx=new double[LL];
    double* cy=new double[LL];
    double* cz=new double[LL];

    int j,k,l,m;
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

//void FindBendZeros(double *magb,double* circulationx,double* circulationy,double* circulationz, vector<knotcurve>& knotcurves,double n, gsl_multimin_fminimizer* minimizerstate)
//{
//    // first thing, clear the knotcurve object before we begin writing a new one
//    knotcurves.clear(); //empty vector with knot curve points
//
//    // initialise the tricubic interpolator for ucvmag
//    likely::TriCubicInterpolator interpolatedb(magb, 1, Lx,Ly,Lz);
//
//    int c =0;
//    bool knotexists = true;
//    vector<int> marked(Lx*Ly*Lz,0);
//
//    while(knotexists)
//    {
//        double   bmax = -1.0; // should always be +ve, so setting it to an initially -ve # means it always gets written to once.
//        int kmax,lmax,mmax;
//        for(int k=0;k<Lx;k++)
//        {
//            for(int l=0; l<Ly; l++)
//            {
//                for(int m=0; m<Lz; m++)   //Central difference
//                {
//                    int n = pt(k,l,m);
//                    if( magb[n] > bmax && marked[n]==0)
//                    {
//                        bmax = magb[n];
//                        kmax = k;
//                        lmax = l;
//                        mmax=m;
//
//                    }
//                }
//            }
//        }
//
//        if(bmax<0.45) knotexists = false;
//
//        if(knotexists)
//        {
//            knotcurves.push_back(knotcurve() );
//            knotcurves[c].knotcurve.push_back(knotpoint());
//            knotcurves[c].knotcurve[0].xcoord=x(kmax);
//            knotcurves[c].knotcurve[0].ycoord=y(lmax);
//            knotcurves[c].knotcurve[0].zcoord=z(mmax);
//
//            int idwn,jdwn,kdwn, modidwn, modjdwn, modkdwn,m,iinc,jinc,kinc;
//            double ucvxs, ucvys, ucvzs, graducvx, graducvy, graducvz, prefactor, xd, yd ,zd, fx, fy, fz, xdiff, ydiff, zdiff;
//            int s=1;
//            bool finish=false;
//            // we will discard the first few points from the knot, using this flag
//            bool burnin=true;
//            // if the curve we are tracing terminates at a boundary, for now we will just strike it from the record, entering a cleanup mode that marks the broken curve as "do not touch"
//            int boundaryhits = 0;
//            // don't allow hits all over the boundary
//            int cooldowntimer = 0;
//            /*calculate local direction of grad u x grad v (the tangent to the knot curve) at point s-1, then move to point s by moving along tangent + unit confinement force*/
//            while (finish==false)
//            {
//
//                /**Find nearest gridpoint**/
//                idwn = (int) ((knotcurves[c].knotcurve[s-1].xcoord/h) - 0.5 + Nx/2.0);
//                jdwn = (int) ((knotcurves[c].knotcurve[s-1].ycoord/h) - 0.5 + Ny/2.0);
//                kdwn = (int) ((knotcurves[c].knotcurve[s-1].zcoord/h) - 0.5 + Nz/2.0);
//                // idwn etc can be off the actual grid , into "ghost" grids around the real one. this is useful for knotcurve tracing over periodic boundaries
//                // but we also need the corresponding real grid positions!
//                modidwn = circularmod(idwn,Nx);
//                modjdwn = circularmod(jdwn,Ny);
//                modkdwn = circularmod(kdwn,Nz);
//
//                // if we have hit a boundary, dont go off grid - rather, set a cleanup flag and stick to the grid egde
//                if(cooldowntimer>0) {cooldowntimer--;}
//                if(BoundaryType==ALLREFLECTING && cooldowntimer ==0)
//                {
//                    if(idwn <=0){idwn=0; boundaryhits++; cooldowntimer=5;}
//                    if(idwn >=Nx-1){idwn=Nx-1; boundaryhits++;cooldowntimer=5;}
//                    if(jdwn <=0){jdwn=0; boundaryhits++;cooldowntimer=5;}
//                    if(jdwn >=Ny-1){jdwn=Ny-1; boundaryhits++;cooldowntimer=5;}
//                    if(kdwn <=0){kdwn=0; boundaryhits++;cooldowntimer=5;}
//                    if(kdwn >=Nz-1){kdwn=Nz-1; boundaryhits++;cooldowntimer=5;}
//                }
//                if(BoundaryType==ZPERIODIC && cooldowntimer==0)
//                {
//                    if(idwn <=0){idwn=0; boundaryhits++;cooldowntimer=5;}
//                    if(idwn >Nx-1){idwn=Nx-1; boundaryhits++;cooldowntimer=5;}
//                    if(jdwn <=0){jdwn=0; boundaryhits++;cooldowntimer=5;}
//                    if(jdwn >Ny-1){jdwn=Ny-1; boundaryhits++;cooldowntimer=5;}
//                }
//
//                ucvxs=0;
//                ucvys=0;
//                ucvzs=0;
//                /*curve to gridpoint down distance*/
//                xd = (knotcurves[c].knotcurve[s-1].xcoord - x(idwn,griddata))/h;
//                yd = (knotcurves[c].knotcurve[s-1].ycoord - y(jdwn,griddata))/h;
//                zd = (knotcurves[c].knotcurve[s-1].zcoord - z(kdwn,griddata))/h;
//                for(m=0;m<8;m++)  //linear interpolation from 8 nearest neighbours
//                {
//                    /* Work out increments*/
//                    iinc = m%2;
//                    jinc = (m/2)%2;
//                    kinc = (m/4)%2;
//                    /*Loop over nearest points*/
//                    i = gridinc(modidwn, iinc, Nx,0);
//                    j = gridinc(modjdwn, jinc, Ny,1);
//                    k = gridinc(modkdwn,kinc, Nz,2);
//                    prefactor = (1-iinc + pow(-1,1+iinc)*xd)*(1-jinc + pow(-1,1+jinc)*yd)*(1-kinc + pow(-1,1+kinc)*zd);
//                    /*interpolate grad u x grad v over nearest points*/
//                    ucvxs += prefactor*ucvx[pt(i,j,k,griddata)];
//                    ucvys += prefactor*ucvy[pt(i,j,k,griddata)];
//                    ucvzs += prefactor*ucvz[pt(i,j,k,griddata)];
//                }
//                double norm = sqrt(ucvxs*ucvxs + ucvys*ucvys + ucvzs*ucvzs);
//                ucvxs = ucvxs/norm; //normalise
//                ucvys = ucvys/norm; //normalise
//                ucvzs = ucvzs/norm; //normalise
//
//                // if we have hit a boundary, we want to back up along the curve instead
//                if(boundaryhits==1)
//                {
//                    ucvxs *= -1 ;
//                    ucvys *= -1 ;
//                    ucvzs *= -1;
//                }
//
//                // okay we have our first guess, move forward in this direction
//                // we actually want to walk in the direction gradv cross gradu - that should be our +ve tangent,
//                // so that the rotation sense of the curve is positive. Get this we - signs below.
//                double testx = knotcurves[c].knotcurve[s-1].xcoord - h*ucvxs;
//                double testy = knotcurves[c].knotcurve[s-1].ycoord - h*ucvys;
//                double testz = knotcurves[c].knotcurve[s-1].zcoord - h*ucvzs;
//
//                // now get the grad at this point
//                idwn = (int) ((testx/h) - 0.5 + Nx/2.0);
//                jdwn = (int) ((testy/h) - 0.5 + Ny/2.0);
//                kdwn = (int) ((testz/h) - 0.5 + Nz/2.0);
//                modidwn = circularmod(idwn,Nx);
//                modjdwn = circularmod(jdwn,Ny);
//                modkdwn = circularmod(kdwn,Nz);
//                graducvx=0;
//                graducvy=0;
//                graducvz=0;
//                /*curve to gridpoint down distance*/
//                xd = (testx - x(idwn,griddata))/h;
//                yd = (testy - y(jdwn,griddata))/h;
//                zd = (testz - z(kdwn,griddata))/h;
//                for(m=0;m<8;m++)  //linear interpolation from 8 nearest neighbours
//                {
//                    /* Work out increments*/
//                    iinc = m%2;
//                    jinc = (m/2)%2;
//                    kinc = (m/4)%2;
//                    /*Loop over nearest points*/
//                    i = gridinc(modidwn, iinc, Nx,0);
//                    j = gridinc(modjdwn, jinc, Ny,1);
//                    k = gridinc(modkdwn,kinc, Nz,2);
//                    prefactor = (1-iinc + pow(-1,1+iinc)*xd)*(1-jinc + pow(-1,1+jinc)*yd)*(1-kinc + pow(-1,1+kinc)*zd);
//                    /*interpolate gradients of |grad u x grad v|*/
//                    graducvx += prefactor*(sqrt(ucvx[pt(gridinc(i,1,Nx,0),j,k,griddata)]*ucvx[pt(gridinc(i,1,Nx,0),j,k,griddata)] + ucvy[pt(gridinc(i,1,Nx,0),j,k,griddata)]*ucvy[pt(gridinc(i,1,Nx,0),j,k,griddata)] + ucvz[pt(gridinc(i,1,Nx,0),j,k,griddata)]*ucvz[pt(gridinc(i,1,Nx,0),j,k,griddata)]) - sqrt(ucvx[pt(gridinc(i,-1,Nx,0),j,k,griddata)]*ucvx[pt(gridinc(i,-1,Nx,0),j,k,griddata)] + ucvy[pt(gridinc(i,-1,Nx,0),j,k,griddata)]*ucvy[pt(gridinc(i,-1,Nx,0),j,k,griddata)] + ucvz[pt(gridinc(i,-1,Nx,0),j,k,griddata)]*ucvz[pt(gridinc(i,-1,Nx,0),j,k,griddata)]))/(2*h);
//                    graducvy += prefactor*(sqrt(ucvx[pt(i,gridinc(j,1,Ny,1),k,griddata)]*ucvx[pt(i,gridinc(j,1,Ny,1),k,griddata)] + ucvy[pt(i,gridinc(j,1,Ny,1),k,griddata)]*ucvy[pt(i,gridinc(j,1,Ny,1),k,griddata)] + ucvz[pt(i,gridinc(j,1,Ny,1),k,griddata)]*ucvz[pt(i,gridinc(j,1,Ny,1),k,griddata)]) - sqrt(ucvx[pt(i,gridinc(j,-1,Ny,1),k,griddata)]*ucvx[pt(i,gridinc(j,-1,Ny,1),k,griddata)] + ucvy[pt(i,gridinc(j,-1,Ny,1),k,griddata)]*ucvy[pt(i,gridinc(j,-1,Ny,1),k,griddata)] + ucvz[pt(i,gridinc(j,-1,Ny,1),k,griddata)]*ucvz[pt(i,gridinc(j,-1,Ny,1),k,griddata)]))/(2*h);
//                    graducvz += prefactor*(sqrt(ucvx[pt(i,j,gridinc(k,1,Nz,2),griddata)]*ucvx[pt(i,j,gridinc(k,1,Nz,2),griddata)] + ucvy[pt(i,j,gridinc(k,1,Nz,2),griddata)]*ucvy[pt(i,j,gridinc(k,1,Nz,2),griddata)] + ucvz[pt(i,j,gridinc(k,1,Nz,2),griddata)]*ucvz[pt(i,j,gridinc(k,1,Nz,2),griddata)]) - sqrt(ucvx[pt(i,j,gridinc(k,-1,Nz,2),griddata)]*ucvx[pt(i,j,gridinc(k,-1,Nz,2),griddata)] + ucvy[pt(i,j,gridinc(k,-1,Nz,2),griddata)]*ucvy[pt(i,j,gridinc(k,-1,Nz,2),griddata)] + ucvz[pt(i,j,gridinc(k,-1,Nz,2),griddata)]*ucvz[pt(i,j,gridinc(k,-1,Nz,2),griddata)]))/(2*h);
//
//                }
//                knotcurves[c].knotcurve.push_back(knotpoint());
//                // one of the vectors in the plane we wish to perfrom our minimisation in
//                fx = (graducvx - (graducvx*ucvxs + graducvy*ucvys + graducvz*ucvzs)*ucvxs);
//                fy = (graducvy - (graducvx*ucvxs + graducvy*ucvys + graducvz*ucvzs)*ucvys);
//                fz = (graducvz - (graducvx*ucvxs + graducvy*ucvys + graducvz*ucvzs)*ucvzs);
//                norm = sqrt(fx*fx + fy*fy + fz*fz);
//                fx = fx/norm;
//                fy = fy/norm;
//                fz = fz/norm;
//
//                // okay we have our direction to perfrom the line minimisation in
//                // the point
//                gsl_vector* v = gsl_vector_alloc (3);
//                gsl_vector_set (v, 0, testx);
//                gsl_vector_set (v, 1, testy);
//                gsl_vector_set (v, 2, testz);
//                // one vector in the plane we with to minimize in
//                gsl_vector* f = gsl_vector_alloc (3);
//                gsl_vector_set (f, 0, fx);
//                gsl_vector_set (f, 1, fy);
//                gsl_vector_set (f, 2, fz);
//                // the ucv vector
//                gsl_vector* ucv = gsl_vector_alloc (3);
//                gsl_vector_set (ucv, 0, ucvxs);
//                gsl_vector_set (ucv, 1, ucvys);
//                gsl_vector_set (ucv, 2, ucvzs);
//                // take a cross product to get the other vector in the plane
//                gsl_vector* b = gsl_vector_alloc (3);
//                cross_product(f,ucv,b);
//                // initial conditions
//                gsl_vector* minimum = gsl_vector_alloc (2);
//                gsl_vector_set (minimum, 0, 0);
//                gsl_vector_set (minimum, 1, 0);
//                struct parameters params; struct parameters* pparams = &params;
//                pparams->ucvmag=&interpolateducvmag;
//                pparams->v = v; pparams->f = f;pparams->b=b;
//                pparams->mygriddata = griddata;
//                // some initial values
//                gsl_multimin_function F;
//                F.n=2;
//                F.f = &my_f;
//                F.params = (void*) pparams;
//                gsl_vector* stepsize = gsl_vector_alloc (2);
//                gsl_vector_set (stepsize, 0, lambda/(8*M_PI));
//                gsl_vector_set (stepsize, 1, lambda/(8*M_PI));
//                gsl_multimin_fminimizer_set (minimizerstate, &F, minimum, stepsize);
//
//                int iter=0;
//                int status =0;
//                double minimizersize=0;
//                do
//                {
//                    iter++;
//                    status = gsl_multimin_fminimizer_iterate(minimizerstate);
//
//                    if (status)
//                        break;
//
//                    minimizersize = gsl_multimin_fminimizer_size (minimizerstate);
//                    status = gsl_multimin_test_size (minimizersize, 1e-2);
//
//                }
//                while (status == GSL_CONTINUE && iter < 5000);
//
//
//                gsl_vector_scale(f,gsl_vector_get(minimizerstate->x, 0));
//                gsl_vector_scale(b,gsl_vector_get(minimizerstate->x, 1));
//                gsl_vector_add(f,b);
//                gsl_vector_add(v,f);
//                knotcurves[c].knotcurve[s].xcoord = gsl_vector_get(v, 0);
//                knotcurves[c].knotcurve[s].ycoord= gsl_vector_get(v, 1);
//                knotcurves[c].knotcurve[s].zcoord= gsl_vector_get(v, 2);
//
//                gsl_vector_free(v);
//                gsl_vector_free(f);
//                gsl_vector_free(b);
//                gsl_vector_free(ucv);
//                gsl_vector_free(stepsize);
//
//                xdiff = knotcurves[c].knotcurve[0].xcoord - knotcurves[c].knotcurve[s].xcoord;     //distance from start/end point
//                ydiff = knotcurves[c].knotcurve[0].ycoord - knotcurves[c].knotcurve[s].ycoord;
//                zdiff = knotcurves[c].knotcurve[0].zcoord - knotcurves[c].knotcurve[s].zcoord;
//
//                if( (boundaryhits==0 && sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff) <h  && s > 10 ) || boundaryhits>1 ||s>5000) finish = true;
//
//                // okay, we just added a point in position s in the vector
//                // if we have a few points in the vector, discard the first few and restart the whole thing - burn it in
//                int newstartingposition =20;
//                if(s==newstartingposition && burnin && (boundaryhits==0))
//                {
//                    knotcurves[c].knotcurve.erase(knotcurves[c].knotcurve.begin(),knotcurves[c].knotcurve.begin()+newstartingposition);
//                    s =0;
//                    burnin =false;
//                }
//
//                s++;
//            }
//            int NP = knotcurves[c].knotcurve.size();  //store number of points in knot curve
//
//            // construct a tube around the knot, to use as an excluded region if we searching for multiple components.
//            double radius = 3;
//            int numiterations = (int)(radius/griddata.h);
//            for(int s=0; s<NP; s++)
//            {
//                int icentral = (int) ((knotcurves[c].knotcurve[s].xcoord/h) - 0.5 + Nx/2.0);
//                int jcentral = (int) ((knotcurves[c].knotcurve[s].ycoord/h) - 0.5 + Ny/2.0);
//                int kcentral = (int) ((knotcurves[c].knotcurve[s].zcoord/h) - 0.5 + Nz/2.0);
//                // construct a ball of radius "radius" around each point in the knotcurve object. we circumscribe it in a cube which is then looped over
//                for(int i =-numiterations;i<=numiterations;i++)
//               {
//                    for(int j=-numiterations ;j<=numiterations;j++)
//                    {
//                        for(int k =-numiterations;k<=numiterations;k++)
//                        {
//                            int modi = circularmod(i+icentral,Nx);
//                            int modj = circularmod(j+jcentral,Ny);
//                            int modk = circularmod(k+kcentral,Nz);
//                            int n = pt(modi,modj,modk,griddata);
//
//                            double dxsq = (x(i+icentral,griddata)-x(icentral,griddata))*(x(i+icentral,griddata)-x(icentral,griddata));
//                            double dysq = (y(j+jcentral,griddata)-y(jcentral,griddata))*(y(j+jcentral,griddata)-y(jcentral,griddata));
//                            double dzsq = (z(k+kcentral,griddata)-z(kcentral,griddata))*(z(k+kcentral,griddata)-z(kcentral,griddata));
//
//                            double r = sqrt(dxsq + dysq + dzsq);
//
//                            if(r < radius )
//                            {
//                                marked[n]=1;
//                            }
//                        }
//                    }
//                }
//            }
//
//            // now comes a lot of curve analysis. but for now Im not going to do this on the boundary curves. 
//            if(boundaryhits==0)
//            {
//
//                /*******Vertex averaging*********/
//
//                double totlength, dl, dx,dy,dz;
//                for(i=0;i<3;i++)   //repeat a couple of times because of end point
//                {
//                    totlength=0;
//                    for(s=0; s<NP; s++)   //Work out total length of curve
//                    {
//                        dx = knotcurves[c].knotcurve[incp(s,1,NP)].xcoord - knotcurves[c].knotcurve[s].xcoord;
//                        dy = knotcurves[c].knotcurve[incp(s,1,NP)].ycoord - knotcurves[c].knotcurve[s].ycoord;
//                        dz = knotcurves[c].knotcurve[incp(s,1,NP)].zcoord - knotcurves[c].knotcurve[s].zcoord;
//                        totlength += sqrt(dx*dx + dy*dy + dz*dz);
//                    }
//                    dl = totlength/NP;
//                    for(s=0; s<NP; s++)    //Move points to have spacing dl
//                    {
//                        dx = knotcurves[c].knotcurve[incp(s,1,NP)].xcoord - knotcurves[c].knotcurve[s].xcoord;
//                        dy = knotcurves[c].knotcurve[incp(s,1,NP)].ycoord - knotcurves[c].knotcurve[s].ycoord;
//                        dz = knotcurves[c].knotcurve[incp(s,1,NP)].zcoord - knotcurves[c].knotcurve[s].zcoord;
//                        double norm = sqrt(dx*dx + dy*dy + dz*dz);
//                        knotcurves[c].knotcurve[incp(s,1,NP)].xcoord = knotcurves[c].knotcurve[s].xcoord + dl*dx/norm;
//                        knotcurves[c].knotcurve[incp(s,1,NP)].ycoord = knotcurves[c].knotcurve[s].ycoord + dl*dy/norm;
//                        knotcurves[c].knotcurve[incp(s,1,NP)].zcoord = knotcurves[c].knotcurve[s].zcoord + dl*dz/norm;
//                    }
//                }
//
//                /*************Curve Smoothing*******************/
//                vector<double> coord(NP);
//                gsl_fft_real_wavetable * real;
//                gsl_fft_halfcomplex_wavetable * hc;
//                gsl_fft_real_workspace * work;
//                work = gsl_fft_real_workspace_alloc (NP);
//                real = gsl_fft_real_wavetable_alloc (NP);
//                hc = gsl_fft_halfcomplex_wavetable_alloc (NP);
//                for(j=1; j<4; j++)
//                {
//                    switch(j)
//                    {
//                        case 1 :
//                            for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].xcoord ; break;
//                        case 2 :
//                            for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].ycoord ; break;
//                        case 3 :
//                            for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].zcoord ; break;
//                    }
//                    double* data = coord.data();
//                    // take the fft
//                    gsl_fft_real_transform (data, 1, NP, real, work);
//                    // 21/11/2016: make our low pass filter. To apply our filter. we should sample frequencies fn = n/Delta N , n = -N/2 ... N/2
//                    // this is discretizing the nyquist interval, with extreme frequency ~1/2Delta.
//                    // to cut out the frequencies of grid fluctuation size and larger we need a lengthscale Delta to
//                    // plug in above. im doing a rough length calc below, this might be overkill.
//                    // at the moment its just a hard filter, we can choose others though.
//                    // compute a rough length to set scale
//                    double filter;
//                    const double cutoff = 2*M_PI*(totlength/(6*lambda));
//                    for (i = 0; i < NP; ++i)
//                    {
//                        filter = 1/sqrt(1+pow((i/cutoff),8));
//                        data[i] *= filter;
//                    };
//                    // transform back
//                    gsl_fft_halfcomplex_inverse (data, 1, NP, hc, work);
//                    switch(j)
//                    {
//                        case 1 :
//                            for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].xcoord = coord[i] ; break;
//                        case 2 :
//                            for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].ycoord = coord[i] ; break;
//                        case 3 :
//                            for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].zcoord = coord[i] ; break;
//                    }
//                }
//
//
//
//                /******************Interpolate direction of grad u for twist calc*******/
//                /**Find nearest gridpoint**/
//                double dxu, dyu, dzu, dxup, dyup, dzup;
//                for(s=0; s<NP; s++)
//                {
//                    idwn = (int) ((knotcurves[c].knotcurve[s].xcoord/h) - 0.5 + Nx/2.0);
//                    jdwn = (int) ((knotcurves[c].knotcurve[s].ycoord/h) - 0.5 + Ny/2.0);
//                    kdwn = (int) ((knotcurves[c].knotcurve[s].zcoord/h) - 0.5 + Nz/2.0);
//                    modidwn = circularmod(idwn,Nx);
//                    modjdwn = circularmod(jdwn,Ny);
//                    modkdwn = circularmod(kdwn,Nz);
//                    if((BoundaryType==ALLREFLECTING) && (idwn<0 || jdwn<0 || kdwn<0 || idwn > Nx-1 || jdwn > Ny-1 || kdwn > Nz-1)) break;
//                    if((BoundaryType==ZPERIODIC) && (idwn<0 || jdwn<0 || idwn > Nx-1 || jdwn > Ny-1 )) break;
//                    dxu=0;
//                    dyu=0;
//                    dzu=0;
//                    /*curve to gridpoint down distance*/
//                    xd = (knotcurves[c].knotcurve[s].xcoord - x(idwn,griddata))/h;
//                    yd = (knotcurves[c].knotcurve[s].ycoord - y(jdwn,griddata))/h;
//                    zd = (knotcurves[c].knotcurve[s].zcoord - z(kdwn,griddata))/h;
//                    for(m=0;m<8;m++)  //linear interpolation of 8 NNs
//                    {
//                        /* Work out increments*/
//                        iinc = m%2;
//                        jinc = (m/2)%2;
//                        kinc = (m/4)%2;
//                        /*Loop over nearest points*/
//                        i = gridinc(modidwn, iinc, Nx,0);
//                        j = gridinc(modjdwn, jinc, Ny,1);
//                        k = gridinc(modkdwn,kinc, Nz,2);
//                        prefactor = (1-iinc + pow(-1,1+iinc)*xd)*(1-jinc + pow(-1,1+jinc)*yd)*(1-kinc + pow(-1,1+kinc)*zd);   //terms of the form (1-xd)(1-yd)zd etc. (interpolation coefficient)
//                        /*interpolate grad u over nearest points*/
//                        dxu += prefactor*0.5*(u[pt(gridinc(i,1,Nx,0),j,k,griddata)] -  u[pt(gridinc(i,-1,Nx,0),j,k,griddata)])/h;  //central diff
//                        dyu += prefactor*0.5*(u[pt(i,gridinc(j,1,Ny,1),k,griddata)] -  u[pt(i,gridinc(j,-1,Ny,1),k,griddata)])/h;
//                        dzu += prefactor*0.5*(u[pt(i,j,gridinc(k,1,Nz,2),griddata)] -  u[pt(i,j,gridinc(k,-1,Nz,2),griddata)])/h;
//                    }
//                    //project du onto perp of tangent direction first
//                    dx = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].xcoord - knotcurves[c].knotcurve[incp(s,-1,NP)].xcoord);   //central diff as a is defined on the points
//                    dy = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].ycoord - knotcurves[c].knotcurve[incp(s,-1,NP)].ycoord);
//                    dz = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].zcoord - knotcurves[c].knotcurve[incp(s,-1,NP)].zcoord);
//                    dxup = dxu - (dxu*dx + dyu*dy + dzu*dz)*dx/(dx*dx+dy*dy+dz*dz);               //Grad u_j * (delta_ij - t_i t_j)
//                    dyup = dyu - (dxu*dx + dyu*dy + dzu*dz)*dy/(dx*dx+dy*dy+dz*dz);
//                    dzup = dzu - (dxu*dx + dyu*dy + dzu*dz)*dz/(dx*dx+dy*dy+dz*dz);
//                    /*Vector a is the normalised gradient of u, should point in direction of max u perp to t*/
//                    double norm = sqrt(dxup*dxup+dyup*dyup+dzup*dzup);
//                    knotcurves[c].knotcurve[s].ax = dxup/norm;
//                    knotcurves[c].knotcurve[s].ay = dyup/norm;
//                    knotcurves[c].knotcurve[s].az = dzup/norm;
//                }
//
//                for(j=1; j<4; j++)
//                {
//                    switch(j)
//                    {
//                        case 1 :
//                            for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].ax ; break;
//                        case 2 :
//                            for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].ay ; break;
//                        case 3 :
//                            for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].az ; break;
//                    }
//                    double* data = coord.data();
//                    // take the fft
//                    gsl_fft_real_transform (data, 1, NP, real, work);
//                    // 21/11/2016: make our low pass filter. To apply our filter. we should sample frequencies fn = n/Delta N , n = -N/2 ... N/2
//                    // this is discretizing the nyquist interval, with extreme frequency ~1/2Delta.
//                    // to cut out the frequencies of grid fluctuation size and larger we need a lengthscale Delta to
//                    // plug in above. im doing a rough length calc below, this might be overkill.
//                    // at the moment its just a hard filter, we can choose others though.
//                    // compute a rough length to set scale
//                    double filter;
//                    const double cutoff = 2*M_PI*(totlength/(6*lambda));
//                    for (i = 0; i < NP; ++i)
//                    {
//                        filter = 1/sqrt(1+pow((i/cutoff),8));
//                        data[i] *= filter;
//                    };
//                    // transform back
//                    gsl_fft_halfcomplex_inverse (data, 1, NP, hc, work);
//                    switch(j)
//                    {
//                        case 1 :
//                            for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].ax= coord[i] ; break;
//                        case 2 :
//                            for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].ay= coord[i] ; break;
//                        case 3 :
//                            for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].az = coord[i] ; break;
//                    }
//                }
//                gsl_fft_real_wavetable_free (real);
//                gsl_fft_halfcomplex_wavetable_free (hc);
//                gsl_fft_real_workspace_free (work);
//
//
//                // CURVE GEOMETRY - get curvatures, torsions, frennet serret frame
//
//
//                NP = knotcurves[c].knotcurve.size();
//                for(s=0; s<NP; s++)
//                {
//                    // forward difference on the tangents
//                    double dx = (knotcurves[c].knotcurve[incp(s,1,NP)].xcoord - knotcurves[c].knotcurve[incp(s,0,NP)].xcoord);
//                    double dy = (knotcurves[c].knotcurve[incp(s,1,NP)].ycoord - knotcurves[c].knotcurve[incp(s,0,NP)].ycoord);
//                    double dz = (knotcurves[c].knotcurve[incp(s,1,NP)].zcoord - knotcurves[c].knotcurve[incp(s,0,NP)].zcoord);
//                    double deltas = sqrt(dx*dx+dy*dy+dz*dz);
//                    knotcurves[c].knotcurve[s].tx = dx/(deltas);
//                    knotcurves[c].knotcurve[s].ty = dy/(deltas);
//                    knotcurves[c].knotcurve[s].tz = dz/(deltas);
//                    knotcurves[c].knotcurve[s].length = deltas;
//                    knotcurves[c].length +=deltas;
//                }
//                for(s=0; s<NP; s++)
//                {
//                    // backwards diff for the normals, amounting to a central diff overall
//                    double nx = 2.0*(knotcurves[c].knotcurve[s].tx-knotcurves[c].knotcurve[incp(s,-1,NP)].tx)/(knotcurves[c].knotcurve[s].length+knotcurves[c].knotcurve[incp(s,-1,NP)].length);
//                    double ny = 2.0*(knotcurves[c].knotcurve[s].ty-knotcurves[c].knotcurve[incp(s,-1,NP)].ty)/(knotcurves[c].knotcurve[s].length+knotcurves[c].knotcurve[incp(s,-1,NP)].length);
//                    double nz = 2.0*(knotcurves[c].knotcurve[s].tz-knotcurves[c].knotcurve[incp(s,-1,NP)].tz)/(knotcurves[c].knotcurve[s].length+knotcurves[c].knotcurve[incp(s,-1,NP)].length);
//                    double curvature = sqrt(nx*nx+ny*ny+nz*nz);
//                    nx /=curvature;
//                    ny /=curvature;
//                    nz /=curvature;
//                    double tx = knotcurves[c].knotcurve[s].tx ;
//                    double ty =  knotcurves[c].knotcurve[s].ty ;
//                    double tz = knotcurves[c].knotcurve[s].tz ;
//                    double bx = ty*nz - tz*ny;
//                    double by = tz*nx - tx*nz;
//                    double bz = tx*ny - ty*nx;
//                    knotcurves[c].knotcurve[s].nx = nx ;
//                    knotcurves[c].knotcurve[s].ny = ny ;
//                    knotcurves[c].knotcurve[s].nz = nz ;
//                    knotcurves[c].knotcurve[s].bx = bx ;
//                    knotcurves[c].knotcurve[s].by = by ;
//                    knotcurves[c].knotcurve[s].bz = bz ;
//                    knotcurves[c].knotcurve[s].curvature = curvature ;
//                }
//                // torsions with a central difference
//                for(s=0; s<NP; s++)
//                {
//                    double bx = knotcurves[c].knotcurve[s].bx;
//                    double by =  knotcurves[c].knotcurve[s].by;
//                    double bz = knotcurves[c].knotcurve[s].bz;
//
//                    double dnxds = 2.0*(knotcurves[c].knotcurve[incp(s,1,NP)].nx-knotcurves[c].knotcurve[incp(s,-1,NP)].nx)/(knotcurves[c].knotcurve[incp(s,1,NP)].length+knotcurves[c].knotcurve[incp(s,-1,NP)].length);
//                    double dnyds = 2.0*(knotcurves[c].knotcurve[incp(s,1,NP)].ny-knotcurves[c].knotcurve[incp(s,-1,NP)].ny)/(knotcurves[c].knotcurve[incp(s,1,NP)].length+knotcurves[c].knotcurve[incp(s,-1,NP)].length);
//                    double dnzds = 2.0*(knotcurves[c].knotcurve[incp(s,1,NP)].nz-knotcurves[c].knotcurve[incp(s,-1,NP)].nz)/(knotcurves[c].knotcurve[incp(s,1,NP)].length+knotcurves[c].knotcurve[incp(s,-1,NP)].length);
//
//                    double torsion = bx*dnxds+by*dnyds+bz*dnzds;
//                    knotcurves[c].knotcurve[s].torsion = torsion ;
//                }
//
//
//                // RIBBON TWIST AND WRITHE
//
//                for(s=0; s<NP; s++)
//                {
//
//                    // twist of this segment
//                    double ds = knotcurves[c].knotcurve[s].length;
//                    double dxds = knotcurves[c].knotcurve[s].tx;
//                    double dyds = knotcurves[c].knotcurve[s].ty;
//                    double dzds = knotcurves[c].knotcurve[s].tz;
//                    double bx = (knotcurves[c].knotcurve[incp(s,1,NP)].ax - knotcurves[c].knotcurve[s].ax)/ds;
//                    double by = (knotcurves[c].knotcurve[incp(s,1,NP)].ay - knotcurves[c].knotcurve[s].ay)/ds;
//                    double bz = (knotcurves[c].knotcurve[incp(s,1,NP)].az - knotcurves[c].knotcurve[s].az)/ds;
//                    knotcurves[c].knotcurve[s].twist = (dxds*(knotcurves[c].knotcurve[s].ay*bz - knotcurves[c].knotcurve[s].az*by) + dyds*(knotcurves[c].knotcurve[s].az*bx - knotcurves[c].knotcurve[s].ax*bz) + dzds*(knotcurves[c].knotcurve[s].ax*by - knotcurves[c].knotcurve[s].ay*bx))/(2*M_PI*sqrt(dxds*dxds + dyds*dyds + dzds*dzds));
//
//                    // "writhe" of this segment. writhe is nonlocal, this is the thing in the integrand over s
//                    knotcurves[c].knotcurve[s].writhe = 0;
//                    for(m=0; m<NP; m++)
//                    {
//                        if(s != m)
//                        {
//                            xdiff = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].xcoord + knotcurves[c].knotcurve[s].xcoord - knotcurves[c].knotcurve[incp(m,1,NP)].xcoord - knotcurves[c].knotcurve[m].xcoord);   //interpolate, consistent with fwd diff
//                            ydiff = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].ycoord + knotcurves[c].knotcurve[s].ycoord - knotcurves[c].knotcurve[incp(m,1,NP)].ycoord - knotcurves[c].knotcurve[m].ycoord);
//                            zdiff = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].zcoord + knotcurves[c].knotcurve[s].zcoord - knotcurves[c].knotcurve[incp(m,1,NP)].zcoord - knotcurves[c].knotcurve[m].zcoord);
//                            double dxdm = (knotcurves[c].knotcurve[incp(m,1,NP)].xcoord - knotcurves[c].knotcurve[m].xcoord)/(ds);
//                            double dydm = (knotcurves[c].knotcurve[incp(m,1,NP)].ycoord - knotcurves[c].knotcurve[m].ycoord)/(ds);
//                            double dzdm = (knotcurves[c].knotcurve[incp(m,1,NP)].zcoord - knotcurves[c].knotcurve[m].zcoord)/(ds);
//                            knotcurves[c].knotcurve[s].writhe += ds*(xdiff*(dyds*dzdm - dzds*dydm) + ydiff*(dzds*dxdm - dxds*dzdm) + zdiff*(dxds*dydm - dyds*dxdm))/(4*M_PI*(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)*sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff));
//                        }
//                    }
//
//                    //Add on writhe, twist
//                    knotcurves[c].writhe += knotcurves[c].knotcurve[s].writhe*ds;
//                    knotcurves[c].twist  += knotcurves[c].knotcurve[s].twist*ds;
//                    // while we are computing the global quantites, get the average position too
//                    knotcurves[c].xavgpos += knotcurves[c].knotcurve[s].xcoord/NP;
//                    knotcurves[c].yavgpos += knotcurves[c].knotcurve[s].ycoord/NP;
//                    knotcurves[c].zavgpos += knotcurves[c].knotcurve[s].zcoord/NP;
//                }
//
//
//                // the ghost grid has been useful for painlessly computing all the above quantities, without worrying about the periodic bc's
//                // but for storage and display, we should put it all in the box
//
//                // (1) construct the proper periodic co-ordinates from our ghost grid
//                double xupperlim = x(griddata.Nx -1 ,griddata);
//                double xlowerlim = x(0,griddata);
//                double deltax = griddata.Nx * h;
//                double yupperlim = y(griddata.Ny -1,griddata);
//                double ylowerlim = y(0,griddata);
//                double deltay = griddata.Ny * h;
//                double zupperlim = z(griddata.Nz -1 ,griddata);
//                double zlowermin = z(0,griddata);
//                double deltaz = griddata.Nz * h;
//                for(s=0; s<NP; s++)
//                {
//                    knotcurves[c].knotcurve[s].modxcoord = knotcurves[c].knotcurve[s].xcoord;
//                    knotcurves[c].knotcurve[s].modycoord = knotcurves[c].knotcurve[s].ycoord;
//                    knotcurves[c].knotcurve[s].modzcoord = knotcurves[c].knotcurve[s].zcoord;
//                    if(knotcurves[c].knotcurve[s].xcoord > xupperlim) {
//                        knotcurves[c].knotcurve[s].modxcoord = knotcurves[c].knotcurve[s].xcoord-deltax;
//                    } ;
//                    if(knotcurves[c].knotcurve[s].xcoord < xlowerlim) {
//                        knotcurves[c].knotcurve[s].modxcoord = knotcurves[c].knotcurve[s].xcoord+deltax;
//                    };
//                    if(knotcurves[c].knotcurve[s].ycoord > yupperlim) {
//                        knotcurves[c].knotcurve[s].modycoord = knotcurves[c].knotcurve[s].ycoord-deltay;
//                    };
//                    if(knotcurves[c].knotcurve[s].ycoord < ylowerlim) {
//                        knotcurves[c].knotcurve[s].modycoord = knotcurves[c].knotcurve[s].ycoord+deltay;
//                    };
//                    if(knotcurves[c].knotcurve[s].zcoord > zupperlim) {
//                        knotcurves[c].knotcurve[s].modzcoord = knotcurves[c].knotcurve[s].zcoord-deltaz;
//                    };
//                    if(knotcurves[c].knotcurve[s].zcoord < zlowermin)
//                    {
//                        knotcurves[c].knotcurve[s].modzcoord = knotcurves[c].knotcurve[s].zcoord+deltaz;
//                    };
//                }
//
//                // (2) standardise the knot such that the top right corner of the bounding box lies in the "actual" grid. This bounding box point may only lie
//                // off grid in the +ve x y z direction.
//                double xmax=knotcurves[c].knotcurve[0].xcoord;
//                double ymax=knotcurves[c].knotcurve[0].ycoord;
//                double zmax=knotcurves[c].knotcurve[0].zcoord;
//                for(s=0; s<NP; s++)
//                {
//                    if(knotcurves[c].knotcurve[s].xcoord>xmax)
//                    {
//                        xmax = knotcurves[c].knotcurve[s].xcoord;
//                    }
//                    if(knotcurves[c].knotcurve[s].ycoord>ymax)
//                    {
//                        ymax = knotcurves[c].knotcurve[s].ycoord;
//                    }
//                    if(knotcurves[c].knotcurve[s].zcoord>zmax)
//                    {
//                        zmax = knotcurves[c].knotcurve[s].zcoord;
//                    }
//                }
//
//                // get how many lattice shifts are needed
//                int xlatticeshift = (int) (round(xmax/(griddata.Nx *griddata.h)));
//                int ylatticeshift = (int) (round(ymax/(griddata.Ny *griddata.h)));
//                int zlatticeshift = (int) (round(zmax/(griddata.Nz *griddata.h)));
//                // perform the shift
//
//                for(int s=0; s<knotcurves[c].knotcurve.size(); s++)
//                {
//                    knotcurves[c].knotcurve[s].xcoord -= (double)(xlatticeshift) * (griddata.Nx *griddata.h);
//                    knotcurves[c].knotcurve[s].ycoord -= (double)(ylatticeshift) * (griddata.Ny *griddata.h);
//                    knotcurves[c].knotcurve[s].zcoord -= (double)(zlatticeshift) * (griddata.Nz *griddata.h);
//                }
//                // now we've done these shifts, we'd better move the knotcurve average position too.
//                knotcurves[c].xavgpos = 0;
//                knotcurves[c].yavgpos = 0;
//                knotcurves[c].zavgpos = 0;
//                for(int s=0; s<knotcurves[c].knotcurve.size(); s++)
//                {
//                    knotcurves[c].xavgpos += knotcurves[c].knotcurve[s].xcoord/NP;
//                    knotcurves[c].yavgpos += knotcurves[c].knotcurve[s].ycoord/NP;
//                    knotcurves[c].zavgpos += knotcurves[c].knotcurve[s].zcoord/NP;
//                }
//                c++;
//            }
//            // if we did hit a boundary, just strike the curve we made from the record. It lives on in the marked array!
//            else
//            {
//                knotcurves.pop_back();
//            }
//        }
//    }
//    // the order of the components within the knotcurves vector is not guaranteed to remain fixed from timestep to timestep. thus, componenet 0 at one timtestep could be
//    // components 1 at the next. the code needs a way of tracking which componenet is which.
//    // at the moment, im doing this by fuzzily comparing summary stats on the components - at this point, the length twist and writhe.
//
//    // these variables have the summary stats from the last timestep
//
//    static vector<double> oldwrithe(knotcurves.size());
//    static vector<double> oldtwist(knotcurves.size());
//    static vector<double> oldlength(knotcurves.size());
//    static vector<double> oldxavgpos(knotcurves.size());
//    static vector<double> oldyavgpos(knotcurves.size());
//    static vector<double> oldzavgpos(knotcurves.size());
//    static bool first = true;
//    static vector<int> permutation(knotcurves.size());
//
//    if(first)
//    {
//        for(int i = 0; i<knotcurves.size();i++)
//        {
//            oldwrithe[i] = knotcurves[i].writhe;
//            oldlength[i] = knotcurves[i].length;
//            oldtwist[i] = knotcurves[i].twist;
//            oldxavgpos[i] = knotcurves[i].xavgpos;
//            oldyavgpos[i] = knotcurves[i].yavgpos;
//            oldzavgpos[i] = knotcurves[i].zavgpos;
//            permutation[i] = i;
//        }
//
//    }
//    else
//    {
//        for(int i = 0; i<knotcurves.size();i++)
//        {
//            double minscore = INFINITY;
//            for(int j = 0; j<knotcurves.size();j++)
//            {
//                double score = fabs((knotcurves[j].length - oldlength[i])/(Nx*h))+fabs(knotcurves[j].writhe - oldwrithe[i]);
//                score += fmod(fabs((knotcurves[j].xavgpos - oldxavgpos[i])),Nx*h)/(Nx*h);
//                score += fmod(fabs((knotcurves[j].yavgpos - oldyavgpos[i])),Ny*h)/(Ny*h);
//                score += fmod(fabs((knotcurves[j].zavgpos - oldzavgpos[i])),Nz*h)/(Nz*h);
//                if(score<minscore) {permutation[i] = j; minscore = score;}
//
//            }
//        }
//
//        // if the permutation isn't valid, print a warning and reset it - it's better to have the curves swap then have them incorrectly both mapped to the same curve.
//        bool isvalidperm = true;
//        for(int i=0;i<knotcurves.size();i++)
//        {
//            std::vector<int>::iterator it = std::find(permutation.begin(),permutation.end(),i);
//            if(it == permutation.end()){isvalidperm=false;}
//        }
//
//        if(!isvalidperm)
//        {
//            cout << "the permutation was not valid. resetting! \n";
//            for(int i = 0; i<knotcurves.size();i++)
//            {
//                permutation[i] = i;
//            }
//
//        }
//
//
//        // apply the permutation to the list of lengths etc. It is now "correct" in the sense that index [0] really is component 0 etc. these labellings are arbitrarlly set at the simulations start and
//        // must be consistently carried forward
//        for(int i = 0; i<knotcurves.size();i++)
//        {
//            oldwrithe[i] = knotcurves[permutation[i]].writhe;
//            oldlength[i] = knotcurves[permutation[i]].length;
//            oldtwist[i] = knotcurves[permutation[i]].twist;
//            oldxavgpos[i] = knotcurves[permutation[i]].xavgpos;
//            oldyavgpos[i] = knotcurves[permutation[i]].yavgpos;
//            oldzavgpos[i] = knotcurves[permutation[i]].zavgpos;
//
//        }
//        // i can't be bothered to do this the smart way.
//        vector<knotcurve> tempknotcurves(knotcurves.size());
//        for(int i = 0; i<knotcurves.size();i++)
//        {
//            tempknotcurves[i]=knotcurves[permutation[i]];
//        }
//        knotcurves=tempknotcurves;
//    }
//    first = false;
//}
//
// some little functions which go back and forth from point to index
int pt(const int k,const  int l,const  int m)       //convert i,j,k to single index
{
    return (m*Ly*Lx+l*Lx+k);
}

double x(int k)
{
    return (k+0.5-Lx/2.0);
}
double y(int l)
{
    return (l+0.5-Ly/2.0);
}
double z(int m)
{
    return (m+0.5-Lz/2.0);
}
