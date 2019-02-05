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

#pragma omp parallel default(none) shared(nx,ny,nz,px,py,pz,hx,hy,hz,hpx,hpy,hpz,n,step,cout)
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
                    writeVTKfiles(n,nx,ny,nz,px,py,pz);       // output VTK files for use with ParaView
                    writeBENDfiles(n, nx,ny,nz);       // output VTK files for use with ParaView
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
                if (rho < RR)
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
                    nx[j] = -sin(M_PI*rho/RR)*sin(phi+theta); // could also use -theta; they seem to be equivalent in energy
                    ny[j] = sin(M_PI*rho/RR)*cos(phi+theta); // could also use -theta
                    nz[j] = -cos(M_PI*rho/RR);

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

#if BC // cheap fix for the Dirichlet boundary conditions
                if (m==0) {xup=j; xdwn=j; yup=j; ydwn=j; zup=j; zdwn=j;}
                if (m==Lz-1) {xup=j; xdwn=j; yup=j; ydwn=j; zup=j; zdwn=j;}
#endif
                //                                // calculate first order derivatives
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

// some little functions which go back and forth from point to index
int pt(const int k,const  int l,const  int m)       //convert i,j,k to single index
{
    return (m*Ly*Lx+l*Lx+k);
}
