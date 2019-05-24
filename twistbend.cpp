/**************************************************************/
/*      Finite difference code for Twist Bend Nematics        */
/*         created November 2018, Gareth Alexander            */
/**************************************************************/

#include "twistbend.h"
#include "SolidAngle/SolidAngle.h"
#include "SolidAngle/InputOutput.h"
#include "SolidAngle/Geometry.h"
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
//includes for the signal processing
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

using namespace std;

int main(int argc, char** argv) 
{

    double* nx=new double[LL];
    double* ny=new double[LL];
    double* nz=new double[LL];

    double* bx=new double[LL];
    double* by=new double[LL];
    double* bz=new double[LL];

    double* bmag=new double[LL];
    double* twist=new double[LL];
    double* FreeEnergy=new double[LL];

    double* tx=new double[LL];
    double* ty=new double[LL];
    double* tz=new double[LL];

    double* hx=new double[LL];
    double* hy=new double[LL];
    double* hz=new double[LL];

    double* px=new double[LL];
    double* py=new double[LL];
    double* pz=new double[LL];

    double* hpx=new double[LL];
    double* hpy=new double[LL];
    double* hpz=new double[LL];

    bool* mask=new bool[LL];

    int n =0;
    startconfig(n, nx,ny,nz,px,py,pz);
    setupmask(mask);
    cout << " Finished the initialisation, starting simulation" << endl;

#pragma omp parallel default(none) shared(nx,ny,nz, hx,hy,hz,px,py,pz,hpx,hpy,hpz, bx,by,bz,bmag, tx,ty,tz,twist, FreeEnergy, n,mask,cout)
    {
        while(n<=Nmax)
        {
#pragma omp single
            {
                if (n%1000==0)
                {
                    cout << "timestep " << n << endl;
                }

                if (n%vtkstepskip==0)
                {
                    cout << "writing VTK files at timestep " << n << endl;
                    computeBendTwistEnergyOrientation(n,nx,ny,nz,bx,by,bz,px,py,pz,bmag,twist,FreeEnergy,tx,ty,tz);
                    writeVTKfiles(n,nx,ny,nz,px,py,pz,bx,by,bz,tx,ty,tz,twist,FreeEnergy);       // output VTK files for use with ParaView
                    // writeCustomPoints(n,nx,ny,nz,bx,by,bz);       // output VTK files for use with ParaView
                }
                if ((n>=curvestarttime) &&(n%curvestepskip==0))
                {
                    cout << "writing the bend zeros at timestep " << n << endl;
                    computeBendTwistEnergyOrientation(n,nx,ny,nz,bx,by,bz,px,py,pz,bmag,twist,FreeEnergy,tx,ty,tz);
                    Link Curve;
                    FindBendZeros(Curve,nx,ny,nz,bx,by,bz, bmag,tx,ty,tz,mask);
                    print_Curve(n,Curve, "vtk_bendzeros");

                }
                n++;
            }
            update(nx,ny,nz,px, py,pz, hx, hy,hz,hpx, hpy,hpz);
        }
    }


    delete nx;
    delete ny;
    delete nz;

    delete px;
    delete py;
    delete pz;


    delete bx;
    delete by;
    delete bz;

    delete bmag;

    delete tx;
    delete ty;
    delete tz;

    delete hx;
    delete hy;
    delete hz;

    delete hpx;
    delete hpy;
    delete hpz;
} // end main

void setupmask(bool* mask)
{
    // currently set up to just mask off the edge of the box
    int buffer =4;
    for(int k=0;k<Lx;k++)
    {
        for(int l=0; l<Ly; l++)
        {
            for(int m=0; m<Lz; m++)
            {
                mask[pt(k,l,m)]=false;
                bool xoverflow = (k < buffer) || (k > Lx - buffer);
                bool yoverflow = (l < buffer) || (l > Ly - buffer);
                bool zoverflow = (m < buffer) || ( m> Lz - buffer);
                if(xoverflow || yoverflow || zoverflow){ mask[pt(k,l,m)] =true;}
            }
        }
    }
}
void startconfig(int & n, double* nx, double* ny,double* nz,double* px, double* py,double* pz)
{
    switch(InitialisationMethod)
    {
        case FROM_FUNCTION:
            {
                // k,l,m  are x,y,z indices, j is 1-d loop index.
                int j,k,l,m;
                int HELICONICAL,HOPF,SQUARE,HOPF2, BENDZERO;
                double k0,l0,m0;

                // define texture to simulate -- 0 or 1
                HELICONICAL = 0; // standard heliconical texture
                HOPF = 0; // Hopf texture
                HOPF2 = 0; // Hopf texture
                SQUARE = 0;
                BENDZERO=1;

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
                    // the raidus of the central core
                    double RR = 0.2*Lx;  // scale
                    // how fat the hopfion should be
                    double RR2 =0.05*Lx;  // scale
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
                            //nx[j] = -sin(M_PI*rho/RR2)*sin(phi-theta); // could also use -theta; they seem to be equivalent in energy
                            //ny[j] = sin(M_PI*rho/RR2)*cos(phi-theta); // could also use -theta
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
                if( BENDZERO == 1)
                {
                    int xup,xdwn,yup,ydwn,zup,zdwn;
                    double rho,theta,phi, rr, norm;    // variables for hopf texture
                    // the radius of the circle 
                    double RR = 0.25*Lx;  // scale
                    // how fat the tubular neighborhood should be
                    double RR2 =0.15*Lx;  // scale

                    // solid angle fn for the loop
                    Link Curve;
                    Curve.Components.resize(1);
                    double maxxin = 0;
                    double maxyin = 0;
                    double maxzin = 0;
                    double minxin = 0;
                    double minyin = 0;
                    double minzin = 0;
                    for(int i=0; i<100;i++)
                    {
                        double xcoord,ycoord,zcoord;
                        knotpoint Point;
                        xcoord = RR*cos((2*M_PI*i)/100.0)+k0;
                        ycoord = RR*sin((2*M_PI*i)/100.0)+l0;
                        zcoord =m0;
                        Point.xcoord = xcoord;
                        Point.ycoord = ycoord;
                        Point.zcoord = zcoord;
                        Curve.Components[0].knotcurve.push_back(Point);
                        if(xcoord>maxxin) maxxin = xcoord;
                        if(ycoord>maxyin) maxyin = ycoord;
                        if(zcoord>maxzin) maxzin = zcoord;
                        if(xcoord<minxin) minxin = xcoord;
                        if(ycoord<minyin) minyin = ycoord;
                        if(zcoord<minzin) minzin = zcoord;
                    }
                    cout << "Filling in the Geometry of the Input Curve \n";
                    ComputeGeometry(Curve);
                    double* omega=new double[LL];
                    OutputScaledKnot(Curve);
                    ComputeSolidAngleAllPoints(Curve,omega);
                    OutputSolidAngle(Curve,omega,"solidangle");

                    for (j=0; j<LL; j++)
                    {
                        // begin with a heliconical director
                        // director field
                        phi = omega[j];
                        
                        nx[j] = (1.0/sqrt(2))*sin(phi);
                        ny[j] = (1.0/sqrt(2))*sin(phi);
                        nz[j]= cos(phi);
                       
                        // polarisation
                        px[j] = -sin(qh*m);
                        py[j] = cos(qh*m);
                        pz[j] = 0.0;

                        // now we do our tubular nbhd

                        // xy plane distance to z axis.
                        rr = sqrt((k-k0)*(k-k0)+(l-l0)*(l-l0));
                        rho = sqrt((rr-RR)*(rr-RR)+(m-m0)*(m-m0));

                        // inside the tube, paste our local model for a bend zero
                        if (rho < RR2)
                        {
                            double phi = (M_PI/2.0) - 0.05;
                            theta = atan2(l-l0,k-k0);

                            // the tubular neighborhood framing system
                            double e30x = -sin(theta);
                            double e30y = cos(theta);
                            double e30z = 0;
                            norm = sqrt(e30x*e30x + e30y*e30y +e30z*e30z);
                            e30x /=norm;
                            e30y /=norm;
                            e30z /=norm;

                            double e10x =cos(theta); 
                            double e10y =sin(theta);
                            double e10z = 0;
                            norm = sqrt(e10x*e10x + e10y*e10y +e10z*e10z);
                            e10x /=norm;
                            e10y /=norm;
                            e10z /=norm;

                            double e20x = 0;
                            double e20y = 0;
                            double e20z = 1;

                            // the director framing
                            double d30x = cos(phi)*e30x + sin(phi)*e20x;
                            double d30y = cos(phi)*e30y + sin(phi)*e20y;
                            double d30z = cos(phi)*e30z + sin(phi)*e20z;

                            double d10x = e10x;
                            double d10y = e10y;
                            double d10z = e10z;

                            double d20x = d30y*d10z - d30z*d10y;
                            double d20y = d30z*d10x - d30x*d10x;
                            double d20z = d30x*d10y - d30y*d10z;

                            // okay, now get the uv coodinates of our xyz point
                            // r- r0
                            double drx = ((k-k0)-(RR/rr)*(k-k0));
                            double dry = ((l-l0)-(RR/rr)*(l-l0));
                            double drz = (m-m0);

                            double u = drx*e10x +dry*e10y + drz*e10z;
                            double v = drx*e20x +dry*e20y + drz*e20z;

                            // these coordinates are in the same dimensions as the tube radius RR2.
                            // we want to rescale the uv coordinates across this tube, say so that u, v reach maximum values of 1/10 on the tube
                            u = u/(RR2);
                            v = v/(RR2);

                            // okay now give the director in the tubular nbhd

                           // double nxcomp=(u + v*v +3*u*u); 
                           // double nycomp =u -2*v*v + 2*u*u;
                           // double nxcomp=u*(1-sqrt(u*u +v*v));
                           // double nycomp=2*u*(1-sqrt(u*u +v*v));
                           // double nzcomp = sqrt(1-nx[j]*nx[j]-ny[j]*ny[j]);
                           nx[j]*=rho/RR2;
                           ny[j]*=rho/RR2;

                           nx[j] = nxcomp*d10x + nycomp*d20x +nzcomp*d30x;
                            ny[j] = nxcomp*d10y + nycomp*d20y +nzcomp*d30y;
                           nz[j]= sqrt(1-nx[j]*nx[j]-ny[j]*ny[j]);

                            px[j] = 0;
                            py[j] = 0;
                            pz[j] = 0.0;

                            // going to try and smooth into the uniform state in a radial fashion
                            
                            /*
                            if(rho>0.01*RR2)
                            {
                                nx[j] = exp(-0.1*(rho/RR2))*nx[j];
                                ny[j] = exp(-0.1*(rho/RR2))*ny[j];
                                nz[j]= sqrt(1-nx[j]*nx[j]-ny[j]*ny[j]);
                            }
                            */
                            
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
                }
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
        case FROM_SOLIDANGLE:
            {
                cout << "Initialising from a twisted solid angle function ...\n";
                InitialiseSystemParameters();
                int j,k,l,m;
                double omega,alpha,R,rho;
                // for initialising the polarisation
                int xup,xdwn,yup,ydwn,zup,zdwn;
                double norm;
                // value of the Hopf invariant
                int h;
                h = 2;
                // initialise the knot
                Link Curve;
                if(0==InitialiseFromFile(Curve))
                {
                    cout << "Filling in the Geometry of the Input Curve \n";
                    ComputeGeometry(Curve);
                    OutputScaledKnot(Curve);
                }
                // initial configuration
                k=l=m=0;
                // radius for the tubular neighbourhood of the curve
                R = 0.28*Lz;

                viewpoint Point;
                omega=0.0; alpha=0.0; // overly cautious!!

                // compute omega and the initial director field
                cout << "starting twisted solid angle calculation ...\n";
                for (j=0; j<LL; j++) {
                    // the far field director is the standard heliconical texture -- right handed
                    nx[j] = sin(thetah)*cos(qh*m);
                    ny[j] = sin(thetah)*sin(qh*m);
                    nz[j] = cos(thetah);
                    // polarisation
                    px[j] = -sin(qh*m);
                    py[j] = cos(qh*m);
                    pz[j] = 0.0;
                    // left handed
                    //      nx[j] *= -1.0;
                    //      px[j] *= -1.0;

                    Point.xcoord = x(k); //1.0*k-Lx/2.0+0.5;
                    Point.ycoord = y(l); //1.0*l-Ly/2.0+0.5;
                    Point.zcoord = z(m); //1.0*m-Lz/2.0+0.5;

                    // we need to calculate the distance to the curve
                    rho = ComputeDistanceOnePoint(Curve,Point);

                    // put a Skyrmion texture inside the tubular neighbourhood
                    if (rho < R) {
                        // here is the solid angle
                        omega = ComputeSolidAngleOnePoint(Curve,Point);
                        // simple thing -- not adapted to distinct components yet!!
                        alpha = ComputeLongitudinalPhase(Curve,Point);

                        // set the director -- fix me!
                        // still need to be able to add Dehn twists to each component separately
                        nx[j] = sin(M_PI*rho/R)*cos(0.5*omega-h*alpha);
                        ny[j] = sin(M_PI*rho/R)*sin(0.5*omega-h*alpha);
                        nz[j] = -cos(M_PI*rho/R);
                    }

                    k++;
                    if (k==Lx) {l++; k=0;}
                    if (l==Ly) {m++; l=0;}
                }
                // set the polarisation to the bend of the director
                k=l=m=0;
                for (j=0; j<LL; j++)
                {
                    Point.xcoord = k; //1.0*k-Lx/2.0+0.5;
                    Point.ycoord = l; //1.0*l-Ly/2.0+0.5;
                    Point.zcoord = m; //1.0*m-Lz/2.0+0.5;

                    // we need to calculate the distance to the curve
                    rho = ComputeDistanceOnePoint(Curve,Point);

                    if (rho < R) {
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
                        norm += 0.000001; // perhaps overly cautious!!
                        px[j] /= norm; py[j] /= norm; pz[j] /= norm;
                    }
                    // deal with the periodic boundaries
                    k++;
                    if (k==Lx) {l++; k=0;}
                    if (l==Ly) {m++; l=0;}
                }

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

void computeBendTwistEnergyOrientation(const int n, const double* nx,const double* ny,const double* nz, double* bx, double* by, double* bz,const double* px, const double* py, const double* pz, double* bmag, double* twist,double* FreeEnergy, double* tx, double* ty, double* tz)
{

    int j;
    // arrays storing the circulaiton, which we will take the curl of
    double* cx=new double[LL];
    double* cy=new double[LL];
    double* cz=new double[LL];

    int k,l,m;
    int xup,xdwn,yup,ydwn,zup,zdwn;
    double Dxnx,Dxny,Dxnz,Dynx,Dyny,Dynz,Dznx,Dzny,Dznz;
    double Dxpx,Dypx,Dzpx,Dxpy,Dypy,Dzpy, Dxpz,Dypz,Dzpz;
    double Dxbx,Dxby,Dxbz,Dybx,Dyby,Dybz,Dzbx,Dzby,Dzbz;
    double Dxcy,Dxcz,Dycx,Dycz,Dzcx,Dzcy;

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

        Dxpx = (px[xup]-px[xdwn])/2.0;
        Dypx = (px[yup]-px[ydwn])/2.0;
        Dzpx = (px[zup]-px[zdwn])/2.0;

        Dxpy = (py[xup]-py[xdwn])/2.0;
        Dypy = (py[yup]-py[ydwn])/2.0;
        Dzpy = (py[zup]-py[zdwn])/2.0;

        Dxpz = (pz[xup]-pz[xdwn])/2.0;
        Dypz = (pz[yup]-pz[ydwn])/2.0;
        Dzpz = (pz[zup]-pz[zdwn])/2.0;


#if BC // one-sided derivates at boundaries
        if (m==0) {
            Dznx = (-3.0*nx[j]+4.0*nx[zup]-nx[zdwn])/2.0;
            Dzny = (-3.0*ny[j]+4.0*ny[zup]-ny[zdwn])/2.0;
            Dznz = (-3.0*nz[j]+4.0*nz[zup]-nz[zdwn])/2.0;

            Dzpx = (-3.0*px[j]+4.0*px[zup]-px[zdwn])/2.0;
            Dzpy = (-3.0*py[j]+4.0*py[zup]-py[zdwn])/2.0;
            Dzpz = (-3.0*pz[j]+4.0*pz[zup]-pz[zdwn])/2.0;
        }
        if (m==Lz-1) {
            Dznx = (3.0*nx[j]-4.0*nx[zdwn]+nx[zup])/2.0;
            Dzny = (3.0*ny[j]-4.0*ny[zdwn]+ny[zup])/2.0;
            Dznz = (3.0*nz[j]-4.0*nz[zdwn]+nz[zup])/2.0;

            Dzpx = (3.0*px[j]-4.0*px[zdwn]+px[zup])/2.0;
            Dzpy = (3.0*py[j]-4.0*py[zdwn]+py[zup])/2.0;
            Dzpz = (3.0*pz[j]-4.0*pz[zdwn]+pz[zup])/2.0;
        }
#endif

        // calculate bend
        bx[j] = nx[j]*Dxnx + ny[j]*Dynx + nz[j]*Dznx;
        by[j] = nx[j]*Dxny + ny[j]*Dyny + nz[j]*Dzny;
        bz[j] = nx[j]*Dxnz + ny[j]*Dynz + nz[j]*Dznz;
        bmag[j] = sqrt(bx[j]*bx[j]+by[j]*by[j]+bz[j]*bz[j]);

        // calculate twist
        twist[j] = nx[j]*(Dynz-Dzny) + ny[j]*(Dznx-Dxnz) + nz[j]*(Dxny-Dynx);

        // calculate free energy
        double gradnsquared =(Dxnx*Dxnx + Dxny*Dxny+ Dxnz*Dxnz) + (Dynx*Dynx + Dyny*Dyny+ Dynz*Dynz)+ (Dznx*Dznx + Dzny*Dzny+ Dznz*Dznz);
        double gradpsquared=(Dxpx*Dxpx + Dxpy*Dxpy+ Dxpz*Dxpz) +( Dypx*Dypx + Dypy*Dypy+ Dypz*Dypz) +( Dzpx*Dzpx + Dzpy*Dzpy+ Dzpz*Dzpz);
        double bdotp = bx[j]*px[j] + by[j]*py[j] + bz[j]*pz[j];
        double pdotp = px[j]*px[j] + py[j]*py[j] + pz[j]*pz[j];

        FreeEnergy[j] = (K/2)*(gradnsquared) - lambda*(bdotp) + (C/2)*(gradpsquared) + (U/4)*(1-pdotp)*(1-pdotp);

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
        Dycx = (cx[yup]-cx[ydwn])/2.0;
        Dzcx = (cx[zup]-cx[zdwn])/2.0;

        Dxcy = (cy[xup]-cy[xdwn])/2.0;
        Dzcy = (cy[zup]-cy[zdwn])/2.0;

        Dxcz = (cz[xup]-cz[xdwn])/2.0;
        Dycz = (cz[yup]-cz[ydwn])/2.0;

#if BC // one-sided derivates at boundaries
        if (m==0) {
            Dzcx = (-3.0*cx[j]+4.0*cx[zup]-cx[zdwn])/2.0;
            Dzcy = (-3.0*cy[j]+4.0*cy[zup]-cy[zdwn])/2.0;
        }
        if (m==Lz-1) {
            Dzcx = (3.0*cx[j]-4.0*cx[zdwn]+cx[zup])/2.0;
            Dzcy = (3.0*cy[j]-4.0*cy[zdwn]+cy[zup])/2.0;
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
    delete cx;
    delete cy;
    delete cz;
}

void FindBendZeros(Link& Curve,double* nx,double* ny,double* nz, double* bx,double* by,double* bz, double *magb,double* tx,double* ty,double* tz, bool* mask)
{
    /*
     *
     * SETTINGS
     *
     */
    // threshold for detecting another bend zero
    const double threshold =0.016;
    // The size of the hemisphere in front of the current point which the code searches for the next point to go to. The data is on a grid of size 1, so it makes no sense for this to be way below 1. some O(1) number here is good.
    const double sphereradius =1.2;
    // How many steps on the hemisphere the GSL minimizer should make. Ive found setting this too high makes the code more jagged. I don't reaaaaallly know why, perhaps being really picky about where the minimum is in "noisy" data causes it to run off.
    const int numgsliterations =3;
    // the initial stepsize for the GSL minimzer. Something below 1, I use something well below it.
    const double initialgslstepsize = 0.05;
    // how close shoud our current point be to the initial point before we terminate and close the loop? Again, something O(1), I don't have an amazing feel.
    const double terminationdistance = 3;
    // when we add the first point, we will be close to the start point and we will hit the termination condition above. To avoid this I just require us to be some number of points along the tracing before termination is allowed. This is the number of points
    const int minnumpoints = 30;
    // specify a lengthscale to kill fluctuations shorter than, for the curve smoothing. Some O(1- 10) numer
    const double filterlength=3;
    // when we do the two push offs, need to specify how far to push. roughly O(1) numbers again.
    double firstpushdist = 1;
    /*
     *
     *
     *
     */

    // initialise the tricubic interpolator for ucvmag
    likely::TriCubicInterpolator interpolatedmagb(magb, 1, Lx,Ly,Lz);

    // an array holding points we have already visited, so we dont keep tracing out the same bend zero
    vector<int> marked(Lx*Ly*Lz,0);


    // GSL initialization
    const gsl_multimin_fminimizer_type *Type;
    Type = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *minimizerstate;
    minimizerstate = gsl_multimin_fminimizer_alloc (Type,2);

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
                    if( magb[n] < bmin && marked[n]==0 && mask[n]==0)
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
            Curve.Components.push_back(knotcurve() );
            Curve.Components[c].knotcurve.push_back(knotpoint());
            Curve.Components[c].knotcurve[0].xcoord=kmin;
            Curve.Components[c].knotcurve[0].ycoord=lmin;
            Curve.Components[c].knotcurve[0].zcoord=mmin;

            int s=1;
            bool lastpoint=false;
            bool finish=false;
            bool BoundaryTerminationCondition=false;
            double flowsign=1.0;
            int bufferhits=0;
            // initially assume the curve we are tracing will close
            Curve.Components[c].closed =true;
            while (!finish)
            {

                // the tangent vector cirlcirculation
                double txs=0;
                double tys=0;
                double tzs=0;

                /**Find nearest gridpoint**/
                int idwn = (int) ((Curve.Components[c].knotcurve[s-1].xcoord));
                int jdwn = (int) ((Curve.Components[c].knotcurve[s-1].ycoord));
                int kdwn = (int) ((Curve.Components[c].knotcurve[s-1].zcoord));
                /*curve to gridpoint down distance*/
                double xd = Curve.Components[c].knotcurve[s-1].xcoord - idwn;
                double yd = Curve.Components[c].knotcurve[s-1].ycoord - jdwn;
                double zd = Curve.Components[c].knotcurve[s-1].zcoord - kdwn;
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
                    double projy = (cartesian[1] - (cartesian[0]*txs + cartesian[1]*tys + cartesian[2]*tzs)*tys);
                    double projz = (cartesian[2] - (cartesian[0]*txs + cartesian[1]*tys + cartesian[2]*tzs)*tzs);
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

                // put them in the Curve
                Curve.Components[c].knotcurve[s-1].tx = txs;
                Curve.Components[c].knotcurve[s-1].ty = tys;
                Curve.Components[c].knotcurve[s-1].tz = tzs;

                Curve.Components[c].knotcurve[s-1].e1x =e1x;
                Curve.Components[c].knotcurve[s-1].e1y = e1y;
                Curve.Components[c].knotcurve[s-1].e1z = e1z;

                Curve.Components[c].knotcurve[s-1].e2x = e2x;
                Curve.Components[c].knotcurve[s-1].e2y = e2y;
                Curve.Components[c].knotcurve[s-1].e2z = e2z;

                if(!lastpoint)
                {
                    //pop a new point on the stack
                    Curve.Components[c].knotcurve.push_back(knotpoint());
                    // okay we have our directions in the plane we want to perfrom the line minimisation in. Time to do it
                    // the point
                    gsl_vector* mypt = gsl_vector_alloc (3);
                    gsl_vector_set (mypt, 0, Curve.Components[c].knotcurve[s-1].xcoord);
                    gsl_vector_set (mypt, 1, Curve.Components[c].knotcurve[s-1].ycoord);
                    gsl_vector_set (mypt, 2, Curve.Components[c].knotcurve[s-1].zcoord);
                    // the tangent
                    // NOTE THIS IS WEIGHTED BY THE FLOWSIGN. Says whether we are flying forwards or backwards along the tangent vector field
                    gsl_vector* t = gsl_vector_alloc (3);
                    gsl_vector_set (t, 0, flowsign*txs);
                    gsl_vector_set (t, 1, flowsign*tys);
                    gsl_vector_set (t, 2, flowsign*tzs);
                    // one vector in the plane we with to minimize in
                    gsl_vector* e1 = gsl_vector_alloc (3);
                    gsl_vector_set (e1, 0, e1x);
                    gsl_vector_set (e1, 1, e1y);
                    gsl_vector_set (e1, 2, e1z);
                    // and the other
                    gsl_vector* e2 = gsl_vector_alloc (3);
                    gsl_vector_set (e2, 0, e2x);
                    gsl_vector_set (e2, 1, e2y);
                    gsl_vector_set (e2, 2, e2z);



                    struct parameters params; struct parameters* pparams = &params;
                    pparams->sphereradius = sphereradius;
                    pparams->ucvmag=&interpolatedmagb;
                    pparams->mypt = mypt; pparams->t = t;pparams->e1=e1; pparams->e2=e2;
                    // settings
                    gsl_multimin_function F;
                    F.n=2;
                    F.f = &my_minimisation_function;
                    F.params = (void*) pparams;
                    // stepsize and origin
                    gsl_vector* stepsize = gsl_vector_alloc (2);
                    gsl_vector_set (stepsize, 0, initialgslstepsize);
                    gsl_vector_set (stepsize, 1, initialgslstepsize);

                    gsl_vector* minimum = gsl_vector_alloc (2);
                    gsl_vector_set (minimum, 0, 0);
                    gsl_vector_set (minimum, 1, 0);

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
                    while (status == GSL_CONTINUE && iter <numgsliterations);

                    double x =gsl_vector_get(minimizerstate->x, 0);
                    double y =gsl_vector_get(minimizerstate->x, 1);
                    gsl_vector_scale(e1,x);
                    gsl_vector_scale(e2,y);
                    gsl_vector_scale(t,sqrt(1-x*x-y*y));
                    gsl_vector_add(mypt,e1);
                    gsl_vector_add(mypt,e2);
                    gsl_vector_add(mypt,t);
                    Curve.Components[c].knotcurve[s].xcoord = gsl_vector_get(mypt, 0);
                    Curve.Components[c].knotcurve[s].ycoord= gsl_vector_get(mypt, 1);
                    Curve.Components[c].knotcurve[s].zcoord= gsl_vector_get(mypt, 2);
                    gsl_vector_free(mypt);
                    gsl_vector_free(e1);
                    gsl_vector_free(e2);
                    gsl_vector_free(t);
                    gsl_vector_free(stepsize);

                    // mark the marked array with this curve
                    marked[pt((int) Curve.Components[c].knotcurve[s].xcoord,(int) Curve.Components[c].knotcurve[s].ycoord,(int) Curve.Components[c].knotcurve[s].zcoord)]=-1;

                    // did we just go into the buffer region around the edge of the box?
                    if(KnotpointInMask(Curve.Components[c].knotcurve[s],mask ))
                    {
                        bufferhits++;
                        if(bufferhits==1)
                        {
                            // if we did, start again from this point on the grid edge, but walking backwards!
                            // strip off everything but the SECOND last point in the vector ,i.e the last point not in the buffer
                            // grab that second last point
                            knotpoint newstartingpoint = Curve.Components[c].knotcurve[s];
                            // erase
                            Curve.Components[c].knotcurve.clear();
                            // put the start point back on
                            Curve.Components[c].knotcurve.push_back(newstartingpoint);
                            // reset s, set the flow sign to reverse
                            s =0;
                            flowsign = -1;
                        }
                        else if(bufferhits==2)
                        {
                            BoundaryTerminationCondition=true;
                            Curve.Components[c].knotcurve.pop_back();
                            Curve.Components[c].closed = false;
                        }
                    }
                }
                // okay we got our new point. Now the question is, when should we terminate this process? Do it by whether we are close to the start point after some nontrivial number of steps
                double xdiff = Curve.Components[c].knotcurve[0].xcoord - Curve.Components[c].knotcurve[s].xcoord;     //distance from start/end point
                double ydiff = Curve.Components[c].knotcurve[0].ycoord - Curve.Components[c].knotcurve[s].ycoord;
                double zdiff = Curve.Components[c].knotcurve[0].zcoord - Curve.Components[c].knotcurve[s].zcoord;
                // when we terminate, get the stats for the last point but dont add another on
                if(lastpoint)finish=true;
                if( (sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff) < terminationdistance && s > minnumpoints )|| BoundaryTerminationCondition || s>2000)
                {
                    lastpoint = true;
                }
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
        c++;
    }
    // clear up
    gsl_multimin_fminimizer_free(minimizerstate);

    // okay, we have the bend zeros. Now do some geometry
    // curve lengths
    for(int c=0; c<Curve.Components.size(); c++)
    {
        int NP = Curve.Components[c].knotcurve.size();
        for(int s=0; s<NP; s++)
        {
            double dx = (Curve.Components[c].knotcurve[(s+1)%NP].xcoord - Curve.Components[c].knotcurve[s].xcoord);
            double dy = (Curve.Components[c].knotcurve[(s+1)%NP].ycoord - Curve.Components[c].knotcurve[s].ycoord);
            double dz = (Curve.Components[c].knotcurve[(s+1)%NP].zcoord - Curve.Components[c].knotcurve[s].zcoord);
            double deltas = sqrt(dx*dx+dy*dy+dz*dz);
            Curve.Components[c].knotcurve[s].length = deltas;
        }
    }

    //get the solid angle framings of each individual curve component. These have zero linking number with the components itself.
    // Note: this isnt the same as the solid angle framing for the whole link. in that case, that framing has SL(k_i) = -sum Lk(k_i,k_j). Im treating each component individually, ignoring the others.
    double* phi = new double[LL];
    for(int c=0; c<Curve.Components.size(); c++)
    {
        Link tempCurve;
        tempCurve.Components.push_back(Curve.Components[c]);
        ComputeSolidAngleAllPoints(tempCurve, phi);
        ComputeSolidAngleFraming(tempCurve,phi);
        Curve.Components[c]=tempCurve.Components[0];
    }
    delete phi;

    // get the bend pushoff
    //  interpolation of the b vector
    likely::TriCubicInterpolator interpolatedbx(bx, 1, Lx,Ly,Lz);
    likely::TriCubicInterpolator interpolatedby(by, 1, Lx,Ly,Lz);
    likely::TriCubicInterpolator interpolatedbz(bz, 1, Lx,Ly,Lz);
    for(int c=0; c<Curve.Components.size(); c++)
    {
        int NP = Curve.Components[c].knotcurve.size();
        for(int s=0; s<NP; s++)
        {
            double tx= Curve.Components[c].knotcurve[s].tx;
            double ty= Curve.Components[c].knotcurve[s].ty;
            double tz= Curve.Components[c].knotcurve[s].tz;

            // step off the curve along the seifert surface, to get a b vector
            double xcoord = Curve.Components[c].knotcurve[s].xcoord+ firstpushdist*Curve.Components[c].knotcurve[s].omegax;
            double ycoord = Curve.Components[c].knotcurve[s].ycoord+ firstpushdist*Curve.Components[c].knotcurve[s].omegay;
            double zcoord = Curve.Components[c].knotcurve[s].zcoord+ firstpushdist*Curve.Components[c].knotcurve[s].omegaz;
            double tempbx = interpolatedbx(xcoord,ycoord,zcoord);
            double tempby = interpolatedby(xcoord,ycoord,zcoord);
            double tempbz = interpolatedbz(xcoord,ycoord,zcoord);

            double norm = sqrt(tempbx*tempbx + tempby*tempby + tempbz*tempbz);
            tempbx = tempbx/norm;
            tempby = tempby/norm;
            tempbz = tempbz/norm;
            Curve.Components[c].knotcurve[s].bx =tempbx;
            Curve.Components[c].knotcurve[s].by =tempby;
            Curve.Components[c].knotcurve[s].bz =tempbz;

            // I also want to project orthogonally and see if that makes a difference
            double projtempbx = (tempbx - (tempbx*tx + tempby*ty + tempbz*tz)*tx);
            double projtempby = (tempby - (tempbx*tx + tempby*ty + tempbz*tz)*ty);
            double projtempbz = (tempbz - (tempbx*tx + tempby*ty + tempbz*tz)*tz);
            norm = sqrt(projtempbx*projtempbx + projtempby*projtempby + projtempbz*projtempbz);
            projtempbx = projtempbx/norm;
            projtempby = projtempby/norm;
            projtempbz = projtempbz/norm;
            Curve.Components[c].knotcurve[s].projbx =projtempbx;
            Curve.Components[c].knotcurve[s].projby =projtempby;
            Curve.Components[c].knotcurve[s].projbz =projtempbz;
        }
    }

    // find the legendrian points
    likely::TriCubicInterpolator interpolatednx(nx, 1, Lx,Ly,Lz);
    likely::TriCubicInterpolator interpolatedny(ny, 1, Lx,Ly,Lz);
    likely::TriCubicInterpolator interpolatednz(nz, 1, Lx,Ly,Lz);
    for(int c=0; c<Curve.Components.size(); c++)
    {
        int NP = Curve.Components[c].knotcurve.size();
        for(int s=0; s<NP; s++)
        {
            double xcoord = Curve.Components[c].knotcurve[s].xcoord;
            double ycoord = Curve.Components[c].knotcurve[s].ycoord;
            double zcoord = Curve.Components[c].knotcurve[s].zcoord;
            double tx= Curve.Components[c].knotcurve[s].tx;
            double ty= Curve.Components[c].knotcurve[s].ty;
            double tz= Curve.Components[c].knotcurve[s].tz;
            double tempnx = interpolatednx(xcoord,ycoord,zcoord);
            double tempny = interpolatedny(xcoord,ycoord,zcoord);
            double tempnz = interpolatednz(xcoord,ycoord,zcoord);
            Curve.Components[c].knotcurve[s].ndott =tx*tempnx+ty*tempny+tz*tempnz;
        }
    }
    // linking and self linking numbers
    /*
    Curve.LinkingMatrix = vector<vector<double>>(Curve.Components.size(),vector<double>(Curve.Components.size(),99));
    for(int i=0; i<Curve.Components.size(); i++)
    {
        // the self linking case 
        int NP = Curve.Components[i].knotcurve.size();
        vector<knotpoint> PushOff(NP);
        double secondpushdist =2;
        for(int s=0; s<NP; s++)
        {
            PushOff[s].xcoord =Curve.Components[i].knotcurve[s].xcoord+secondpushdist*Curve.Components[i].knotcurve[s].bx;
            PushOff[s].ycoord =Curve.Components[i].knotcurve[s].ycoord+secondpushdist*Curve.Components[i].knotcurve[s].by;
            PushOff[s].zcoord =Curve.Components[i].knotcurve[s].zcoord+secondpushdist*Curve.Components[i].knotcurve[s].bz;
        }
        Curve.LinkingMatrix[i][i]=LinkingNumber(Curve.Components[i].knotcurve,PushOff);
        for(int j=0; j<Curve.Components.size(); j++)
        {
            // self linking
            if(i!=j)
            {
                // linking
               Curve.LinkingMatrix[i][j]= LinkingNumber(Curve.Components[i].knotcurve,Curve.Components[j].knotcurve);
            }
        }
    }
    */
    // smooth the bend zero for presentation purposes
    CurveSmoothing(Curve,filterlength);
}

// curve smoothing via a low pass filter
void CurveSmoothing(Link& Curve, int filterlength) 
{
    // for now this is just a moving average filter. In the future it could be a blackman windowed sinc etc. who knows!?
    for(int c=0; c<Curve.Components.size(); c++)
    {

        int NP = Curve.Components[c].knotcurve.size();
        vector<double> paddedcoord(NP+(filterlength-1),0);
        vector<double> smoothedpaddedcoord(NP+filterlength-1,0);
        int startindex = (filterlength-1)/2;
        int endindex = NP+startindex-1;
        for(int k=1; k<4; k++)
        {
            // reset the padded arrays
            for(int i=0; i<paddedcoord.size(); i++)
            {
                paddedcoord[i] =  0;
                smoothedpaddedcoord[i] = 0;
            }

            switch(k)
            {
                case 1 :
                    {for(int i=0; i<NP; i++) {paddedcoord[i+startindex] =  Curve.Components[c].knotcurve[i].xcoord ;} break;}
                case 2 :
                    {for(int i=0; i<NP; i++) {paddedcoord[i+startindex] =  Curve.Components[c].knotcurve[i].ycoord ;} break;}
                case 3 :
                    {for(int i=0; i<NP; i++) {paddedcoord[i+startindex] =  Curve.Components[c].knotcurve[i].zcoord ;} break;}
            }

            // im going to pad the ends of the array, not with zeros, but with the end values repeated
            for(int i=0;i<startindex;i++){paddedcoord[i]=paddedcoord[startindex];}
            for(int i=endindex;i<paddedcoord.size();i++){paddedcoord[i]=paddedcoord[endindex];}

            // okay apply the filter
            for(int i=startindex; i<=endindex; i++)
            {
                for(int j =-(filterlength-1)/2; j <=(filterlength-1)/2;j++)
                {
                    smoothedpaddedcoord[i] += (1/((double)filterlength))*paddedcoord[i+j];
                }
            }
            // copy back in
            switch(k)
            {
                case 1 :
                    {for(int i=startindex; i<=endindex; i++) {Curve.Components[c].knotcurve[i-startindex].xcoord = smoothedpaddedcoord[i] ;} break;}
                case 2 :
                    {for(int i=startindex; i<=endindex; i++) {Curve.Components[c].knotcurve[i-startindex].ycoord = smoothedpaddedcoord[i] ;} break;}
                case 3 :
                    { for(int i=startindex; i<=endindex; i++) {Curve.Components[c].knotcurve[i-startindex].zcoord = smoothedpaddedcoord[i] ;} break;}
            }

        }
    }
}

double my_minimisation_function(const gsl_vector* minimum, void* params)
{
    struct parameters* myparameters = (struct parameters *) params;
    double sphereradius = myparameters->sphereradius;
    likely::TriCubicInterpolator* interpolateducvmag = myparameters->ucvmag;
    gsl_vector* temppt = gsl_vector_alloc (3);
    gsl_vector* tempt = gsl_vector_alloc (3);
    gsl_vector* tempe1 = gsl_vector_alloc (3);
    gsl_vector* tempe2 = gsl_vector_alloc (3);
    gsl_vector_memcpy (temppt,myparameters->mypt);
    gsl_vector_memcpy (tempt,myparameters->t);
    gsl_vector_memcpy (tempe1,myparameters->e1);
    gsl_vector_memcpy (tempe2,myparameters->e2);

    // the minimisation happens on a hemisphere about our test point, chosen to be where in the half space where t points. the minima coordinates are x-y coords in this space. the point is given by x e1 + y e2 + sqrt(1- x^2 -y^2)t
    double x= gsl_vector_get (minimum, 0);
    double y= gsl_vector_get (minimum, 1);
    if(1-x*x-y*y <0){cout << "below zero!";}
    gsl_vector_scale(tempe1,sphereradius*x);
    gsl_vector_scale(tempe2,sphereradius*y);
    gsl_vector_scale(tempt,sphereradius*sqrt(1-x*x-y*y));
    gsl_vector_add(temppt,tempe1);
    gsl_vector_add(temppt,tempe2);
    gsl_vector_add(temppt,tempt);
    double px = gsl_vector_get(temppt, 0);
    double py = gsl_vector_get(temppt, 1);
    double pz = gsl_vector_get(temppt, 2);
    gsl_vector_free(temppt);
    gsl_vector_free(tempt);
    gsl_vector_free(tempe1);
    gsl_vector_free(tempe2);

    double value = ((*interpolateducvmag)(px,py,pz));
    return value;
}

bool KnotpointInMask(const knotpoint knotpoint, const bool* mask)
{
    // our mask is defined on gridpoints. We define being "in the mask" as being inside a cube with 1's on all its verices

    // the (x,y,z) value of the bottom left gridpoint corresponding to the cube our point is inside
    int idwn = (int) (knotpoint.xcoord);
    int jdwn = (int) (knotpoint.ycoord);
    int kdwn = (int) (knotpoint.zcoord);

    bool inmask = mask[pt(idwn,jdwn,kdwn)];
    for(int iinc=0;iinc<=1;iinc++)
    {
        for(int jinc=0;jinc<=1;jinc++)
        {
            for(int kinc=0;kinc<=1;kinc++)
            {
                inmask = inmask && mask[pt(idwn+iinc,jdwn+jinc,kdwn+kinc)];
            }
        }
    }
    return inmask;
}

// some little functions which go back and forth from point to index
int pt(const int k,const  int l,const  int m)       //convert i,j,k to single index
{
    return (m*Ly*Lx+l*Lx+k);
}

int incp(int i, int p, int N)    //increment i with p for periodic boundary
{
    if(i+p<0) return (N+i+p);
    else return ((i+p)%N);
}

int mod(int i, int N)    //my own mod fn
{
    if(i<0) return (N+i);
    else return ((i)%N);
}
int circularmod(int i, int N)    // mod i by N in a ciruclar fashion, ie wrapping around both in the +ve and -ve directions
{
    if(i<0) return (N - ((-i)%N))%N;
    else return i%N;
}

double x(int i)
{
    return (i+0.5-Lx/2.0);
}

double y(int j)
{
    return (j+0.5-Ly/2.0);
}

double z(int k)
{
    return (k+0.5-Lz/2.0);
}

double LinkingNumber(vector<knotpoint>& Curve1, vector<knotpoint>& Curve2)
{
    double dLinkingNumber=0;
    int NP1 = Curve1.size();
    int NP2 = Curve2.size();
    for(int s1=0; s1<NP1; s1++)
    {
        double x1 = Curve1[s1].xcoord;
        double y1 = Curve1[s1].ycoord;
        double z1 = Curve1[s1].zcoord;
        double dx1 = (Curve1[(s1+1)%NP1].xcoord - Curve1[s1].xcoord);
        double dy1 = (Curve1[(s1+1)%NP1].ycoord - Curve1[s1].ycoord);
        double dz1 = (Curve1[(s1+1)%NP1].zcoord - Curve1[s1].zcoord);
        for(int s2=0; s2<NP2; s2++)
        {
            double x2 = Curve2[s2].xcoord;
            double y2 = Curve2[s2].ycoord;
            double z2 = Curve2[s2].zcoord;
            double dx2 = (Curve2[(s2+1)%NP2].xcoord - Curve2[s2].xcoord);
            double dy2 = (Curve2[(s2+1)%NP2].ycoord - Curve2[s2].ycoord);
            double dz2 = (Curve2[(s2+1)%NP2].zcoord - Curve2[s2].zcoord);

            double num = (x2-x1)*(dy1*dz2-dz1*dy2)+(y2-y1)*(dz1*dx2-dx1*dz2)+(z2-z1)*(dx1*dy2-dy1*dx2);
            double denom = std::pow( (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1) , 1.5);
            dLinkingNumber += (num/denom);
        }
    }
    return dLinkingNumber;
}


