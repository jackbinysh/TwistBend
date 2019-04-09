#include "SolidAngle.h"    
#include "Geometry.h"
#include "InputOutput.h"
#include "Constants.h"
#include "../twistbend.h"
#include "../TriCubicInterpolator.h"
#include <math.h>
#include <cmath>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

double ComputeSolidAngleOnePoint(const Link& Curve, const viewpoint& View)
{
    double totalomega = 0.0;
    for(int i=0; i<Curve.Components.size(); i++)
    {
        int NP = Curve.Components[i].knotcurve.size();

        // the user's settings for n_infty
        double ninftyx = Initialninftyx;
        double ninftyy = Initialninftyy;
        double ninftyz = Initialninftyz;
        double mag = sqrt(ninftyx*ninftyx+ninftyy*ninftyy+ninftyz*ninftyz);
        if(std::abs(mag-1)>0.0001)
        {
            cout << "the magnitude of your n_infty vector isn't close to 1 - did you normalise? I've done it for you but best double check it's really the direction you want \n";
        }
        ninftyx /= mag;
        ninftyy /= mag;
        ninftyz /= mag;

        // can we proceed to the calculation? Initially the answer is no
        bool thresholdokay = false;
        // how many times have we needed to pick a new vector?
        int numtimesthresholdexceeded = 0;
        while ( thresholdokay == false)
        {
            double ndotnmin = 1.0;
            double ndotnmax = -1.0;
            int smin;
            for (int s=0; s<NP; s++)
            {
                double viewx = Curve.Components[i].knotcurve[s].xcoord - View.xcoord;
                double viewy = Curve.Components[i].knotcurve[s].ycoord - View.ycoord;
                double viewz = Curve.Components[i].knotcurve[s].zcoord - View.zcoord;
                double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);
                double ndotninfty = (viewx*ninftyx+viewy*ninftyy+viewz*ninftyz)/dist;
                if (ndotninfty<ndotnmin) {ndotnmin = ndotninfty; smin = s;}
                if (ndotninfty>ndotnmax) {ndotnmax = ndotninfty;}
            }

            thresholdokay = true;
            // is ninfty an okay choice ?
            if (ndotnmin < Threshold_ninf_dot_n)
            {
                // if i sent ninfty -> -ninfty would it be okay?
                if (ndotnmax > -Threshold_ninf_dot_n)
                {
                    // problem case - we need to choose a genuinely new direction. We do so from a list of hardcoded vectors.
                    // in practice I've never hit the end of this list, but if we do, just give up.
                    thresholdokay = false;
                    ninftyx = hardcodedvectors[3*numtimesthresholdexceeded];
                    ninftyy = hardcodedvectors[3*numtimesthresholdexceeded+1];
                    ninftyz = hardcodedvectors[3*numtimesthresholdexceeded+2];

                    numtimesthresholdexceeded ++;
                    if(numtimesthresholdexceeded == numhardcodedvectors) {thresholdokay = true;} // give up jack
                }
                else
                {
                    // flippings fine, do that
                    ninftyx = -ninftyx;
                    ninftyy = -ninftyy;
                    ninftyz = -ninftyz;
                }
            }
        }

        // okay we have an acceptable n_infty - now do the integration. Simple trapezium rule quadrature
        double Integral = 0;
        for (int s=0; s<NP; s++)
        {
            double viewx = Curve.Components[i].knotcurve[s].xcoord - View.xcoord;
            double viewy = Curve.Components[i].knotcurve[s].ycoord - View.ycoord;
            double viewz = Curve.Components[i].knotcurve[s].zcoord - View.zcoord;
            double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);
            double ndotninfty = viewx*ninftyx + viewy*ninftyy + viewz*ninftyz;
            double tx = Curve.Components[i].knotcurve[s].tx;
            double ty = Curve.Components[i].knotcurve[s].ty;
            double tz = Curve.Components[i].knotcurve[s].tz;
            double ds = 0.5*(Curve.Components[i].knotcurve[s].length+Curve.Components[i].knotcurve[incp(s,-1,NP)].length);
            Integral += (ds/dist)*(ninftyz*(ty*viewx-tx*viewy)+ninftyx*(tz*viewy-ty*viewz)+ninftyy*(tx*viewz-tz*viewx))/(dist + ndotninfty);
        }

        totalomega += Integral;
    }
    while(totalomega>4*M_PI) totalomega -= 4*M_PI;
    while(totalomega<0) totalomega += 4*M_PI;
    return totalomega;
}

void ComputeSolidAngleAllPoints(const Link& Curve,double* omega)
{

#pragma omp parallel default(none) shared(omega,Curve,cout)
    {
        viewpoint Point;
        int progressbarcounter =0;

#pragma omp for
        for (int k=0; k<Lz; k++)
        {
            for (int j=0; j<Ly; j++)
            {
                for (int i=0; i<Lx; i++)
                {
                    Point.xcoord = i;
                    Point.ycoord = j ;
                    Point.zcoord = k;

                    double SolidAngle = ComputeSolidAngleOnePoint(Curve,Point);

                    int n = pt(i,j,k);
                    omega[n]=SolidAngle;
                }
            }
            // progress bar!
#ifdef _OPENMP
            if(omp_get_thread_num() == 0)
            {
                if ((  ( ( ((double)k)*omp_get_num_threads() )/Lz ) * 100 ) > progressbarcounter)
                {
                    cout<< progressbarcounter << "%" << endl;
                    progressbarcounter += 10;
                }
            }
#endif
#ifndef _OPENMP
            if (( ((double)(k)/Lz)*100 )> progressbarcounter )
            {
                cout<< progressbarcounter << "%" << endl;
                progressbarcounter += 10;
            }
#endif
        }
    }
}
// find the difference between the framing we have given the curve, and the framing provided by the solid angle function.
void ComputeSolidAngleFraming(Link&  Curve, double* phi)
{
    // settings
    int NumTestPoints=100;
    std::vector< double > testphis(NumTestPoints);
    // radius of test circle
    double r = 3;
    double phi0 = 0;

    // Im going to want phi values around the filament. for this, we construct an interpolator object
    likely::TriCubicInterpolator interpolatedphi(phi,1,Lx,Ly,Lz);

    for(int i=0; i<Curve.Components.size(); i++)
    {
        int NP = Curve.Components[i].knotcurve.size();

        for (int s=0; s<NP; s++)
        {
            // grab the point, compute the framing
            double xcoord= Curve.Components[i].knotcurve[s].xcoord;
            double ycoord= Curve.Components[i].knotcurve[s].ycoord;
            double zcoord= Curve.Components[i].knotcurve[s].zcoord;
            double ax= Curve.Components[i].knotcurve[s].e1x;
            double ay= Curve.Components[i].knotcurve[s].e1y;
            double az= Curve.Components[i].knotcurve[s].e1z;
            double tx= Curve.Components[i].knotcurve[s].tx;
            double ty= Curve.Components[i].knotcurve[s].ty;
            double tz= Curve.Components[i].knotcurve[s].tz;
            double tcax = ty*az-tz*ay;
            double tcay = tz*ax-tx*az;
            double tcaz = tx*ay-ty*ax;

            // okay, we have the frame at this point. lets walk in a small circle
            // around the filament and just find the rough value of phi closest to phi0.
            for(int q=0;q<NumTestPoints;q++)
            {
                double theta = ( ((double)q) / (double(NumTestPoints)) )*2*M_PI;
                double vx = r*(cos(theta)*ax + sin(theta)*tcax);
                double vy = r*(cos(theta)*ay + sin(theta)*tcay);
                double vz = r*(cos(theta)*az + sin(theta)*tcaz);
                testphis[q]= interpolatedphi(xcoord+vx,ycoord+vy,zcoord+vz);
            }

            // okay we have our test values. now find the one cloest to -2*pi and the one closest to 2*pi -
            // aka the min and max. I want to exclude this region from the search, as its where the cut is.
            double minphi = 20;
            int minq = -1;
            for(int q=0;q<NumTestPoints;q++)
            {
                if(testphis[q]<minphi)
                {
                    minphi = testphis[q];
                    minq = q;
                }
            }
            double maxphi = -1;
            int maxq = -1;
            for(int q=0;q<NumTestPoints;q++)
            {
                if(testphis[q]>maxphi)
                {
                    maxphi = testphis[q];
                    maxq = q;
                }
            }

            // which region should we exclude? we need to pick one of two arcs of a circle. we choose the
            // one with the smaller length.

            if( mod((maxq - minq),NumTestPoints) > mod((minq - maxq),NumTestPoints))
            {
                for(int q=maxq;q<maxq + mod((minq - maxq),NumTestPoints);q++){testphis[mod(q,NumTestPoints)] = INFINITY;}
            }
            if( mod((maxq - minq),NumTestPoints) < mod((minq - maxq),NumTestPoints))
            {
                for(int q=minq;q<minq + mod((maxq - minq),NumTestPoints);q++){testphis[mod(q,NumTestPoints)] = INFINITY;}
            }

            minphi = -10;
            double finalvx=0 ;
            double finalvy=0;
            double finalvz=0 ;
            for(int q=0;q<NumTestPoints;q++)
            {
                if(fabs(testphis[q]-phi0) < fabs(minphi-phi0))
                {
                    minphi = testphis[q];
                    double theta = ( ((double)q) / (double(NumTestPoints)) )*2*M_PI;
                    double vx = r*(cos(theta)*ax + sin(theta)*tcax);
                    double vy = r*(cos(theta)*ay + sin(theta)*tcay);
                    double vz = r*(cos(theta)*az + sin(theta)*tcaz);
                    finalvx =vx;
                    finalvy =vy;
                    finalvz =vz;
                }
            }
            // normalize so the output is a unit vector - our search circle may have a radius which is not 1.
            double norm = sqrt(finalvx*finalvx + finalvy*finalvy+ finalvz*finalvz);
            finalvx = finalvx/norm;
            finalvy = finalvy/norm;
            finalvz = finalvz/norm;

            Curve.Components[i].knotcurve[s].omegax = finalvx;
            Curve.Components[i].knotcurve[s].omegay = finalvy;
            Curve.Components[i].knotcurve[s].omegaz = finalvz;
        }
    }
}

double ComputeDistanceOnePoint(const Link& Curve, const viewpoint& View)
{
    double mindist = 10000.0;
    for(int i=0; i<Curve.Components.size(); i++)
    {
        int NP = Curve.Components[i].knotcurve.size();
        for (int s=0; s<NP; s++)
        {
            double viewx = Curve.Components[i].knotcurve[s].xcoord - View.xcoord;
            double viewy = Curve.Components[i].knotcurve[s].ycoord - View.ycoord;
            double viewz = Curve.Components[i].knotcurve[s].zcoord - View.zcoord;
            double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);
            if (dist < mindist) {
                mindist = dist;
            }
        }
    }
    return mindist;
}

// not adapted to distinct components yet!!
double ComputeLongitudinalPhase(const Link& Curve, const viewpoint& View)
{
    double phase = 0.0;
    double mindist = 10000.0;
    for(int i=0; i<Curve.Components.size(); i++)
    {
        int NP = Curve.Components[i].knotcurve.size();
        for (int s=0; s<NP; s++)
        {
            double viewx = Curve.Components[i].knotcurve[s].xcoord - View.xcoord;
            double viewy = Curve.Components[i].knotcurve[s].ycoord - View.ycoord;
            double viewz = Curve.Components[i].knotcurve[s].zcoord - View.zcoord;
            double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);
            if (dist < mindist) {
                mindist = dist;
                phase = 2.0*M_PI*s/NP;
            }
        }
    }
    return phase;
}
