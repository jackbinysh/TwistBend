#include "SolidAngleLink.h"
#include "TriCubicInterpolator.h"
#include <math.h>
#include <string.h>
#include <omp.h>

// computes the solid angle the link presents from a viewpoint in its complement
// this version sweeps the solid angle out by moving the link in from being initially asymptotically far away
double SolidAngleCalc(const Link& Curve, const viewpoint& View)
{
    double totalomega = 0;
    double mindist = 1000.0; // scroll waves
    for(int i=0; i<Curve.Components.size(); i++)
    {
        double Integral = 0;
        int NP = Curve.Components[i].knotcurve.size();

        // first define and check a choice of asymptotic direction
        double ndotnmin = 1.0;
        double ndotnmax = -1.0;
        int smin;
        // define the asymptotic direction -- ninfty -- z-axis by default (in lower half space)
        double ninftyx = 0.0;
        double ninftyy = 0.0;
        double ninftyz = 1.0;
        if (View.zcoord>0) {ninftyz = -1.0;} // minus z in the upper half space
        for (int s=0; s<NP; s++)
        {
            // define the view vector -- n = (Curve - View)/|Curve - View|
            double viewx = Curve.Components[i].knotcurve[s].xcoord - View.xcoord;
            double viewy = Curve.Components[i].knotcurve[s].ycoord - View.ycoord;
            double viewz = Curve.Components[i].knotcurve[s].zcoord - View.zcoord;
            double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);
            double ndotninfty = viewz*ninftyz/dist;
            if (ndotninfty<ndotnmin) {ndotnmin = ndotninfty; smin = s;}
            if (ndotninfty>ndotnmax) {ndotnmax = ndotninfty;}
        }
        if (ndotnmin < -0.98) // check if a threshold is exceeded -- value can be changed
        {
            if (ndotnmax < 0.98) {ninftyz = -ninftyz;} // flip direction
            else                                       // unless another threshold is exceeded -- value can be changed
            {
                ninftyz = 0.0;
                ninftyx = Curve.Components[i].knotcurve[smin].ty;    // set an orthogonal direction -- not guaranteed to be a good choice
                ninftyy = -Curve.Components[i].knotcurve[smin].tx;
                double norm = sqrt(ninftyx*ninftyx + ninftyy*ninftyy);
                ninftyx /= norm;
                ninftyy /= norm;
                //	  if (View.xcoord>0) {ninftyx = -ninftyx;} // could be picky about signs -- old code here, beware !!
            }
        }
        for (int s=0; s<NP; s++)
        {
            // define the view vector -- n = (Curve - View)/|Curve - View|
            double viewx = Curve.Components[i].knotcurve[s].xcoord - View.xcoord;
            double viewy = Curve.Components[i].knotcurve[s].ycoord - View.ycoord;
            double viewz = Curve.Components[i].knotcurve[s].zcoord - View.zcoord;
            double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);
            double ndotninfty = viewx*ninftyx + viewy*ninftyy + viewz*ninftyz;
            double tx = Curve.Components[i].knotcurve[s].tx;
            double ty = Curve.Components[i].knotcurve[s].ty;
            double tz = Curve.Components[i].knotcurve[s].tz;
            // trapezium rule quadrature
            double ds = 0.5*(Curve.Components[i].knotcurve[s].length+Curve.Components[i].knotcurve[incp(s,-1,NP)].length);
            // and here's the integrand
            Integral += (ds/dist)*(ninftyz*(ty*viewx-tx*viewy)+ninftyx*(tz*viewy-ty*viewz)+ninftyy*(tx*viewz-tz*viewx))/(dist + ndotninfty);
        }

        totalomega += Integral;
    }
    return totalomega;
}

void ComputeSolidAngle(double* phi, const Link& Curve)
{
    double SolidAngle;
    viewpoint Point;
    int i,j,k;
    for(i=0; i<Lx; i++)
    {
        for(j=0; j<Ly; j++)
        {
            for(k=0; k<Lz; k++)
            {
                int n = pt(i,j,k);
                Point.xcoord = i;
                Point.ycoord = j ;
                Point.zcoord = k;

                SolidAngle = SolidAngleCalc(Curve,Point);
                // put in the interval [-2pi,2pi]
                while(SolidAngle>2*M_PI) SolidAngle -= 4*M_PI;
                while(SolidAngle<-2*M_PI) SolidAngle += 4*M_PI;
                phi[n]= SolidAngle;

            }
        }
    }

}

// find the difference between the framing we have given the curve, and the framing provided by the solid angle function.
void ComputeSolidAngleFraming(double* phi,Link&  Curve)
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
            double mintheta = -10;
            double finalvx=0 ;
            double finalvy=0;
            double finalvz=0 ;
            for(int q=0;q<NumTestPoints;q++)
            {
                if(fabs(testphis[q]-phi0) < fabs(minphi-phi0))
                {
                    minphi = testphis[q];
                    double theta = ( ((double)q) / (double(NumTestPoints)) )*2*M_PI;
                    mintheta = theta;
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

