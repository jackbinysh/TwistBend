/* Fitzhugh-Nagumo reaction diffusion simulation with arbitrary vortex lines
   OPENMP VERSION
   Created by Carl Whitfield
   Last modified 03/01/17

   Operational order of the code:
   1) The code takes as input an stl file (defined in knot_filename) which defines an orientable surface with a boundary.
   2) This surface is scaled to fill a box of size xmax x ymax x zmax.
   3) A nunmerical integral is performed to calculate a phase field (phi_calc) on the 3D grid which winds around the boundary of the surface.
   4) This is then used to initialise the Fitzhugh-Nagumo set of partial differential equations such that on initialisation (uv_initialise):
   u = 2cos(phi) - 0.4        v = sin(phi) - 0.4
   The pde's used are
   dudt = (u - u^3/3 - v)/epsilon + Del^2 u
   dvdt = epsilon*(u + beta - gam v)
   5) The update method is Runge-Kutta fourth order (update_uv) unless RK4 is set to 0, otherwise Euler forward method is used.
   6) A parametric curve for the knot is found at each unit T


   The parameters epsilon, beta and gam are set to give rise to scroll waves (see Sutcliffe, Winfree, PRE 2003 and Maucher Sutcliffe PRL 2016) which eminate from a closed curve, on initialisation this curve corresponds to the boundary of the surface in the stl file.

   See below for various options to start the code from previous outputted data.*/
#include "FN_Knot.h"    //contains user defined variables for the simulation, and the parameters used 
#include "Initialisation.h"    //contains user defined variables for the simulation, and the parameters used
#include "TriCubicInterpolator.h"    //contains user defined variables for the simulation, and the parameters used
#include "ReadingWriting.h"    //contains user defined variables for the simulation, and the parameters used
#include <omp.h>
#include <math.h>
#include <string.h>
#include <algorithm>
//includes for the signal processing
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

int main (void)
{
    Griddata griddata;
    griddata.Nx = initialNx;
    griddata.Ny = initialNy;
    griddata.Nz = initialNz;
    griddata.h = initialh;
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    // all major allocations are here
    // the main data storage arrays, contain info associated with the grid
    vector<double>phi(Nx*Ny*Nz);  //scalar potential
    vector<double>u(Nx*Ny*Nz);
    vector<double>v(Nx*Ny*Nz);
    vector<double>ucvx(Nx*Ny*Nz);
    vector<double>ucvy(Nx*Ny*Nz);
    vector<double>ucvz(Nx*Ny*Nz);
    vector<double>ucvmag(Nx*Ny*Nz);// mod(grad u cross grad v)
    vector<double>ku(4*Nx*Ny*Nz);
    vector<double>kv(4*Nx*Ny*Nz);
    // objects to hold information about the knotcurve we find, andthe surface we read in
    vector<knotcurve > knotcurves; // a structure containing some number of knot curves, each curve a list of knotpoints
    vector<knotcurve > knotcurvesold; // a structure containing some number of knot curves, each curve a list of knotpoints
    vector<triangle> knotsurface;    //structure for storing knot surface coordinates
    // sensor points we output u values at
    viewpoint sensorpoint;
    // GSL initialization
    const gsl_multimin_fminimizer_type *Type;
    gsl_multimin_fminimizer *minimizerstate;
    Type = gsl_multimin_fminimizer_nmsimplex2;
    minimizerstate = gsl_multimin_fminimizer_alloc (Type,2);

    // setting things from globals
    int starttime = 0;
    int FrequentKnotplotPrintIteration = (int)(FrequentKnotplotPrintTime/dtime);
    int VelocityKnotplotPrintIteration = (int)(VelocityKnotplotPrintTime/dtime);
    int InitialSkipIteration = (int)(InitialSkipTime/dtime);
    int UVPrintIteration = (int)(UVPrintTime/dtime);
    sensorpoint.xcoord = sensorxcoord ;
    sensorpoint.ycoord = sensorycoord ;
    sensorpoint.zcoord = sensorzcoord ;
    // initialising timers
    time_t rawtime;
    time (&rawtime);
    struct tm * timeinfo;

    // INITIALISATION

    switch(option)
    {
        case FROM_UV_FILE:
            {
                cout << "Reading input file...\n";
                if(uvfile_read(u,v,ku,kv, ucvx,ucvy,ucvz,ucvmag,griddata)){return 1;}
                // get the start time -  we hack this together as so:
                // the filename looks like uv_plotxxx.vtk, we want the xxx. so we find the t, find the ., and grab everyting between
                string number = B_filename.substr(B_filename.find('t')+1,B_filename.find('.')-B_filename.find('t')-1);
                starttime = atoi(number.c_str());
                break;
            }
        case FROM_FUNCTION:
            {
                phi_calc_manual(phi,griddata);
                cout << "Calculating u and v...\n";
                uv_initialise(phi,u,v,griddata);
                break;
            }
        case FROM_SURFACE_FILE:
            {
                init_from_surface_file(knotsurface);
                phi_calc_surface(phi,knotsurface,griddata);
                cout << "Calculating u and v...\n";
                uv_initialise(phi,u,v,griddata);
                break;
            }
        case FROM_CURVE_FILE:
            {
                Link Curve;
                InitialiseFromFile(Curve);
                cout << "calculating the solid angle..." << endl;
                phi_calc_curve(phi,Curve,griddata);
                cout << "Calculating u and v...\n";
                uv_initialise(phi,u,v,griddata);
            }

    }

    // UPDATE
    cout << "Updating u and v...\n";

    double CurrentTime = starttime;
    int CurrentIteration = (int)(CurrentTime/dtime);
#pragma omp parallel default(none) shared (u,v,ku,kv,ucvx, CurrentIteration,InitialSkipIteration,FrequentKnotplotPrintIteration,UVPrintIteration,VelocityKnotplotPrintIteration,ucvy, ucvz,ucvmag,cout, rawtime, starttime, timeinfo,CurrentTime, knotcurves,knotcurvesold,minimizerstate,griddata,sensorpoint)
    {
        while(CurrentTime <= TTime)
        {
#pragma omp single
            {

                // its useful to have an oppurtunity to print the knotcurve, without doing the velocity tracking, whihc doesnt work too well if we go more frequenclty
                // than a cycle
                if( ( CurrentIteration >= InitialSkipIteration ) && ( CurrentIteration%FrequentKnotplotPrintIteration==0) )
                {
                    cout << "T = " << CurrentTime << endl;
                    time (&rawtime);
                    timeinfo = localtime (&rawtime);
                    cout << "current time \t" << asctime(timeinfo) << "\n";
                    crossgrad_calc(u,v,ucvx,ucvy,ucvz,ucvmag,griddata); //find Grad u cross Grad v

                    find_knot_properties(ucvx,ucvy,ucvz,ucvmag,u,knotcurves,CurrentTime,minimizerstate ,griddata);      //find knot curve and twist and writhe
                    print_knot(CurrentTime, knotcurves, griddata);

                    print_sensor_point(CurrentTime,sensorpoint,u,griddata);
                }

                // run the curve tracing, and find the velocity of the one we previously stored, then print that previous one
                if( ( CurrentIteration > InitialSkipIteration ) && ( CurrentIteration%VelocityKnotplotPrintIteration==0) )
                {
                    crossgrad_calc(u,v,ucvx,ucvy,ucvz,ucvmag,griddata); //find Grad u cross Grad v

                    find_knot_properties(ucvx,ucvy,ucvz,ucvmag,u,knotcurves,CurrentTime,minimizerstate ,griddata);      //find knot curve and twist and writhe
                    if(!knotcurvesold.empty())
                    {
                        find_knot_velocity(knotcurves,knotcurvesold,griddata,VelocityKnotplotPrintTime);
                        print_knot(CurrentTime - VelocityKnotplotPrintTime , knotcurvesold, griddata);
                    }
                    knotcurvesold = knotcurves;

                }

                // print the UV, and ucrossv data
                if(CurrentIteration%UVPrintIteration==0)
                {
                    crossgrad_calc(u,v,ucvx,ucvy,ucvz,ucvmag,griddata); //find Grad u cross Grad v
                    print_uv(u,v,ucvx,ucvy,ucvz,ucvmag,CurrentTime,griddata);
                }
                //though its useful to have a double time, we want to be careful to avoid double round off accumulation in the timer
                CurrentIteration++;
                CurrentTime  = ((double)(CurrentIteration) * dtime);
            }
            uv_update(u,v,ku,kv,griddata);
        }
    }
    return 0;
}


void uv_initialise(vector<double>&phi, vector<double>&u, vector<double>&v, const Griddata& griddata)
{
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    int n;

    for(n=0; n<Nx*Ny*Nz; n++)
    {
        // turns out we want -solidangle in order that the initialisation has the same sense - "clockwise vortex rotation in time" as our initialisation curve.
        double minusphi = -phi[n];
        u[n] = (2*cos(minusphi) - 0.4);
        v[n] = (sin(minusphi) - 0.4);
    }
}

void crossgrad_calc( vector<double>&u, vector<double>&v, vector<double>&ucvx, vector<double>&ucvy, vector<double>&ucvz, vector<double>&ucvmag,const Griddata& griddata)
{
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    double h = griddata.h;
    int i,j,k,n,kup,kdown;
    double dxu,dyu,dzu,dxv,dyv,dzv;
    for(i=0;i<Nx;i++)
    {
        for(j=0; j<Ny; j++)
        {
            for(k=0; k<Nz; k++)   //Central difference
            {
                kup = gridinc(k,1,Nz,2);
                kdown = gridinc(k,-1,Nz,2);
                dxu = 0.5*(u[pt(gridinc(i,1,Nx,0),j,k,griddata)]-u[pt(gridinc(i,-1,Nx,0),j,k,griddata)])/h;
                dxv = 0.5*(v[pt(gridinc(i,1,Nx,0),j,k,griddata)]-v[pt(gridinc(i,-1,Nx,0),j,k,griddata)])/h;
                dyu = 0.5*(u[pt(i,gridinc(j,1,Ny,1),k,griddata)]-u[pt(i,gridinc(j,-1,Ny,1),k,griddata)])/h;
                dyv = 0.5*(v[pt(i,gridinc(j,1,Ny,1),k,griddata)]-v[pt(i,gridinc(j,-1,Ny,1),k,griddata)])/h;
                dzu = 0.5*(u[pt(i,j,kup,griddata)]-u[pt(i,j,kdown,griddata)])/h;
                dzv = 0.5*(v[pt(i,j,kup,griddata)]-v[pt(i,j,kdown,griddata)])/h;
                //          dxu =(-u[pt(gridinc(i,2,Nx,0),j,k,griddata)]+8*u[pt(gridinc(i,1,Nx,0),j,k,griddata)]-8*u[pt(gridinc(i,-1,Nx,0),j,k,griddata)]+u[pt(gridinc(i,-2,Nx,0),j,k,griddata)])/(12*h);
                //          dxv =(-v[pt(gridinc(i,2,Nx,0),j,k,griddata)]+8*v[pt(gridinc(i,1,Nx,0),j,k,griddata)]-8*v[pt(gridinc(i,-1,Nx,0),j,k,griddata)]+v[pt(gridinc(i,-2,Nx,0),j,k,griddata)])/(12*h);
                //          dyu =(-u[pt(gridinc(j,2,Ny,1),j,k,griddata)]+8*u[pt(gridinc(j,1,Ny,1),j,k,griddata)]-8*u[pt(gridinc(j,-1,Ny,1),j,k,griddata)]+u[pt(gridinc(j,-2,Ny,1),j,k,griddata)])/(12*h);
                //          dyv =(-v[pt(gridinc(j,2,Ny,1),j,k,griddata)]+8*v[pt(gridinc(j,1,Ny,1),j,k,griddata)]-8*v[pt(gridinc(j,-1,Ny,1),j,k,griddata)]+v[pt(gridinc(j,-2,Ny,1),j,k,griddata)])/(12*h);
                //          dzu =(-u[pt(gridinc(k,2,Nz,2),j,k,griddata)]+8*u[pt(gridinc(k,1,Nz,2),j,k,griddata)]-8*u[pt(gridinc(k,-1,Nz,2),j,k,griddata)]+u[pt(gridinc(k,-2,Nz,2),j,k,griddata)])/(12*h);
                //          dzv =(-v[pt(gridinc(k,2,Nz,2),j,k,griddata)]+8*v[pt(gridinc(k,1,Nz,2),j,k,griddata)]-8*v[pt(gridinc(k,-1,Nz,2),j,k,griddata)]+v[pt(gridinc(k,-2,Nz,2),j,k,griddata)])/(12*h);
                n = pt(i,j,k,griddata);
                ucvx[n] = dyu*dzv - dzu*dyv;
                ucvy[n] = dzu*dxv - dxu*dzv;    //Grad u cross Grad v
                ucvz[n] = dxu*dyv - dyu*dxv;
                ucvmag[n] = sqrt(ucvx[n]*ucvx[n] + ucvy[n]*ucvy[n] + ucvz[n]*ucvz[n]);
            }
        }
    }
}

void find_knot_properties( vector<double>&ucvx, vector<double>&ucvy, vector<double>&ucvz, vector<double>& ucvmag,vector<double>&u,vector<knotcurve>& knotcurves,double t, gsl_multimin_fminimizer* minimizerstate, const Griddata& griddata)
{
    // first thing, clear the knotcurve object before we begin writing a new one
    knotcurves.clear(); //empty vector with knot curve points

    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    double h = griddata.h;

    // initialise the tricubic interpolator for ucvmag
    likely::TriCubicInterpolator interpolateducvmag(ucvmag, h, Nx,Ny,Nz);

    int c =0;
    bool knotexists = true;
    vector<int> marked(Nx*Ny*Nz,0);

    while(knotexists)
    {
        double   ucvmax = -1.0; // should always be +ve, so setting it to an initially -ve # means it always gets written to once.
        int n,i,j,k,imax,jmax,kmax;
        for(i=0;i<Nx;i++)
        {
            for(j=0; j<Ny; j++)
            {
                for(k=0; k<Nz; k++)   //Central difference
                {
                    n = pt(i,j,k,griddata);
                    if( ucvmag[n] > ucvmax && marked[n]==0)
                    {
                        ucvmax = ucvmag[n];
                        imax = i;
                        jmax = j;
                        kmax=k;

                    }
                }
            }
        }

        if(ucvmax<0.45) knotexists = false;

        if(knotexists)
        {
            knotcurves.push_back(knotcurve() );
            knotcurves[c].knotcurve.push_back(knotpoint());
            knotcurves[c].knotcurve[0].xcoord=x(imax,griddata);
            knotcurves[c].knotcurve[0].ycoord=y(jmax,griddata);
            knotcurves[c].knotcurve[0].zcoord=z(kmax,griddata);

            int idwn,jdwn,kdwn, modidwn, modjdwn, modkdwn,m,iinc,jinc,kinc;
            double ucvxs, ucvys, ucvzs, graducvx, graducvy, graducvz, prefactor, xd, yd ,zd, fx, fy, fz, xdiff, ydiff, zdiff;
            int s=1;
            bool finish=false;
            // we will discard the first few points from the knot, using this flag
            bool burnin=true;
            // if the curve we are tracing terminates at a boundary, for now we will just strike it from the record, entering a cleanup mode that marks the broken curve as "do not touch"
            int boundaryhits = 0;
            // don't allow hits all over the boundary
            int cooldowntimer = 0;
            /*calculate local direction of grad u x grad v (the tangent to the knot curve) at point s-1, then move to point s by moving along tangent + unit confinement force*/
            while (finish==false)
            {

                /**Find nearest gridpoint**/
                idwn = (int) ((knotcurves[c].knotcurve[s-1].xcoord/h) - 0.5 + Nx/2.0);
                jdwn = (int) ((knotcurves[c].knotcurve[s-1].ycoord/h) - 0.5 + Ny/2.0);
                kdwn = (int) ((knotcurves[c].knotcurve[s-1].zcoord/h) - 0.5 + Nz/2.0);
                // idwn etc can be off the actual grid , into "ghost" grids around the real one. this is useful for knotcurve tracing over periodic boundaries
                // but we also need the corresponding real grid positions!
                modidwn = circularmod(idwn,Nx);
                modjdwn = circularmod(jdwn,Ny);
                modkdwn = circularmod(kdwn,Nz);

                // if we have hit a boundary, dont go off grid - rather, set a cleanup flag and stick to the grid egde
                if(cooldowntimer>0) {cooldowntimer--;}
                if(BoundaryType==ALLREFLECTING && cooldowntimer ==0)
                {
                    if(idwn <=0){idwn=0; boundaryhits++; cooldowntimer=5;}
                    if(idwn >=Nx-1){idwn=Nx-1; boundaryhits++;cooldowntimer=5;}
                    if(jdwn <=0){jdwn=0; boundaryhits++;cooldowntimer=5;}
                    if(jdwn >=Ny-1){jdwn=Ny-1; boundaryhits++;cooldowntimer=5;}
                    if(kdwn <=0){kdwn=0; boundaryhits++;cooldowntimer=5;}
                    if(kdwn >=Nz-1){kdwn=Nz-1; boundaryhits++;cooldowntimer=5;}
                }
                if(BoundaryType==ZPERIODIC && cooldowntimer==0)
                {
                    if(idwn <=0){idwn=0; boundaryhits++;cooldowntimer=5;}
                    if(idwn >Nx-1){idwn=Nx-1; boundaryhits++;cooldowntimer=5;}
                    if(jdwn <=0){jdwn=0; boundaryhits++;cooldowntimer=5;}
                    if(jdwn >Ny-1){jdwn=Ny-1; boundaryhits++;cooldowntimer=5;}
                }

                ucvxs=0;
                ucvys=0;
                ucvzs=0;
                /*curve to gridpoint down distance*/
                xd = (knotcurves[c].knotcurve[s-1].xcoord - x(idwn,griddata))/h;
                yd = (knotcurves[c].knotcurve[s-1].ycoord - y(jdwn,griddata))/h;
                zd = (knotcurves[c].knotcurve[s-1].zcoord - z(kdwn,griddata))/h;
                for(m=0;m<8;m++)  //linear interpolation from 8 nearest neighbours
                {
                    /* Work out increments*/
                    iinc = m%2;
                    jinc = (m/2)%2;
                    kinc = (m/4)%2;
                    /*Loop over nearest points*/
                    i = gridinc(modidwn, iinc, Nx,0);
                    j = gridinc(modjdwn, jinc, Ny,1);
                    k = gridinc(modkdwn,kinc, Nz,2);
                    prefactor = (1-iinc + pow(-1,1+iinc)*xd)*(1-jinc + pow(-1,1+jinc)*yd)*(1-kinc + pow(-1,1+kinc)*zd);
                    /*interpolate grad u x grad v over nearest points*/
                    ucvxs += prefactor*ucvx[pt(i,j,k,griddata)];
                    ucvys += prefactor*ucvy[pt(i,j,k,griddata)];
                    ucvzs += prefactor*ucvz[pt(i,j,k,griddata)];
                }
                double norm = sqrt(ucvxs*ucvxs + ucvys*ucvys + ucvzs*ucvzs);
                ucvxs = ucvxs/norm; //normalise
                ucvys = ucvys/norm; //normalise
                ucvzs = ucvzs/norm; //normalise

                // if we have hit a boundary, we want to back up along the curve instead
                if(boundaryhits==1)
                {
                    ucvxs *= -1 ;
                    ucvys *= -1 ;
                    ucvzs *= -1;
                }

                // okay we have our first guess, move forward in this direction
                // we actually want to walk in the direction gradv cross gradu - that should be our +ve tangent,
                // so that the rotation sense of the curve is positive. Get this we - signs below.
                double testx = knotcurves[c].knotcurve[s-1].xcoord - h*ucvxs;
                double testy = knotcurves[c].knotcurve[s-1].ycoord - h*ucvys;
                double testz = knotcurves[c].knotcurve[s-1].zcoord - h*ucvzs;

                // now get the grad at this point
                idwn = (int) ((testx/h) - 0.5 + Nx/2.0);
                jdwn = (int) ((testy/h) - 0.5 + Ny/2.0);
                kdwn = (int) ((testz/h) - 0.5 + Nz/2.0);
                modidwn = circularmod(idwn,Nx);
                modjdwn = circularmod(jdwn,Ny);
                modkdwn = circularmod(kdwn,Nz);
                graducvx=0;
                graducvy=0;
                graducvz=0;
                /*curve to gridpoint down distance*/
                xd = (testx - x(idwn,griddata))/h;
                yd = (testy - y(jdwn,griddata))/h;
                zd = (testz - z(kdwn,griddata))/h;
                for(m=0;m<8;m++)  //linear interpolation from 8 nearest neighbours
                {
                    /* Work out increments*/
                    iinc = m%2;
                    jinc = (m/2)%2;
                    kinc = (m/4)%2;
                    /*Loop over nearest points*/
                    i = gridinc(modidwn, iinc, Nx,0);
                    j = gridinc(modjdwn, jinc, Ny,1);
                    k = gridinc(modkdwn,kinc, Nz,2);
                    prefactor = (1-iinc + pow(-1,1+iinc)*xd)*(1-jinc + pow(-1,1+jinc)*yd)*(1-kinc + pow(-1,1+kinc)*zd);
                    /*interpolate gradients of |grad u x grad v|*/
                    graducvx += prefactor*(sqrt(ucvx[pt(gridinc(i,1,Nx,0),j,k,griddata)]*ucvx[pt(gridinc(i,1,Nx,0),j,k,griddata)] + ucvy[pt(gridinc(i,1,Nx,0),j,k,griddata)]*ucvy[pt(gridinc(i,1,Nx,0),j,k,griddata)] + ucvz[pt(gridinc(i,1,Nx,0),j,k,griddata)]*ucvz[pt(gridinc(i,1,Nx,0),j,k,griddata)]) - sqrt(ucvx[pt(gridinc(i,-1,Nx,0),j,k,griddata)]*ucvx[pt(gridinc(i,-1,Nx,0),j,k,griddata)] + ucvy[pt(gridinc(i,-1,Nx,0),j,k,griddata)]*ucvy[pt(gridinc(i,-1,Nx,0),j,k,griddata)] + ucvz[pt(gridinc(i,-1,Nx,0),j,k,griddata)]*ucvz[pt(gridinc(i,-1,Nx,0),j,k,griddata)]))/(2*h);
                    graducvy += prefactor*(sqrt(ucvx[pt(i,gridinc(j,1,Ny,1),k,griddata)]*ucvx[pt(i,gridinc(j,1,Ny,1),k,griddata)] + ucvy[pt(i,gridinc(j,1,Ny,1),k,griddata)]*ucvy[pt(i,gridinc(j,1,Ny,1),k,griddata)] + ucvz[pt(i,gridinc(j,1,Ny,1),k,griddata)]*ucvz[pt(i,gridinc(j,1,Ny,1),k,griddata)]) - sqrt(ucvx[pt(i,gridinc(j,-1,Ny,1),k,griddata)]*ucvx[pt(i,gridinc(j,-1,Ny,1),k,griddata)] + ucvy[pt(i,gridinc(j,-1,Ny,1),k,griddata)]*ucvy[pt(i,gridinc(j,-1,Ny,1),k,griddata)] + ucvz[pt(i,gridinc(j,-1,Ny,1),k,griddata)]*ucvz[pt(i,gridinc(j,-1,Ny,1),k,griddata)]))/(2*h);
                    graducvz += prefactor*(sqrt(ucvx[pt(i,j,gridinc(k,1,Nz,2),griddata)]*ucvx[pt(i,j,gridinc(k,1,Nz,2),griddata)] + ucvy[pt(i,j,gridinc(k,1,Nz,2),griddata)]*ucvy[pt(i,j,gridinc(k,1,Nz,2),griddata)] + ucvz[pt(i,j,gridinc(k,1,Nz,2),griddata)]*ucvz[pt(i,j,gridinc(k,1,Nz,2),griddata)]) - sqrt(ucvx[pt(i,j,gridinc(k,-1,Nz,2),griddata)]*ucvx[pt(i,j,gridinc(k,-1,Nz,2),griddata)] + ucvy[pt(i,j,gridinc(k,-1,Nz,2),griddata)]*ucvy[pt(i,j,gridinc(k,-1,Nz,2),griddata)] + ucvz[pt(i,j,gridinc(k,-1,Nz,2),griddata)]*ucvz[pt(i,j,gridinc(k,-1,Nz,2),griddata)]))/(2*h);

                }
                knotcurves[c].knotcurve.push_back(knotpoint());
                // one of the vectors in the plane we wish to perfrom our minimisation in
                fx = (graducvx - (graducvx*ucvxs + graducvy*ucvys + graducvz*ucvzs)*ucvxs);
                fy = (graducvy - (graducvx*ucvxs + graducvy*ucvys + graducvz*ucvzs)*ucvys);
                fz = (graducvz - (graducvx*ucvxs + graducvy*ucvys + graducvz*ucvzs)*ucvzs);
                norm = sqrt(fx*fx + fy*fy + fz*fz);
                fx = fx/norm;
                fy = fy/norm;
                fz = fz/norm;

                // okay we have our direction to perfrom the line minimisation in
                // the point
                gsl_vector* v = gsl_vector_alloc (3);
                gsl_vector_set (v, 0, testx);
                gsl_vector_set (v, 1, testy);
                gsl_vector_set (v, 2, testz);
                // one vector in the plane we with to minimize in
                gsl_vector* f = gsl_vector_alloc (3);
                gsl_vector_set (f, 0, fx);
                gsl_vector_set (f, 1, fy);
                gsl_vector_set (f, 2, fz);
                // the ucv vector
                gsl_vector* ucv = gsl_vector_alloc (3);
                gsl_vector_set (ucv, 0, ucvxs);
                gsl_vector_set (ucv, 1, ucvys);
                gsl_vector_set (ucv, 2, ucvzs);
                // take a cross product to get the other vector in the plane
                gsl_vector* b = gsl_vector_alloc (3);
                cross_product(f,ucv,b);
                // initial conditions
                gsl_vector* minimum = gsl_vector_alloc (2);
                gsl_vector_set (minimum, 0, 0);
                gsl_vector_set (minimum, 1, 0);
                struct parameters params; struct parameters* pparams = &params;
                pparams->ucvmag=&interpolateducvmag;
                pparams->v = v; pparams->f = f;pparams->b=b;
                pparams->mygriddata = griddata;
                // some initial values
                gsl_multimin_function F;
                F.n=2;
                F.f = &my_f;
                F.params = (void*) pparams;
                gsl_vector* stepsize = gsl_vector_alloc (2);
                gsl_vector_set (stepsize, 0, lambda/(8*M_PI));
                gsl_vector_set (stepsize, 1, lambda/(8*M_PI));
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
                    status = gsl_multimin_test_size (minimizersize, 1e-2);

                }
                while (status == GSL_CONTINUE && iter < 5000);


                gsl_vector_scale(f,gsl_vector_get(minimizerstate->x, 0));
                gsl_vector_scale(b,gsl_vector_get(minimizerstate->x, 1));
                gsl_vector_add(f,b);
                gsl_vector_add(v,f);
                knotcurves[c].knotcurve[s].xcoord = gsl_vector_get(v, 0);
                knotcurves[c].knotcurve[s].ycoord= gsl_vector_get(v, 1);
                knotcurves[c].knotcurve[s].zcoord= gsl_vector_get(v, 2);

                gsl_vector_free(v);
                gsl_vector_free(f);
                gsl_vector_free(b);
                gsl_vector_free(ucv);
                gsl_vector_free(stepsize);

                xdiff = knotcurves[c].knotcurve[0].xcoord - knotcurves[c].knotcurve[s].xcoord;     //distance from start/end point
                ydiff = knotcurves[c].knotcurve[0].ycoord - knotcurves[c].knotcurve[s].ycoord;
                zdiff = knotcurves[c].knotcurve[0].zcoord - knotcurves[c].knotcurve[s].zcoord;

                if( (boundaryhits==0 && sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff) <h  && s > 10 ) || boundaryhits>1 ||s>5000) finish = true;

                // okay, we just added a point in position s in the vector
                // if we have a few points in the vector, discard the first few and restart the whole thing - burn it in
                int newstartingposition =20;
                if(s==newstartingposition && burnin && (boundaryhits==0))
                {
                    knotcurves[c].knotcurve.erase(knotcurves[c].knotcurve.begin(),knotcurves[c].knotcurve.begin()+newstartingposition);
                    s =0;
                    burnin =false;
                }

                s++;
            }
            int NP = knotcurves[c].knotcurve.size();  //store number of points in knot curve

            // construct a tube around the knot, to use as an excluded region if we searching for multiple components.
            double radius = 3;
            int numiterations = (int)(radius/griddata.h);
            for(int s=0; s<NP; s++)
            {
                int icentral = (int) ((knotcurves[c].knotcurve[s].xcoord/h) - 0.5 + Nx/2.0);
                int jcentral = (int) ((knotcurves[c].knotcurve[s].ycoord/h) - 0.5 + Ny/2.0);
                int kcentral = (int) ((knotcurves[c].knotcurve[s].zcoord/h) - 0.5 + Nz/2.0);
                // construct a ball of radius "radius" around each point in the knotcurve object. we circumscribe it in a cube which is then looped over
                for(int i =-numiterations;i<=numiterations;i++)
               {
                    for(int j=-numiterations ;j<=numiterations;j++)
                    {
                        for(int k =-numiterations;k<=numiterations;k++)
                        {
                            int modi = circularmod(i+icentral,Nx);
                            int modj = circularmod(j+jcentral,Ny);
                            int modk = circularmod(k+kcentral,Nz);
                            int n = pt(modi,modj,modk,griddata);

                            double dxsq = (x(i+icentral,griddata)-x(icentral,griddata))*(x(i+icentral,griddata)-x(icentral,griddata));
                            double dysq = (y(j+jcentral,griddata)-y(jcentral,griddata))*(y(j+jcentral,griddata)-y(jcentral,griddata));
                            double dzsq = (z(k+kcentral,griddata)-z(kcentral,griddata))*(z(k+kcentral,griddata)-z(kcentral,griddata));

                            double r = sqrt(dxsq + dysq + dzsq);

                            if(r < radius )
                            {
                                marked[n]=1;
                            }
                        }
                    }
                }
            }

            // now comes a lot of curve analysis. but for now Im not going to do this on the boundary curves. 
            if(boundaryhits==0)
            {

                /*******Vertex averaging*********/

                double totlength, dl, dx,dy,dz;
                for(i=0;i<3;i++)   //repeat a couple of times because of end point
                {
                    totlength=0;
                    for(s=0; s<NP; s++)   //Work out total length of curve
                    {
                        dx = knotcurves[c].knotcurve[incp(s,1,NP)].xcoord - knotcurves[c].knotcurve[s].xcoord;
                        dy = knotcurves[c].knotcurve[incp(s,1,NP)].ycoord - knotcurves[c].knotcurve[s].ycoord;
                        dz = knotcurves[c].knotcurve[incp(s,1,NP)].zcoord - knotcurves[c].knotcurve[s].zcoord;
                        totlength += sqrt(dx*dx + dy*dy + dz*dz);
                    }
                    dl = totlength/NP;
                    for(s=0; s<NP; s++)    //Move points to have spacing dl
                    {
                        dx = knotcurves[c].knotcurve[incp(s,1,NP)].xcoord - knotcurves[c].knotcurve[s].xcoord;
                        dy = knotcurves[c].knotcurve[incp(s,1,NP)].ycoord - knotcurves[c].knotcurve[s].ycoord;
                        dz = knotcurves[c].knotcurve[incp(s,1,NP)].zcoord - knotcurves[c].knotcurve[s].zcoord;
                        double norm = sqrt(dx*dx + dy*dy + dz*dz);
                        knotcurves[c].knotcurve[incp(s,1,NP)].xcoord = knotcurves[c].knotcurve[s].xcoord + dl*dx/norm;
                        knotcurves[c].knotcurve[incp(s,1,NP)].ycoord = knotcurves[c].knotcurve[s].ycoord + dl*dy/norm;
                        knotcurves[c].knotcurve[incp(s,1,NP)].zcoord = knotcurves[c].knotcurve[s].zcoord + dl*dz/norm;
                    }
                }

                /*************Curve Smoothing*******************/
                vector<double> coord(NP);
                gsl_fft_real_wavetable * real;
                gsl_fft_halfcomplex_wavetable * hc;
                gsl_fft_real_workspace * work;
                work = gsl_fft_real_workspace_alloc (NP);
                real = gsl_fft_real_wavetable_alloc (NP);
                hc = gsl_fft_halfcomplex_wavetable_alloc (NP);
                for(j=1; j<4; j++)
                {
                    switch(j)
                    {
                        case 1 :
                            for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].xcoord ; break;
                        case 2 :
                            for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].ycoord ; break;
                        case 3 :
                            for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].zcoord ; break;
                    }
                    double* data = coord.data();
                    // take the fft
                    gsl_fft_real_transform (data, 1, NP, real, work);
                    // 21/11/2016: make our low pass filter. To apply our filter. we should sample frequencies fn = n/Delta N , n = -N/2 ... N/2
                    // this is discretizing the nyquist interval, with extreme frequency ~1/2Delta.
                    // to cut out the frequencies of grid fluctuation size and larger we need a lengthscale Delta to
                    // plug in above. im doing a rough length calc below, this might be overkill.
                    // at the moment its just a hard filter, we can choose others though.
                    // compute a rough length to set scale
                    double filter;
                    const double cutoff = 2*M_PI*(totlength/(6*lambda));
                    for (i = 0; i < NP; ++i)
                    {
                        filter = 1/sqrt(1+pow((i/cutoff),8));
                        data[i] *= filter;
                    };
                    // transform back
                    gsl_fft_halfcomplex_inverse (data, 1, NP, hc, work);
                    switch(j)
                    {
                        case 1 :
                            for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].xcoord = coord[i] ; break;
                        case 2 :
                            for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].ycoord = coord[i] ; break;
                        case 3 :
                            for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].zcoord = coord[i] ; break;
                    }
                }



                /******************Interpolate direction of grad u for twist calc*******/
                /**Find nearest gridpoint**/
                double dxu, dyu, dzu, dxup, dyup, dzup;
                for(s=0; s<NP; s++)
                {
                    idwn = (int) ((knotcurves[c].knotcurve[s].xcoord/h) - 0.5 + Nx/2.0);
                    jdwn = (int) ((knotcurves[c].knotcurve[s].ycoord/h) - 0.5 + Ny/2.0);
                    kdwn = (int) ((knotcurves[c].knotcurve[s].zcoord/h) - 0.5 + Nz/2.0);
                    modidwn = circularmod(idwn,Nx);
                    modjdwn = circularmod(jdwn,Ny);
                    modkdwn = circularmod(kdwn,Nz);
                    if((BoundaryType==ALLREFLECTING) && (idwn<0 || jdwn<0 || kdwn<0 || idwn > Nx-1 || jdwn > Ny-1 || kdwn > Nz-1)) break;
                    if((BoundaryType==ZPERIODIC) && (idwn<0 || jdwn<0 || idwn > Nx-1 || jdwn > Ny-1 )) break;
                    dxu=0;
                    dyu=0;
                    dzu=0;
                    /*curve to gridpoint down distance*/
                    xd = (knotcurves[c].knotcurve[s].xcoord - x(idwn,griddata))/h;
                    yd = (knotcurves[c].knotcurve[s].ycoord - y(jdwn,griddata))/h;
                    zd = (knotcurves[c].knotcurve[s].zcoord - z(kdwn,griddata))/h;
                    for(m=0;m<8;m++)  //linear interpolation of 8 NNs
                    {
                        /* Work out increments*/
                        iinc = m%2;
                        jinc = (m/2)%2;
                        kinc = (m/4)%2;
                        /*Loop over nearest points*/
                        i = gridinc(modidwn, iinc, Nx,0);
                        j = gridinc(modjdwn, jinc, Ny,1);
                        k = gridinc(modkdwn,kinc, Nz,2);
                        prefactor = (1-iinc + pow(-1,1+iinc)*xd)*(1-jinc + pow(-1,1+jinc)*yd)*(1-kinc + pow(-1,1+kinc)*zd);   //terms of the form (1-xd)(1-yd)zd etc. (interpolation coefficient)
                        /*interpolate grad u over nearest points*/
                        dxu += prefactor*0.5*(u[pt(gridinc(i,1,Nx,0),j,k,griddata)] -  u[pt(gridinc(i,-1,Nx,0),j,k,griddata)])/h;  //central diff
                        dyu += prefactor*0.5*(u[pt(i,gridinc(j,1,Ny,1),k,griddata)] -  u[pt(i,gridinc(j,-1,Ny,1),k,griddata)])/h;
                        dzu += prefactor*0.5*(u[pt(i,j,gridinc(k,1,Nz,2),griddata)] -  u[pt(i,j,gridinc(k,-1,Nz,2),griddata)])/h;
                    }
                    //project du onto perp of tangent direction first
                    dx = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].xcoord - knotcurves[c].knotcurve[incp(s,-1,NP)].xcoord);   //central diff as a is defined on the points
                    dy = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].ycoord - knotcurves[c].knotcurve[incp(s,-1,NP)].ycoord);
                    dz = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].zcoord - knotcurves[c].knotcurve[incp(s,-1,NP)].zcoord);
                    dxup = dxu - (dxu*dx + dyu*dy + dzu*dz)*dx/(dx*dx+dy*dy+dz*dz);               //Grad u_j * (delta_ij - t_i t_j)
                    dyup = dyu - (dxu*dx + dyu*dy + dzu*dz)*dy/(dx*dx+dy*dy+dz*dz);
                    dzup = dzu - (dxu*dx + dyu*dy + dzu*dz)*dz/(dx*dx+dy*dy+dz*dz);
                    /*Vector a is the normalised gradient of u, should point in direction of max u perp to t*/
                    double norm = sqrt(dxup*dxup+dyup*dyup+dzup*dzup);
                    knotcurves[c].knotcurve[s].ax = dxup/norm;
                    knotcurves[c].knotcurve[s].ay = dyup/norm;
                    knotcurves[c].knotcurve[s].az = dzup/norm;
                }

                for(j=1; j<4; j++)
                {
                    switch(j)
                    {
                        case 1 :
                            for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].ax ; break;
                        case 2 :
                            for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].ay ; break;
                        case 3 :
                            for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].az ; break;
                    }
                    double* data = coord.data();
                    // take the fft
                    gsl_fft_real_transform (data, 1, NP, real, work);
                    // 21/11/2016: make our low pass filter. To apply our filter. we should sample frequencies fn = n/Delta N , n = -N/2 ... N/2
                    // this is discretizing the nyquist interval, with extreme frequency ~1/2Delta.
                    // to cut out the frequencies of grid fluctuation size and larger we need a lengthscale Delta to
                    // plug in above. im doing a rough length calc below, this might be overkill.
                    // at the moment its just a hard filter, we can choose others though.
                    // compute a rough length to set scale
                    double filter;
                    const double cutoff = 2*M_PI*(totlength/(6*lambda));
                    for (i = 0; i < NP; ++i)
                    {
                        filter = 1/sqrt(1+pow((i/cutoff),8));
                        data[i] *= filter;
                    };
                    // transform back
                    gsl_fft_halfcomplex_inverse (data, 1, NP, hc, work);
                    switch(j)
                    {
                        case 1 :
                            for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].ax= coord[i] ; break;
                        case 2 :
                            for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].ay= coord[i] ; break;
                        case 3 :
                            for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].az = coord[i] ; break;
                    }
                }
                gsl_fft_real_wavetable_free (real);
                gsl_fft_halfcomplex_wavetable_free (hc);
                gsl_fft_real_workspace_free (work);


                // CURVE GEOMETRY - get curvatures, torsions, frennet serret frame


                NP = knotcurves[c].knotcurve.size();
                for(s=0; s<NP; s++)
                {
                    // forward difference on the tangents
                    double dx = (knotcurves[c].knotcurve[incp(s,1,NP)].xcoord - knotcurves[c].knotcurve[incp(s,0,NP)].xcoord);
                    double dy = (knotcurves[c].knotcurve[incp(s,1,NP)].ycoord - knotcurves[c].knotcurve[incp(s,0,NP)].ycoord);
                    double dz = (knotcurves[c].knotcurve[incp(s,1,NP)].zcoord - knotcurves[c].knotcurve[incp(s,0,NP)].zcoord);
                    double deltas = sqrt(dx*dx+dy*dy+dz*dz);
                    knotcurves[c].knotcurve[s].tx = dx/(deltas);
                    knotcurves[c].knotcurve[s].ty = dy/(deltas);
                    knotcurves[c].knotcurve[s].tz = dz/(deltas);
                    knotcurves[c].knotcurve[s].length = deltas;
                    knotcurves[c].length +=deltas;
                }
                for(s=0; s<NP; s++)
                {
                    // backwards diff for the normals, amounting to a central diff overall
                    double nx = 2.0*(knotcurves[c].knotcurve[s].tx-knotcurves[c].knotcurve[incp(s,-1,NP)].tx)/(knotcurves[c].knotcurve[s].length+knotcurves[c].knotcurve[incp(s,-1,NP)].length);
                    double ny = 2.0*(knotcurves[c].knotcurve[s].ty-knotcurves[c].knotcurve[incp(s,-1,NP)].ty)/(knotcurves[c].knotcurve[s].length+knotcurves[c].knotcurve[incp(s,-1,NP)].length);
                    double nz = 2.0*(knotcurves[c].knotcurve[s].tz-knotcurves[c].knotcurve[incp(s,-1,NP)].tz)/(knotcurves[c].knotcurve[s].length+knotcurves[c].knotcurve[incp(s,-1,NP)].length);
                    double curvature = sqrt(nx*nx+ny*ny+nz*nz);
                    nx /=curvature;
                    ny /=curvature;
                    nz /=curvature;
                    double tx = knotcurves[c].knotcurve[s].tx ;
                    double ty =  knotcurves[c].knotcurve[s].ty ;
                    double tz = knotcurves[c].knotcurve[s].tz ;
                    double bx = ty*nz - tz*ny;
                    double by = tz*nx - tx*nz;
                    double bz = tx*ny - ty*nx;
                    knotcurves[c].knotcurve[s].nx = nx ;
                    knotcurves[c].knotcurve[s].ny = ny ;
                    knotcurves[c].knotcurve[s].nz = nz ;
                    knotcurves[c].knotcurve[s].bx = bx ;
                    knotcurves[c].knotcurve[s].by = by ;
                    knotcurves[c].knotcurve[s].bz = bz ;
                    knotcurves[c].knotcurve[s].curvature = curvature ;
                }
                // torsions with a central difference
                for(s=0; s<NP; s++)
                {
                    double bx = knotcurves[c].knotcurve[s].bx;
                    double by =  knotcurves[c].knotcurve[s].by;
                    double bz = knotcurves[c].knotcurve[s].bz;

                    double dnxds = 2.0*(knotcurves[c].knotcurve[incp(s,1,NP)].nx-knotcurves[c].knotcurve[incp(s,-1,NP)].nx)/(knotcurves[c].knotcurve[incp(s,1,NP)].length+knotcurves[c].knotcurve[incp(s,-1,NP)].length);
                    double dnyds = 2.0*(knotcurves[c].knotcurve[incp(s,1,NP)].ny-knotcurves[c].knotcurve[incp(s,-1,NP)].ny)/(knotcurves[c].knotcurve[incp(s,1,NP)].length+knotcurves[c].knotcurve[incp(s,-1,NP)].length);
                    double dnzds = 2.0*(knotcurves[c].knotcurve[incp(s,1,NP)].nz-knotcurves[c].knotcurve[incp(s,-1,NP)].nz)/(knotcurves[c].knotcurve[incp(s,1,NP)].length+knotcurves[c].knotcurve[incp(s,-1,NP)].length);

                    double torsion = bx*dnxds+by*dnyds+bz*dnzds;
                    knotcurves[c].knotcurve[s].torsion = torsion ;
                }


                // RIBBON TWIST AND WRITHE

                for(s=0; s<NP; s++)
                {

                    // twist of this segment
                    double ds = knotcurves[c].knotcurve[s].length;
                    double dxds = knotcurves[c].knotcurve[s].tx;
                    double dyds = knotcurves[c].knotcurve[s].ty;
                    double dzds = knotcurves[c].knotcurve[s].tz;
                    double bx = (knotcurves[c].knotcurve[incp(s,1,NP)].ax - knotcurves[c].knotcurve[s].ax)/ds;
                    double by = (knotcurves[c].knotcurve[incp(s,1,NP)].ay - knotcurves[c].knotcurve[s].ay)/ds;
                    double bz = (knotcurves[c].knotcurve[incp(s,1,NP)].az - knotcurves[c].knotcurve[s].az)/ds;
                    knotcurves[c].knotcurve[s].twist = (dxds*(knotcurves[c].knotcurve[s].ay*bz - knotcurves[c].knotcurve[s].az*by) + dyds*(knotcurves[c].knotcurve[s].az*bx - knotcurves[c].knotcurve[s].ax*bz) + dzds*(knotcurves[c].knotcurve[s].ax*by - knotcurves[c].knotcurve[s].ay*bx))/(2*M_PI*sqrt(dxds*dxds + dyds*dyds + dzds*dzds));

                    // "writhe" of this segment. writhe is nonlocal, this is the thing in the integrand over s
                    knotcurves[c].knotcurve[s].writhe = 0;
                    for(m=0; m<NP; m++)
                    {
                        if(s != m)
                        {
                            xdiff = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].xcoord + knotcurves[c].knotcurve[s].xcoord - knotcurves[c].knotcurve[incp(m,1,NP)].xcoord - knotcurves[c].knotcurve[m].xcoord);   //interpolate, consistent with fwd diff
                            ydiff = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].ycoord + knotcurves[c].knotcurve[s].ycoord - knotcurves[c].knotcurve[incp(m,1,NP)].ycoord - knotcurves[c].knotcurve[m].ycoord);
                            zdiff = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].zcoord + knotcurves[c].knotcurve[s].zcoord - knotcurves[c].knotcurve[incp(m,1,NP)].zcoord - knotcurves[c].knotcurve[m].zcoord);
                            double dxdm = (knotcurves[c].knotcurve[incp(m,1,NP)].xcoord - knotcurves[c].knotcurve[m].xcoord)/(ds);
                            double dydm = (knotcurves[c].knotcurve[incp(m,1,NP)].ycoord - knotcurves[c].knotcurve[m].ycoord)/(ds);
                            double dzdm = (knotcurves[c].knotcurve[incp(m,1,NP)].zcoord - knotcurves[c].knotcurve[m].zcoord)/(ds);
                            knotcurves[c].knotcurve[s].writhe += ds*(xdiff*(dyds*dzdm - dzds*dydm) + ydiff*(dzds*dxdm - dxds*dzdm) + zdiff*(dxds*dydm - dyds*dxdm))/(4*M_PI*(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)*sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff));
                        }
                    }

                    //Add on writhe, twist
                    knotcurves[c].writhe += knotcurves[c].knotcurve[s].writhe*ds;
                    knotcurves[c].twist  += knotcurves[c].knotcurve[s].twist*ds;
                    // while we are computing the global quantites, get the average position too
                    knotcurves[c].xavgpos += knotcurves[c].knotcurve[s].xcoord/NP;
                    knotcurves[c].yavgpos += knotcurves[c].knotcurve[s].ycoord/NP;
                    knotcurves[c].zavgpos += knotcurves[c].knotcurve[s].zcoord/NP;
                }


                // the ghost grid has been useful for painlessly computing all the above quantities, without worrying about the periodic bc's
                // but for storage and display, we should put it all in the box

                // (1) construct the proper periodic co-ordinates from our ghost grid
                double xupperlim = x(griddata.Nx -1 ,griddata);
                double xlowerlim = x(0,griddata);
                double deltax = griddata.Nx * h;
                double yupperlim = y(griddata.Ny -1,griddata);
                double ylowerlim = y(0,griddata);
                double deltay = griddata.Ny * h;
                double zupperlim = z(griddata.Nz -1 ,griddata);
                double zlowermin = z(0,griddata);
                double deltaz = griddata.Nz * h;
                for(s=0; s<NP; s++)
                {
                    knotcurves[c].knotcurve[s].modxcoord = knotcurves[c].knotcurve[s].xcoord;
                    knotcurves[c].knotcurve[s].modycoord = knotcurves[c].knotcurve[s].ycoord;
                    knotcurves[c].knotcurve[s].modzcoord = knotcurves[c].knotcurve[s].zcoord;
                    if(knotcurves[c].knotcurve[s].xcoord > xupperlim) {
                        knotcurves[c].knotcurve[s].modxcoord = knotcurves[c].knotcurve[s].xcoord-deltax;
                    } ;
                    if(knotcurves[c].knotcurve[s].xcoord < xlowerlim) {
                        knotcurves[c].knotcurve[s].modxcoord = knotcurves[c].knotcurve[s].xcoord+deltax;
                    };
                    if(knotcurves[c].knotcurve[s].ycoord > yupperlim) {
                        knotcurves[c].knotcurve[s].modycoord = knotcurves[c].knotcurve[s].ycoord-deltay;
                    };
                    if(knotcurves[c].knotcurve[s].ycoord < ylowerlim) {
                        knotcurves[c].knotcurve[s].modycoord = knotcurves[c].knotcurve[s].ycoord+deltay;
                    };
                    if(knotcurves[c].knotcurve[s].zcoord > zupperlim) {
                        knotcurves[c].knotcurve[s].modzcoord = knotcurves[c].knotcurve[s].zcoord-deltaz;
                    };
                    if(knotcurves[c].knotcurve[s].zcoord < zlowermin)
                    {
                        knotcurves[c].knotcurve[s].modzcoord = knotcurves[c].knotcurve[s].zcoord+deltaz;
                    };
                }

                // (2) standardise the knot such that the top right corner of the bounding box lies in the "actual" grid. This bounding box point may only lie
                // off grid in the +ve x y z direction.
                double xmax=knotcurves[c].knotcurve[0].xcoord;
                double ymax=knotcurves[c].knotcurve[0].ycoord;
                double zmax=knotcurves[c].knotcurve[0].zcoord;
                for(s=0; s<NP; s++)
                {
                    if(knotcurves[c].knotcurve[s].xcoord>xmax)
                    {
                        xmax = knotcurves[c].knotcurve[s].xcoord;
                    }
                    if(knotcurves[c].knotcurve[s].ycoord>ymax)
                    {
                        ymax = knotcurves[c].knotcurve[s].ycoord;
                    }
                    if(knotcurves[c].knotcurve[s].zcoord>zmax)
                    {
                        zmax = knotcurves[c].knotcurve[s].zcoord;
                    }
                }

                // get how many lattice shifts are needed
                int xlatticeshift = (int) (round(xmax/(griddata.Nx *griddata.h)));
                int ylatticeshift = (int) (round(ymax/(griddata.Ny *griddata.h)));
                int zlatticeshift = (int) (round(zmax/(griddata.Nz *griddata.h)));
                // perform the shift

                for(int s=0; s<knotcurves[c].knotcurve.size(); s++)
                {
                    knotcurves[c].knotcurve[s].xcoord -= (double)(xlatticeshift) * (griddata.Nx *griddata.h);
                    knotcurves[c].knotcurve[s].ycoord -= (double)(ylatticeshift) * (griddata.Ny *griddata.h);
                    knotcurves[c].knotcurve[s].zcoord -= (double)(zlatticeshift) * (griddata.Nz *griddata.h);
                }
                // now we've done these shifts, we'd better move the knotcurve average position too.
                knotcurves[c].xavgpos = 0;
                knotcurves[c].yavgpos = 0;
                knotcurves[c].zavgpos = 0;
                for(int s=0; s<knotcurves[c].knotcurve.size(); s++)
                {
                    knotcurves[c].xavgpos += knotcurves[c].knotcurve[s].xcoord/NP;
                    knotcurves[c].yavgpos += knotcurves[c].knotcurve[s].ycoord/NP;
                    knotcurves[c].zavgpos += knotcurves[c].knotcurve[s].zcoord/NP;
                }
                c++;
            }
            // if we did hit a boundary, just strike the curve we made from the record. It lives on in the marked array!
            else
            {
                knotcurves.pop_back();
            }
        }
    }
    // the order of the components within the knotcurves vector is not guaranteed to remain fixed from timestep to timestep. thus, componenet 0 at one timtestep could be
    // components 1 at the next. the code needs a way of tracking which componenet is which.
    // at the moment, im doing this by fuzzily comparing summary stats on the components - at this point, the length twist and writhe.

    // these variables have the summary stats from the last timestep

    static vector<double> oldwrithe(knotcurves.size());
    static vector<double> oldtwist(knotcurves.size());
    static vector<double> oldlength(knotcurves.size());
    static vector<double> oldxavgpos(knotcurves.size());
    static vector<double> oldyavgpos(knotcurves.size());
    static vector<double> oldzavgpos(knotcurves.size());
    static bool first = true;
    static vector<int> permutation(knotcurves.size());

    if(first)
    {
        for(int i = 0; i<knotcurves.size();i++)
        {
            oldwrithe[i] = knotcurves[i].writhe;
            oldlength[i] = knotcurves[i].length;
            oldtwist[i] = knotcurves[i].twist;
            oldxavgpos[i] = knotcurves[i].xavgpos;
            oldyavgpos[i] = knotcurves[i].yavgpos;
            oldzavgpos[i] = knotcurves[i].zavgpos;
            permutation[i] = i;
        }

    }
    else
    {
        for(int i = 0; i<knotcurves.size();i++)
        {
            double minscore = INFINITY;
            for(int j = 0; j<knotcurves.size();j++)
            {
                double score = fabs((knotcurves[j].length - oldlength[i])/(Nx*h))+fabs(knotcurves[j].writhe - oldwrithe[i]);
                score += fmod(fabs((knotcurves[j].xavgpos - oldxavgpos[i])),Nx*h)/(Nx*h);
                score += fmod(fabs((knotcurves[j].yavgpos - oldyavgpos[i])),Ny*h)/(Ny*h);
                score += fmod(fabs((knotcurves[j].zavgpos - oldzavgpos[i])),Nz*h)/(Nz*h);
                if(score<minscore) {permutation[i] = j; minscore = score;}

            }
        }

        // if the permutation isn't valid, print a warning and reset it - it's better to have the curves swap then have them incorrectly both mapped to the same curve.
        bool isvalidperm = true;
        for(int i=0;i<knotcurves.size();i++)
        {
            std::vector<int>::iterator it = std::find(permutation.begin(),permutation.end(),i);
            if(it == permutation.end()){isvalidperm=false;}
        }

        if(!isvalidperm)
        {
            cout << "the permutation was not valid. resetting! \n";
            for(int i = 0; i<knotcurves.size();i++)
            {
                permutation[i] = i;
            }

        }


        // apply the permutation to the list of lengths etc. It is now "correct" in the sense that index [0] really is component 0 etc. these labellings are arbitrarlly set at the simulations start and
        // must be consistently carried forward
        for(int i = 0; i<knotcurves.size();i++)
        {
            oldwrithe[i] = knotcurves[permutation[i]].writhe;
            oldlength[i] = knotcurves[permutation[i]].length;
            oldtwist[i] = knotcurves[permutation[i]].twist;
            oldxavgpos[i] = knotcurves[permutation[i]].xavgpos;
            oldyavgpos[i] = knotcurves[permutation[i]].yavgpos;
            oldzavgpos[i] = knotcurves[permutation[i]].zavgpos;

        }
        // i can't be bothered to do this the smart way.
        vector<knotcurve> tempknotcurves(knotcurves.size());
        for(int i = 0; i<knotcurves.size();i++)
        {
            tempknotcurves[i]=knotcurves[permutation[i]];
        }
        knotcurves=tempknotcurves;
    }
    first = false;
}

void find_knot_velocity(const vector<knotcurve>& knotcurves,vector<knotcurve>& knotcurvesold,const Griddata& griddata,const double deltatime)
{
    for(int c=0;c<knotcurvesold.size();c++)
    {

        int NP = knotcurves[c].knotcurve.size();
        int NPold = knotcurvesold[c].knotcurve.size();

        for(int s = 0; s< knotcurvesold[c].knotcurve.size(); s++)
        {
            double IntersectionFraction =-1;
            std::vector<double> IntersectionPoint(3);
            std::vector<double> ClosestIntersection(3);
            double closestdistancesquare = knotcurvesold[c].length;
            for(int t = 0 ; t<knotcurves[c].knotcurve.size();t++)
            {
                int intersection = 0;
                intersection = intersect3D_SegmentPlane( knotcurves[c].knotcurve[t%NP], knotcurves[c].knotcurve[(t+1)%NP], knotcurvesold[c].knotcurve[s%NPold], knotcurvesold[c].knotcurve[(s+1)%NPold], IntersectionFraction, IntersectionPoint );
                if(intersection ==1)
                {
                    double intersectiondistancesquare = (IntersectionPoint[0] - knotcurvesold[c].knotcurve[s].xcoord )*(IntersectionPoint[0] - knotcurvesold[c].knotcurve[s].xcoord )+ (IntersectionPoint[1] - knotcurvesold[c].knotcurve[s].ycoord )*(IntersectionPoint[1] - knotcurvesold[c].knotcurve[s].ycoord )+ (IntersectionPoint[2] - knotcurvesold[c].knotcurve[s].zcoord )*(IntersectionPoint[2] - knotcurvesold[c].knotcurve[s].zcoord );
                    if(intersectiondistancesquare < closestdistancesquare)
                    {
                        closestdistancesquare = intersectiondistancesquare;
                        ClosestIntersection[0] = IntersectionPoint[0];
                        ClosestIntersection[1] = IntersectionPoint[1];
                        ClosestIntersection[2] = IntersectionPoint[2];
                    }

                }
            }
            // work out velocity and twist rate
            knotcurvesold[c].knotcurve[s].vx = (ClosestIntersection[0] - knotcurvesold[c].knotcurve[s].xcoord )/ deltatime;
            knotcurvesold[c].knotcurve[s].vy = (ClosestIntersection[1] - knotcurvesold[c].knotcurve[s].ycoord )/ deltatime;
            knotcurvesold[c].knotcurve[s].vz = (ClosestIntersection[2] - knotcurvesold[c].knotcurve[s].zcoord )/ deltatime;
            // for convenience, lets also output the decomposition into normal and binormal
            double vdotn = knotcurvesold[c].knotcurve[s].nx*knotcurvesold[c].knotcurve[s].vx+knotcurvesold[c].knotcurve[s].ny*knotcurvesold[c].knotcurve[s].vy+knotcurvesold[c].knotcurve[s].nz*knotcurvesold[c].knotcurve[s].vz;
            double vdotb = knotcurvesold[c].knotcurve[s].bx*knotcurvesold[c].knotcurve[s].vx+knotcurvesold[c].knotcurve[s].by*knotcurvesold[c].knotcurve[s].vy+knotcurvesold[c].knotcurve[s].bz*knotcurvesold[c].knotcurve[s].vz;

            knotcurvesold[c].knotcurve[s].vdotnx = vdotn * knotcurvesold[c].knotcurve[s].nx ;
            knotcurvesold[c].knotcurve[s].vdotny = vdotn * knotcurvesold[c].knotcurve[s].ny ;
            knotcurvesold[c].knotcurve[s].vdotnz = vdotn * knotcurvesold[c].knotcurve[s].nz ;
            knotcurvesold[c].knotcurve[s].vdotbx = vdotb * knotcurvesold[c].knotcurve[s].bx ;
            knotcurvesold[c].knotcurve[s].vdotby = vdotb * knotcurvesold[c].knotcurve[s].by ;
            knotcurvesold[c].knotcurve[s].vdotbz = vdotb * knotcurvesold[c].knotcurve[s].bz ;
        }
    }
}
void uv_update(vector<double>&u, vector<double>&v,  vector<double>&ku, vector<double>&kv,const Griddata& griddata)
{
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    const double h = griddata.h;
    int i,j,k,l,n,kup,kdown,iup,idown,jup,jdown;
    double D2u;
    const int arraysize = Nx*Ny*Nz;


    // some constants we will use over and over below:
    const double sixth = 1.0/6.0;
    const double ONETHIRD = 1.0/3.0;
    const double oneoverepsilon = 1.0/epsilon;
    const double oneoverhsq = 1.0/(h*h);
    // first loop. get k1, store (in testun] and testv[n], the value u[n]+h/2k1)
#pragma omp for 
    for(i=0;i<Nx;i++)
    {
        for(j=0; j<Ny; j++)
        {
            for(k=0; k<Nz; k++)   //Central difference
            {
                n = pt(i,j,k,griddata);
                kup = gridinc(k,1,Nz,2);
                kdown = gridinc(k,-1,Nz,2);
                D2u = oneoverhsq*(u[pt(gridinc(i,1,Nx,0),j,k,griddata)] + u[pt(gridinc(i,-1,Nx,0),j,k,griddata)] + u[pt(i,gridinc(j,1,Ny,1),k,griddata)] + u[pt(i,gridinc(j,-1,Ny,1),k,griddata)] + u[pt(i,j,kup,griddata)] + u[pt(i,j,kdown,griddata)] - 6.0*u[n]);
                ku[n] = oneoverepsilon*(u[n] - (ONETHIRD*u[n])*(u[n]*u[n]) - v[n]) + D2u;
                kv[n] = epsilon*(u[n] + beta - gam*v[n]);
            }
        }
    }
    // 2nd and 3rd loops
    double inc ;
    for(l=1;l<=3;l++)  //u and v update for each fractional time step
    {
        switch (l)
        {
            case 1:
                {
                    inc=0.5;   //add k1 to uv and add to total k
                }
                break;

            case 2:
                {
                    inc=0.5 ;   //add k1 to uv and add to total k
                }
                break;
            case 3:
                {
                    inc=1 ;   //add k1 to uv and add to total k
                }
                break;
        }
#pragma omp for 
        for(i=0;i<Nx;i++)
        {
            for(j=0; j<Ny; j++)
            {
                for(k=0; k<Nz; k++)   //Central difference
                {
                    n = pt(i,j,k,griddata);

                    iup = pt(gridinc(i,1,Nx,0),j,k,griddata);
                    idown =pt(gridinc(i,-1,Nx,0),j,k,griddata);
                    jup = pt(i,gridinc(j,1,Ny,1),k,griddata);
                    jdown =pt(i,gridinc(j,-1,Ny,1),k,griddata);
                    kup = pt(i,j,gridinc(k,1,Nz,2),griddata);
                    kdown = pt(i,j,gridinc(k,-1,Nz,2),griddata);
                    double currentu = u[n] + dtime*inc*ku[(l-1)*arraysize+n];
                    double currentv = v[n] + dtime*inc*kv[(l-1)*arraysize+n];

                    D2u = oneoverhsq*((u[iup]+dtime*inc*ku[(l-1)*arraysize+iup]) + (u[idown]+dtime*inc*ku[(l-1)*arraysize+idown]) +(u[jup]+dtime*inc*ku[(l-1)*arraysize+jup]) +(u[jdown]+dtime*inc*ku[(l-1)*arraysize+jdown]) + (u[kup]+dtime*inc*ku[(l-1)*arraysize+kup]) + (u[kdown]+dtime*inc*ku[(l-1)*arraysize+kdown])- 6.0*(currentu));


                    ku[arraysize*l+n] = oneoverepsilon*(currentu - (ONETHIRD*currentu)*(currentu*currentu) - currentv) + D2u;
                    kv[arraysize*l+n] = epsilon*(currentu + beta - gam*currentv);
                }
            }
        }
    }
#pragma omp for 
    for(n=0;n<Nx*Ny*Nz;n++)
    {

        u[n] = u[n] + dtime*sixth*(ku[n]+2*ku[arraysize+n]+2*ku[2*arraysize+n]+ku[3*arraysize+n]);
        v[n] = v[n] + dtime*sixth*(kv[n]+2*kv[arraysize+n]+2*kv[2*arraysize+n]+kv[3*arraysize+n]);
    }

}

/*************************File reading and writing*****************************/

int intersect3D_SegmentPlane( knotpoint SegmentStart, knotpoint SegmentEnd, knotpoint PlaneSegmentStart, knotpoint PlaneSegmentEnd, double& IntersectionFraction, std::vector<double>& IntersectionPoint )
{
    double ux = SegmentEnd.xcoord - SegmentStart.xcoord ;
    double uy = SegmentEnd.ycoord - SegmentStart.ycoord ;
    double uz = SegmentEnd.zcoord - SegmentStart.zcoord ;

    double wx= SegmentStart.xcoord - PlaneSegmentStart.xcoord ;
    double wy = SegmentStart.ycoord - PlaneSegmentStart.ycoord ;
    double wz = SegmentStart.zcoord - PlaneSegmentStart.zcoord ;

    double nx= PlaneSegmentEnd.xcoord  - PlaneSegmentStart.xcoord ;
    double ny = PlaneSegmentEnd.ycoord  - PlaneSegmentStart.ycoord ;
    double nz = PlaneSegmentEnd.zcoord  - PlaneSegmentStart.zcoord ;

    double D = nx*ux+ ny*uy + nz*uz;
    double N = - (nx*wx+ ny*wy + nz*wz);

    if (fabs(D) < 0.01)
    {           // segment is parallel to plane
        if (N == 0)                      // segment lies in plane
            return 2;
        else
            return 0;                    // no intersection
    }

    double sI = N / D;
    if (sI < 0 || sI > 1)
        return 0;                        // no intersection


    IntersectionFraction = sI;
    IntersectionPoint[0] = SegmentStart.xcoord + sI * ux;
    IntersectionPoint[1] = SegmentStart.ycoord + sI * uy;
    IntersectionPoint[2] = SegmentStart.zcoord + sI * uz;
    return 1;
}

double my_f(const gsl_vector* minimum, void* params)
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

    double value = -1*((*interpolateducvmag)(px,py,pz));
    return value;
}
void cross_product(const gsl_vector *u, const gsl_vector *v, gsl_vector *product)
{
    double p1 = gsl_vector_get(u, 1)*gsl_vector_get(v, 2)
        - gsl_vector_get(u, 2)*gsl_vector_get(v, 1);

    double p2 = gsl_vector_get(u, 2)*gsl_vector_get(v, 0)
        - gsl_vector_get(u, 0)*gsl_vector_get(v, 2);

    double p3 = gsl_vector_get(u, 0)*gsl_vector_get(v, 1)
        - gsl_vector_get(u, 1)*gsl_vector_get(v, 0);

    gsl_vector_set(product, 0, p1);
    gsl_vector_set(product, 1, p2);
    gsl_vector_set(product, 2, p3);
}
void rotatedisplace(double& xcoord, double& ycoord, double& zcoord, const double theta, const double ux,const double uy,const double uz)
{

    double xprime = ( cos(theta) + (ux*ux)*(1-cos(theta)) )*xcoord + ( ux*uy*(1-cos(theta)) - uz*sin(theta) )*ycoord + ( ux*uz*(1-cos(theta)) + uy*sin(theta) )*zcoord;
    double yprime =(uy*ux*(1-cos(theta)) + uz*sin(theta) )*xcoord + ( cos(theta) + (uy*uy)*(1-cos(theta)) )*ycoord + ( uy*uz*(1-cos(theta)) - ux*sin(theta)  )*zcoord;
    double zprime = (uz*ux*(1-cos(theta)) - uy*sin(theta) )*xcoord + ( uz*uy*(1-cos(theta)) + ux*sin(theta)  )*ycoord + ( cos(theta) + (uz*uz)*(1-cos(theta)) )*zcoord;

    xcoord = xprime;
    ycoord = yprime;
    zcoord = zprime;

}
int circularmod(int i, int N)    // mod i by N in a cirucler fashion, ie wrapping around both in the +ve and -ve directions
{
    if(i<0) return (N - ((-i)%N))%N;
    else return i%N;
}
// inlined functions for incrementing things respecting boundaries
int incp(int i, int p, int N)    //increment i with p for periodic boundary
{
    if(i+p<0) return (N+i+p);
    else return ((i+p)%N);
}

int incw(int i, int p, int N)    //increment with reflecting boundary between -1 and 0 and N-1 and N
{
    if(i+p<0) return (-(i+p+1));
    if(i+p>N-1) return (2*N-(i+p+1));
    return (i+p);
}

int incabsorb(int i, int p, int N)    //increment with reflecting boundary between -1 and 0 and N-1 and N
{
    if(i+p<0) return (0);
    if(i+p>N-1) return (N-1);
    return (i+p);
}
// this function is specifically designed to incremenet, in the direction specified, respecting the boundary conditions, which are global enums
int gridinc(int i, int p, int N, int direction )    //increment with reflecting boundary between -1 and 0 and N-1 and N
{

    if(BoundaryType == ALLREFLECTING)
    {
        return incw(i,p,N);
    }

    if(BoundaryType == ALLPERIODIC)
    {
        return incp(i,p,N);
    }

    if(BoundaryType == ZPERIODIC)
    {
        if(direction ==2) return incp(i,p,N);
        else return incw(i,p,N);
    }
    return 0;
}
double x(int i,const Griddata& griddata)
{
    return (i+0.5-griddata.Nx/2.0)*griddata.h;
}
double y(int i,const Griddata& griddata)
{
    return (i+0.5-griddata.Ny/2.0)*griddata.h;
}
double z(int i,const Griddata& griddata)
{
    return (i+0.5-griddata.Nz/2.0)*griddata.h;
}
int pt( int i,  int j,  int k,const Griddata& griddata)       //convert i,j,k to single index
{
    return (i*griddata.Ny*griddata.Nz+j*griddata.Nz+k);
}

int coordstopt(double x, double y, double z, Griddata&griddata)
{
    double h = griddata.h;

    double doublei = (x/h) + (griddata.Nx/2.0)-0.5;
    // now, C++ truncates towards 0. We have a positive number here, which we want to round to the nearest integer. We do this by adding 0.5 and casting
    int i = (int)(doublei + 0.5);

    double doublej = (y/h) + (griddata.Ny/2.0)-0.5;
    int j = (int)(doublej + 0.5);

    double doublek = (z/h) + (griddata.Nz/2.0)-0.5;
    int k = (int)(doublek + 0.5);

    int n = pt(i,j,k,griddata);
    return n;

}

int sign(int i)
{
    if(i==0) return 0;
    else return i/abs(i);
}
