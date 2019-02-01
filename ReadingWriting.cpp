#include "ReadingWriting.h"
#include <string.h>

int file_read_ASCII( double *nx,double *ny, double *nz, double *px, double *py,double*pz)
{
    string temp,buff;
    stringstream ss;
    ifstream fin (director_filename.c_str());
    int i,j,k,n;

    for(i=0;i<10;i++)
    {
        if(fin.good())
        {
            if(getline(fin,buff)) temp = buff;
        }
        else
        {
            cout << "Something went wrong!\n";
            return 1;
        }
    }

    for(k=0; k<Lz; k++)
    {
        for(j=0; j<Ly; j++)
        {
            for(i=0; i<Lx; i++)
            {
                n=pt(i,j,k);
                ss.clear();
                ss.str("");
                if(fin.good())
                {
                    if(getline(fin,buff))
                    {
                        ss << buff;
                        ss >> nx[n] >> ny[n] >>nz[n] ;
                    }
                }
                else
                {
                    cout << "Something went wrong!\n";
                    return 1;
                }
            }
        }
    }
    fin.close();

    string temp2,buff2;
    stringstream ss2;
    ifstream fin2 (polarisation_filename.c_str());
    for(i=0;i<10;i++)
    {
        if(fin2.good())
        {
            if(getline(fin2,buff2)) temp2 = buff2;
        }
        else
        {
            cout << "Something went wrong!\n";
            return 1;
        }
    }

    for(k=0; k<Lz; k++)
    {
        for(j=0; j<Ly; j++)
        {
            for(i=0; i<Lx; i++)
            {
                n=pt(i,j,k);
                ss2.clear();
                ss2.str("");
                if(fin2.good())
                {
                    if(getline(fin2,buff2))
                    {
                        ss2 << buff2;
                        ss2 >> px[n] >> py[n] >>pz[n] ;
                    }
                }
                else
                {
                    cout << "Something went wrong!\n";
                    return 1;
                }
            }
        }
    }
    fin2.close();
    return 0;
}

int file_read( double *nx,double *ny, double *nz, double *px, double *py,double*pz)
{
    file_read_ASCII(nx,ny,nz,px,py,pz);
    return 0;
}
/**********************************************************************/
void writeVTKfiles(double n,double *nx,double *ny,double *nz,double *px,double *py,double *pz) 
{
    int j;

    char vtk_director[200],vtk_polarisation[200];
    ofstream output,output2;
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
void writeBENDfiles(double n, double* nx,double* ny,double* nz) 
{
    ofstream output;
    char vtk_BEND[200];
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

