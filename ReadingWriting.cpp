#include "ReadingWriting.h"
#include <string.h>

int file_read_ASCII( double *nx,double *ny, double *nz, double *px, double *py,double*pz)
{
    string temp,buff;
    stringstream ss;
    ifstream fin (director_filename.c_str());
    int i,j,k,n;

    for(i=0;i<9;i++)
    {
        if(fin.good())
        {
            if(getline(fin,buff)) temp = buff;
            size_t found =temp.find("DIMENSIONS");
            if(found !=std::string::npos)
            {
                temp.replace(found,10, "");
                int tempLx,tempLy,tempLz;
                ss << temp;
                ss >> tempLx >> tempLy >>tempLz ;
               if((tempLx!= Lx) ||(tempLy!= Ly)||(tempLz!= Lz))
               {
                    cout << "Box dimensions dont match\n";
                    return 1;
               }
            }
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
    for(i=0;i<9;i++)
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

int file_read(double *nx,double *ny, double *nz, double *px, double *py,double* pz)
{
    file_read_ASCII(nx,ny,nz,px,py,pz);
    return 0;
}
/**********************************************************************/
void writeVTKfiles(const int n, const double *nx,const double *ny,const double *nz,const double *px,const double *py,const double *pz,const double *bx,const double *by,const double *bz,const double *tx,const double *ty,const double *tz )
{
    int j;

    char vtk_data[200];
    ofstream output;
    sprintf(vtk_data,"%svtk_data_%d.vtk",prefix,n);
    output.open(vtk_data);
    output.precision(12);

    // header data for the VTK file
    output << "# vtk DataFile Version 3.0" << endl;
    output << "simulation data" << endl;
    output << "ASCII" << endl;
    output << "DATASET STRUCTURED_POINTS" << endl;
    output << "DIMENSIONS " << Lx << " " << Ly << " " << Lz << endl;
    output << "ASPECT_RATIO 1 1 1" << endl;
    output << "ORIGIN 0 0 0" << endl;
    output << "POINT_DATA " << LL << endl;
    output << "VECTORS n double" << endl; // output the director field
    for (j=0; j<LL; j++) {
        output << nx[j] << " " << ny[j] << " " << nz[j] << endl;
    }
    output << "VECTORS p double" << endl; // output the polarisation field
    for (j=0; j<LL; j++) {
        output << px[j] << " " << py[j] << " " << pz[j] << endl;
    }
    output << "VECTORS b double" << endl; // output the polarisation field
    for (j=0; j<LL; j++) {
        output << bx[j] << " " << by[j] << " " << bz[j] << endl;
    }
    output << "VECTORS curlcirc double" << endl; // output the polarisation field
    for (j=0; j<LL; j++) {
        output << tx[j] << " " << ty[j] << " " << tz[j] << endl;
    }
    output.close();
} // end writeVTKfiles
