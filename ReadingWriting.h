#include "twistbend.h"
using namespace std;

#ifndef READINGWRITING_H
#define READINGWRITING_H


int file_read( double *nx,double *ny, double *nz, double *px, double *py,double*pz);
int file_read_ASCII( double *nx,double *ny, double *nz, double *px, double *py,double*pz);

void writeVTKfiles(const int n, const double *nx,const double *ny,const double *nz,const double *px,const double *py,const double *pz,const double *bx,const double *by,const double *bz,const double *tx,const double *ty,const double *tz, const double* twist ,const double* FreeEnergy);
void print_Curve( double t, Link& Curve, string Name);
void writeCustomPoints(double t,double *nx,double *ny, double *nz,double *bx,double *by,double *bz);

#endif //READINGWRITING_H
