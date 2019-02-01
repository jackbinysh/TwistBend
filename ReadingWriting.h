#include "twistbend.h"
using namespace std;

#ifndef READINGWRITING_H
#define READINGWRITING_H


int file_read( double *nx,double *ny, double *nz, double *px, double *py,double*pz);
int file_read_ASCII( double *nx,double *ny, double *nz, double *px, double *py,double*pz);

void writeVTKfiles(int n,double *nx,double *ny,double *nz,double *px,double *py,double *pz) ;
void writeBENDfiles(int n, double* nx,double* ny,double* nz) ;

#endif //READINGWRITING_H
