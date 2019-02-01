#include "twistbend.h"
using namespace std;

#ifndef READINGWRITING_H
#define READINGWRITING_H


int file_read(double n, double *nx,double *ny, double *nz, double *px, double *py,double*pz);
int file_read_ASCII(double n, double *nx,double *ny, double *nz, double *px, double *py,double*pz);

void writeVTKfiles(void);
void writeBENDfiles(void);

#endif //READINGWRITING_H
