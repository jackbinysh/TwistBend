#include "twistbend.h"
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>
#include <time.h>
using namespace std;

#ifndef Solid_Angle_Link_H
#define Solid_Angle_Link_H

struct viewpoint
{
    double xcoord;
    double ycoord;
    double zcoord;
};

/***********************Functions for outputting the solid angle*************************/

double SolidAngleCalc(const Link& Curve, const viewpoint& View);
void ComputeSolidAngle(double *phi, const struct Link& Curve);
void ComputeSolidAngleFraming(double* phi, Link &Curve);

#endif
