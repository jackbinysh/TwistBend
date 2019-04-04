/*
    Copyright 2018
    Jack Binysh <j.binysh@warwick.ac.uk>,
    Gareth Alexander <g.p.alexadner@warwick.ack.uk>
    This software is provided under a BSD license. See LICENSE.txt.

    This file contains:
    -the datastuctures used to encode the link, and the viewpoint to look at it from.
    - a list of hardcoded vectors we succesively pick n_infty from if our initial choice was no good
    - the main functions to actually compute the solid angle function.
    - a few little helper math/logic functions

*/

#ifndef SOLIDANGLE_H
#define SOLIDANGLE_H

#include "../twistbend.h"
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


// A list of hardcoded directions to succesively try of the user inputted one gives n.n_infty exceeding the threshold.
// This is faster than actually selecting at random and sufficient for our purposes. And we don't have to worry about random numbers + threads = state that differs run to run
// read in threes:
const double hardcodedvectors [] = {1,0,0  ,  0,1,0  ,  0,0,1   ,  0,0.70710678118,0.70710678118  ,  0.70710678118,0,0.70710678118  ,  0, 0.70710678118,0.70710678118};
const int numhardcodedvectors = 6;


// Runs through our 3D grid omega, and computes the solid angle at each point.
void ComputeSolidAngleAllPoints(const Link& Curve, double *omega);

// All the real work of the code happens in this function - given a viewpoint, give me the solid angle subtended by the knot at that viewpoint.
// this function implements the formula for this solid angle, \omega({\bf x}) = \int \frac{{\bf n}_{\infty} \times {\bf n} }{1+{\bf n}_\infty \cdot {\bf n}} \cdot \mathrm{d}{\bf n} ,
// with a user defined threshold checking whether we are too close to a surface of discontinuity. If we are, we first swap from n_infty -> -n_infty and try again. If that doesn't work
// we pick a new one at random, and keep trying. Its output is standardised to the range 0 to 4 pi
double ComputeSolidAngleOnePoint(const Link& Curve, const viewpoint& View);

// new function for computing minimum distance to the curve
double ComputeDistanceOnePoint(const Link& Curve, const viewpoint& View);
// new function for computing longitudinal phase along the curve -- not adapted to components yet!!
double ComputeLongitudinalPhase(const Link& Curve, const viewpoint& View);
// function to get the solid angle framing
void ComputeSolidAngleFraming(Link&  Curve, double *phi);


#endif
