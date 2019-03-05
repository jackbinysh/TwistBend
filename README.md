# Bend Zeros in Twist Bend Nematics
------

Some general settings for the simulation are found in twistbend.h, but currently all internal settings for how the curve tracing works are inside the findbendzeros function in twistbend.cpp --- something to change in the future!

-----
#Settings in twistbend.h
- Nmax: number of timesteps to run
- dt: the timestep
- Lx, Ly,Lz: gridspacings. the code is nondimensionalised so dx = 1.
- K, thetah,qh,lambda, C, U: Twist bend parameters
- BC: use fixed or periodic BC's for the simulation. Currently the simulation works fine either way, the curve tracking wont work if it hits a boundary of any kind!
- vtkstepskip: How often to print the 3D grid data
-curvestepskip: How often to run curve extraction and print the curves (they take no memoery at all so this can be frequent).
-curvestarttime: When to start doing the tracing. Only start trying to trace when initial conditions have been smoothed away; im using around 3000 at the moment.

if you want to read things in:

- InitialisationMethod : read in from a vtk file or start fresh?
- director_filename etc: give the name from the file to read in from
- starttime: shift the 0 of time to agree with what timestamp the file you read in has.

----
# Initialisation Settings in twistbend.cpp

- You pick what texture you want to initialise in startconfig(): it contains two hopf textures, skyrmion lattice etc start configs. I have only really tested the HOPF initialisations.

----
# Whats going on inside FindBendZeros()

All the action happens inside FindBendZeros. The datatype "Link" which goes into this function is defined in twistbend.h; its just a vector of curves (one for each link component), themselves vectors of knotpoints. Knotpoints are another custom datatype defined in twistbend.h. They are just coordinates with additional data on them --- curve tangent vector, framing vector etc.

## list of settings, their location
All the settings are at the top of the function, each commented by their function. Its reproduced below:
-     threshold for detecting another bend zero (const double threshold) 
-     The size of the hemisphere in front of the current point which the code searches for the next point to go to. The data is on a grid of size 1, so it makes no sense for this to be way below 1. some O(1) number here is good. (const double sphereradius)
-     How many steps on the hemisphere the GSL minimizer should make. Ive found setting this too high makes the code more jagged. I don't reaaaaallly know why, perhaps being really picky about where the minimum is in "noisy" data causes it to run off. ( const int numgsliterations;)
-     the initial stepsize for the GSL minimzer. Something below 1, I use something well below it. ( const double initialgslstepsize;)
-     how close shoud our current point be to the initial point before we terminate and close the loop? Again, something O(1), I don't have an amazing feel. (const double terminationdistance;)
-     when we add the first point, we will be close to the start point and we will hit the termination condition above. To avoid this I just require us to be some number of points along the tracing before termination is allowed. This is the number of points (const int minnumpoints;)
-     when we do the two push offs, need to specify how far to push. roughly O(1) numbers again.  (double firstpushdist; double secondpushdist;)

## Outline of how the code works

To trace a curve, we find a point on the grid with bend magnitude below "threshold". thats the first point in a curve inside our Link object. To get the next point, we search a hemisphere in "in front of us" i.e a hemisphere centred on the orienting field vector at the current point. We do this search to stop us flowing off the bend zero, as we would do (numerically) just flowing along the integral curves of the orienting field. The size of this hemisphere is given by "spherradius". various minimiser settings are given by "numgsliterations, initialgslstepsize). Termination conditions on when this search stops are given by "terminationdistance, minnumpoints".

After a curve is traced out, the code marks every point connected to the curve by a path of values all below "threshold" --- it grabs the component of the volume defined by the condition (magnitudeb<"threshold") which is connected to the curve we just traced. The code then doesnt look here again for the start of a new curve. This is to get around the fact that magb ballooned outwards from the curves when we looked in paraview.

After this is done, the pushoffs happen, and theres some lowpass filtering which doesnt work at the mmoment.


