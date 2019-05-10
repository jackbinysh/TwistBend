CXX=g++
CXXFLAGS= -O3 -fopenmp  
LDLIBS=  -lgsl -lgslcblas -lm  
LDFLAGS = -O3 -fopenmp 
OBJS= TriCubicInterpolator.o ReadingWriting.o twistbend.o  SolidAngle/SolidAngle.o SolidAngle/InputOutput.o SolidAngle/Geometry.o SolidAngle/Constants.o 
DEPS= TriCubicInterpolator.h ReadingWriting.h twistbend.h  SolidAngle/SolidAngle.h SolidAngle/InputOutput.h SolidAngle/Geometry.h SolidAngle/Constants.h 

%.o: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

all:Code clean

Code:$(OBJS)
	$(CXX) -o twistbend $(OBJS) $(LDLIBS) $(LDFLAGS)

.PHONY: clean

clean:
	rm -f SolidAngle/*.o
	rm -f *.o


