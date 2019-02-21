
	# All settings here for use on Benson Muite's laptop

	#define the complier
	COMPILER = gfortran
	# compilation settinmpngs, optimization, precision, parallelization
	FLAGS = -Ofast 

	# libraries
	LIBS =  -L. liblapack.a	libfftw3.a librefblas.a libfftw3_omp.a 
	

# source list for main program
SOURCES =  xyperiodic.f95 savedata.f90 inverse.f90 dispmodule.o omp_lib.f90

test: $(SOURCES)
		${COMPILER} -o navier.out $(FLAGS) $(SOURCES) $(LIBS)

clean:
	rm *.o
