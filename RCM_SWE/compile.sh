# delete first old compilations of modules
if [ -f *.o ]; then
  rm *.o *.mod
fi

# compile modules

gfortran -c -pg arrays.f90
gfortran -c -pg workspace.f90
gfortran -c -pg constants.f90
gfortran -c -pg additional.f90
gfortran -c -pg input.f90
gfortran -c -pg fluid_solver.f90 
gfortran -c -pg output.f90 


# compile main program 
gfortran -g -fbounds-check -pg -o main main.f90 additional.o arrays.o workspace.o constants.o fluid_solver.o output.o input.o 
# remove module compilations
rm *.o *.mod
