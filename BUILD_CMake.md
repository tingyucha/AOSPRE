# Build APAR Observing Simulator with CMake

1. Install `cmake`. On CentOS 8, install `cmake3`. 
2. Create a directory for your build. For the sake of example, we'll assume the build directory is in the root of the code checkout. Change any instance of `..` with the path to where you checked out the code if your situation is different.
3. `cd build`
4. `cmake [-DPARALLEL_AOS=On] ..`
4a. Add `-DPARALLEL_AOS=On` if you wish to build it with MPI. An installed instance of MPI is required.
5. `make`