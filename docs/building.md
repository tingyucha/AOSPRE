# Build the AOSPRE executable from source

## Prerequisites

For important considerations before trying to build AOSPRE, please see the ["Preliminaries and Prerequisites"](preliminaries.md) documentation page.

## Obtain AOSPRE code

### Clone the AOSPRE GitHub repository

The source repository for AOSPRE is located on GitHub at [https://github.com/NCAR/AOSPRE/](https://github.com/NCAR/AOSPRE/).

It may be obtained by running the following command:

`git clone https://github.com/NCAR/AOSPRE.git` (HTTP)

`git clone git@github.com:NCAR/AOSPRE.git` (SSH)

The source code exists in the directory `AOSPRE/code/embed-crsim` and is written in Fortran.

## Configuring the Build
### Configuring with Makefile

After downloading the code via `git clone` as above, and insuring that required prerequisites are installed, change your working directory to the source code directory:

`cd AOSPRE/code/embed-crsim`

This directory contains the Fortran source code for AOSPRE.  In addition, the file *`Makefile`* has instructions for the `make` utility to compile the source code.

Open the *`Makefile`* with your preferred text editor.  Select options most similar to your environment (by removing the leading *`#`* character and spaces), or create a new set of options for your environment following the examples already in the *`Makefile`*.  The following options are used by the *`Makefile`*, and may need to be customized for your particular computing environment:

**`FC`** specifies the command that invokes the Fortran compiler.  Common values are `gfortran` (for the GCC Fortran compiler), `ifort` (for the Intel Fortran compiler), and `nvfortran` (for the NVidia Fortran compiler).  For multi-processing builds, MPICH and OpenMPI provide specific wrappers to invoke the Fortran compiler (*e.g.* `mpif90` or `mpifort`).  Depending on the computing environment, full path names to these compiler commands may be required.

**`LD`** specifies the comand that invokes the loader (linker) used in the final stage of building the executable, and is usually the same command used for the Fortran compiler (**`FC`** above).
  
**`FORTRAN_FREEFORM`** specifies the particular compiler flags needed to recognize free-form Fortran source code style.  These flags are compiler dependent; examples for GCC, Intel, and Nvidia Fortran compilers are given in the default Makefile.

**`NETCDF_INC`** specifies the *`include`* directory where the NetCDF Fortran module file resizes.  This file is often called *`netcdf.mod`*.  The full path to the directory that contains the *`netcdf.mod`* file must be provided here.  If the `nf-config` command is available (available with a well-configured NetCDF installation), the command `nf-config --fflags` may provide the appropriate directory.

The default Makefile assumes that the NetCDF *`lib`* and *`include`* directories are subdirectories of the same parent directory, specified by environment variable **`NETCDF`**. Particular installations of the NetCDF software may not follow that convention, so users may need to explicitly set the NetCDF *`include`* directory here.

**`NETCDF_LIB`** specifies the *`lib`* directory where the NetCDF library files reside. The full path to the directory that contains the NetCDF library files must be provided here. If the `nf-config` command is available (available with a well-configured NetCDF installation), the command `nf-config --flibs` may include the appropriate directory.

The default Makefile assumes that the NetCDF *`lib`* and *`include`* directories are subdirectories of the same parent directory, specified by environment variable **`NETCDF`**. Particular installations of the NetCDF software may not follow that convention, so users may need to explicitly set the NetCDF *`lib`* directory here.

**`FFLAGS`** specifies additional options that may be required for the Fortran compiler.  These options are often compiler-dependent.  The Makefile is set up so that **`FFLAGS`** will inherit values from **`NETCDF_INC`** and **`FORTRAN_FREEFORM`**.  Additional flags may include, for example, optimization flags.

**`LDLIBS`** specifies options instructing the loader to load specific libraries.  The Makefile is set up for **`LDLIBS`** to inherit the directory path specified in **`NETCDF_LIB`**, while explicitly referencing the *netcdff* library.  

**`DEBUG`** specifies options to include additional debugging capabilities for the executable.  These options result in a substantially slower program execution, so they are usually turned off (commented out) unless needed for debugging purposes.

<!--
**`LDFLAGS`** may be used to specify additional options to the loader.  If not needed it may be left blank.
-->

### Configuring with CMake

After cloning the repo, create a build directory and `cd` to that directory.

```bash
mkdir build
cd build
```

Call CMake with the desired options

```bash
cmake [-DPARALLEL=On] [-DCMAKE_BUILD_TYPE=Debug] <source_repo>
```

`<source_repo>` should point to the root of the AOSPRE source repo.

`-DPARALLEL=On` will attempt to find and link against an implementation of MPI.

`-DCMAKE_BUILD_TYPE=Debug` will build a debug version of the executable. This is significantly slower, but easier to debug. Other possible options besides `Debug` are `RelWithDebInfo` (optimizations on, but debug info included), `MinSizeRel` (smallest executable size), and `Release` (default, optimizations and no debug info)

Several environment variables can also be set to affect how CMake configures the build. Such environment variables can either be set using shell syntax (e.g. `FOO=bar cmake ..`) or with the `-D` option (e.g. `cmake -DFOO=bar ..`). Following are some environment variables of note:

**`FC`** sets the location of the Fortran compiler. Note: If this is set to `mpifort` or similar MPI wrapper, the `-DPARALLEL=On` option does still need to be set.

**`NETCDF_ROOT`** is the root of the NetCDF installation. This might need to be set if your NetCDF location is in a non-standard location.

**`MPI_HOME`** is the root of the MPI installation. If `PARALLEL` is turned on and CMake has trouble finding your installation of MPI, setting this variable should help.

## Building

Invoke build system command to compile.

```bash
make
```

The resulting executable should be called *`aospre`*

### Cleaning

```bash
make clean
```

This removes all artifacts creating during the build such as object files, module files, and the executable. This may be useful if the compiler has changed, or if there were errors in a previous attempt to build the code.
