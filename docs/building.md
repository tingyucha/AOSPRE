# Build the AOSPRE executable from source

## Prerequisites

> For important considerations before trying to build AOSPRE, please see
  the ["Preliminaries and Prerequisites"](preliminaries.md)
  documentation page.

## Download AOSPRE code

> ### Clone the AOSPRE github repository

> > The source code for the APAR Observing Simulator may be acquired by cloning the github reposoitory:
> >
> > <code>git clone https://github.com/NCAR/AOSPRE.git</code>
> >
> > this creates local directory
<nobr><i>AOSPRE</i></nobr>.  The source code is
found in directory
<nobr><i>AOSPRE/code/embed-crsim</i></nobr>.  The
code is written in Fortran.

> ### (alternate) Download and unpack an AOSPRE tarball.

> > [download tgz file and untar]

## Configure the Makefile

After downloading the code via <code>git clone</code> as above, and
insuring that required prerequisites are installed, change your
working directory to the source code directory:

<code>cd AOSPRE/code/embed-crsim</code>

This directory contains the Fortran source code for AOSPRE.  In addition,
the file _Makefile_ has instructions for the _make_ utility to compile
the source code.

Open the _Makefile_ with your preferred text editor.  Select options
most similar to your environment (by removing the leading <nobr>_#_
character</nobr> and spaces), or create a new set of options for your
environment following the examples already in the _Makefile_.  The
following options are used by the _Makefile_, and may need to be
customized for your particular computing environment:

### FC

  > <em>FC</em> specifies the command that invokes the Fortran
  compiler.  Common values are <code>gfortran</code> (for the GCC
  Fortran compiler), <code>ifort</code> (for the Intel Fortran
  compiler), and <code>nvfortran</code> (for the NVidia Fortran
  compiler).  For multi-processing builds, MPICH and OpenMPI provide
  specific wrappers to invoke the Fortran compiler (<i>e.g.</i>,
  <code>mpif90</code> or <code>mpifort</code>).  Depending on the
  computing environment, full path names to these compiler commands
  may be required.

### LD
  
  > <em>LD</em> specifies the comand that invokes the loader (linker)
  used in the final stage of building the executable, and is usually
  the same command used for the Fortran compiler (<em>FC</em> above).
  
### FORTRAN_FREEFORM

  > <em>FORTRAN_FREEFORM</em> specifies the particular compiler flags
    needed to recognize free-form Fortran source code style.  These
    flags are compiler dependent; examples for GCC, Intel, and Nvidia
    Fortran compilers are given in the default Makefile.

### NETCDF_INC

  > <em>NETCDF_INC</em> specifies the <i>include</i> directory where
    the NetCDF Fortran module file resizes.  This file is often called
    <i>netcdf.mod</i>.  The full path to the directory that contains
    the <i>netcdf.mod</i> file must be provided here.  If the
    <code><nobr>nf-config</nobr></code> command is available (available with a
    well-configured NetCDF installation), 
    the command <code><nobr>nf-config --fflags</nobr></code> may provide the appropriate directory.

  > The default Makefile assumes that the NetCDF <i>lib</i> and
    <i>include</i> directories are subdirectories of the same parent
    directory, specified by environment variable <i>NETCDF</i>.
    Particular installations of the NetCDF software may not follow
    that convention, so users may need to explicitly set the NetCDF
    <i>include</i> directory here.

### NETCDF_LIB

  > <em>NETCDF_LIB</em> specifies the <i>lib</i> directory where the
    NetCDF library files reside.  The full path to the directory that
    contains the NetCDF library files must be provided here.  If the
    <code><nobr>nf-config</nobr></code> command is available (available with a
    well-configured NetCDF installation), 
    the command <code><nobr>nf-config --flibs</nobr></code> may include the appropriate directory.

  > The default Makefile assumes that the NetCDF <i>lib</i> and
    <i>include</i> directories are subdirectories of the same parent
    directory, specified by environment variable <i>NETCDF</i>.
    Particular installations of the NetCDF software may not follow
    that convention, so users may need to explicitly set the NetCDF
    <i>lib</i> directory here.

### FFLAGS

  > <em>FFLAGS</em> specifies additional options that may be required
  for the Fortran compiler.  These options are often
  compiler-dependent.  The Makefile is set up so that <i>FFLAGS</i> will
  inherit values from <i>NETCDF_INC</i> and <i>FORTRAN_FREEFORM</i>.  Additional
  flags may include, for example, optimization flags.

### LDLIBS

  > <em>LDLIBS</em> specifies options instructing the loader to load
    specific libraries.  The Makefile is set up for <i>LDLIBS</i> to
    inherit the directory path specified in <i>NETCDF_LIB</i>, while
    explicitly referencing the <i>netcdff</i> library.  

### DEBUG

  > <em>DEBUG</em> specifies options to include additional debugging
  capabilities for the executable.  These options result in a
  substantially slower program execution, so they are usually turned
  off (commented out) unless needed for debugging purposes.

<!--

### LDFLAGS

  > <em>LDFLAGS</em> may be used to specify additional options to the
  loader.  If not needed it may be left blank.

-->

## Invoke the <i>make</i> command

> ### to compile

>> To compile, invoke the _make_ command:
>>
>> <code>make</code>
>>
>> The resulting build should create the executable file _a.out_.

> ### to clean the build directory

>> Invoking <code>make clean</code> removes all compiled code,
   including object files, Fortran module files, and the AOSPRE
   executable.  This may be useful if the compiler has changed, or if
   there were errors in a previous attempt to build the code.
