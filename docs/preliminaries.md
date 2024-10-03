# Preliminaries and Prerequisites

Several topics must be considered before attempting to build and run
AOS.

## Linux computing environment recommended

> AOS has been developed and tested primarily in Linux environments.
  Other environments (e.g., MacOS, Docker/Singularity/Apptainer
  containers) may work for the knowledgeable and motivated user, but
  are not currently supported.

> Substantial desktop or laptop systems may have the resources for
  running relatively small, low-resolution experiments through AOS.
  More detailed, realistic AOS experiments will likely require
  powerful servers or supercomputer allocations.

## Fortran

> The core code for AOS is Fortran, and has been tested with Fortran
  compilers from GCC, Intel, and NVidia.
 
> The GCC gfortran compiler is commonly available through a Linux
  distribution's package manager (e.g., by installing package
  "gcc-gfortran"); Intel and NVidia both have Fortran compilers freely
  available for academic or personal use.  Other modern Fortran
  compilers may work as well, but have not been tested.

## NetCDF C and Fortran libraries:

> AOS relies on the NetCDF C and Fortran libraries, software developed
  by UCAR/Unidata <nobr>(<a href="http://doi.org/10.5065/D6H70CW6">http://doi.org/10.5065/D6H70CW6</a>).</nobr>

> Many Linux distributions include NetCDF libraries as an optional
  package (e.g., from EPEL libraries of Centos and Alma Linux
  distributions).  Users may also build NetCDF from source with code
  available from the above URL.  The NetCDF Fortran library should be
  compiled with the same compiler as will be used to build AOS.
  NetCDF libraries installed from a Linux software repository are most
  likely compiled with the GCC gfortran compiler.

## CR-SIM

> AOS uses lookup tables distributed with the Cloud Resolving Model
  Radar Simulator (CR-SIM).  See documentation page ["Acquire the
  CR-SIM tables"](acquiretables.md) for details on how to acquire
  these files.

## WRF or CM1 output files

> AOS requires model output files from WRF or CM1.  Best results come
  from high resolution (smaller than 1-km grid spacing) and high
  output frequency (smaller than 1 minute).  This implies a large data
  volume (on the order of gigabytes or even terabytes) which cannot be
  included in this AOS distribution.  Very small sample model output
  files, enough for a radar sweep or two, are included in the test
  directory.

  >> ###Morrison microphysics used in NWP model simulations

>> AOS currently works only with model output files from WRF or CM1
simulations that have used the Morrison microphysics scheme,
specifically mp_physics option 10 in the WRF namelist.  Other schemes
will likely be supported with later releases of AOS.  </li></ul>

## Optional software packages:

Additionally, users may find some optional packages useful:

>> ##MPICH or OpenMPI (_optional_)

   >> MPICH (<a href="https://www.mpich.org">www.mpich.org</a>) and
   OpenMPI (<a href="https://www.open-mpi.org/">ww.open-mpi.org</a>)
   are implementations of the Message Passing Interface (MPI)
   standard.  AOS uses MPI in its (optional) multi-processor
   configuration.  These software libraries are commonly available
   through your Linux distribution's package manager, or may be
   acquired through their respective web sites.

>> ## ncview (_optional_)

   >> ncview is commonly used as a quick-look tool for netCDF
   format files (which AOS produces). It is available at <a
   href="https://cirrus.ucsd.edu/ncview">https://cirrus.ucsd.edu/ncview</a>
   or may be available through your Linux distribution's package
   manager.

>> ## LROSE (_optional_)

   >> LROSE is a collection of software for working with Lidar,
   Radar, and Profiler datasets.  AOS creates files which are
   compatible with several utilities from the LROSE toolkit.  LROSE may
   be found at <a href="http://lrose.net">http://lrose.net</a>.
   
