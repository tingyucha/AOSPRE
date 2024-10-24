# Configuration

A number of configuration steps are required before running AOSPRE.
These steps are outlined below.  It may be helpful also to review the
the ["Preliminaries and Prerequisites"](preliminaries.md)
documentation page.

#########################################################################

## Collect weather model output files

It is up to the user to either produce their own model simulation
output or obtain output from another source. Without it, the AOSPRE
cannot run. It is recommended to use a model output grid spacing no
larger than 1 km, as higher resolutions are better able to represent
atmospheric features and variability at the scales needed for AOSPRE to
produce realistic simulated radar scans.

Make sure the model output files are in a location accessible from the
AOSPRE run directory (i.e., not on a remote server).  The files will be
found at runtime by AOSPRE using the "wrf_glob_pattern"
[namelist](namelist.md) option.  Check to to see that the
"wrf_glob_pattern" option will find the desired model output files.

#########################################################################

## Design a flight plan

The flight path through the simulation is set by
[namelist](namelist.md) options "flight_waypoints_x",
"flight_waypoints_y", and "flight_waypoints_z".  The x and y
coordinates are specified in by model grid point index, and the z
coordinate may be specified in pressure (hPa) or height MSL (m)
depending on the "flight_level_coordinate" [namelist](namelist.md)
option.  The flight path is traversed by straight line paths between
successive (x,y,z) coordinates.

Designing a flight plan can be one of the more difficult parts of
using AOSPRE.  Users may find the ncview utility useful for examining the
NWP simulation files and finding gridpoint coordinates at appropriate
times and locations.  The "flight_planner" shell script (in the
"scripts" directory) might also be useful.

Once a flight plan is designed, users may want to run AOSPRE in a mode to
create a few surveillance PPI scans along the flight path, and use
these scans to verify that the flight plan is suitable for sampling
the regions of interest in the NWP simulation, before running more
extensive and compute-intensive flights and scanning strategies.

#########################################################################

## Configure CR-SIM options

CR-SIM options are configured in the CR-SIM configuration file.  The
file name for the CR-SIM configuration file is referenced in the AOSPRE
"CRSIM_Config" [namelist](namelist.md) option.

Only specific settings in the CR-SIM configuration file are used by
AOSPRE.

#########################################################################


## Configure AOSPRE namelist

* See ["Namelist"](namelist.md) documentation page.







#########################################################################

## Configure scanning tables

AOSPRE scanning tables are organized in a two-column format text file and
contain information about the primary axis of rotation, sweep mode,
and the information contained in each column (i.e., rotation and tilt
angles). The basic scan files use Z as the primary axis of rotation
with beams pointing at 240&deg; and 275&deg; relative to the aircraft
fuselage. Tilt angles for the RHI scan are between -53&deg; and
53.5&deg; at intervals of 1.5&deg; for a total of 72 beams per
rotation angle. The surveillance scan mode uses a tilt of 0&deg; with
a full azimuthal rotation at intervals of 1.5&deg; for a total of 241
beams.

The user can create their own scan tables as desired.

Sample scanning tables include:

* Scanning_Table_LHS_PPI_sector_Z (for surveillance mode)
* Scanning_Table_LHS_Multiplexing_RHI_Zprime_ordered (for left panel RHI)
* Scanning_Table_RHS_Multiplexing_RHI_Zprime_ordered (for right panel RHI)

The scanning table file name is referenced in the AOSPRE namelist.

#########################################################################

## Program execution

Usage

*executable* *namelist* [*namelist* [*namelist*]]

MPICH usage

mpirun -n <np> *executable* *namelist* [*namelist* [*namelist*]]

#########################################################################

## Output

# CfRadial-format NetCDF files

CfRadial format is described in CfRadial document.


