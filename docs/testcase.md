# Test case mini-tutorial

A small case has been set up for a quick and simple test run. Input files,
configuration files, and namelists have all been set up for this test
case.

From the top-level <code><nobr>AOSPRE</nobr></code>
directory, move to the working directory for the test case:

> <code>cd test</code>

This directory will be your working directory for running the test case.

Several components are required for a successful run of AOSPRE.  All of
these components have been prepared in a simple configuration for this
test case.  Do a listing of the current working directory:

> <code>ls -l</code>

The various files required for a run of AOSPRE are described below.

## WRF-formatted model output files

The principal data from which simulated radar sweep will be
constructed are contained in model output files from a high-resolution
weather model.  In this example, the CM1 model running with 500-m grid
spacing simulated an idealized supercell.  For this example, three
vertical slices from the full model run have been extracted at two
output times.  This is enough for a couple of RHI sweeps as a
demonstration.  These data are found in the "wrfout" files in this directory:

> <code>wrfout_test_4650.nc<br></code>
> <code>wrfout_test_4660.nc</code>

The <code>4650</code> and <code>4660</code> numerals in these sample
filenames indicate the number of seconds since model initialization.
You can use the NetCDF utility program <code>ncdump</code> to examine
the contents of these model output files, and the <code>ncview</code>
graphics utility program (if available) to get a quick look at the
meteorological fields in these files.

## Scanning tables

AOSPRE relies on user-provided scanning tables to define the beams and
sweeps that will be used to sample the atmospheric state as
represented in the model output files.  For this example, these files are:

> <code>scanning_bot.txt</code><br>
> <code>scanning_lhs.txt</code><br>
> <code>scanning_rhs.txt</code><br>
> <code>scanning_top.txt</code>

where "bot", "lhs", "rhs", and "top" indicate the bottom panel of the
APAR platform, the left-hand (port) side panel, the right-hand
(starboard) side panel, and the top panel.

The bulk of these scanning table files is a series of beam directions,
expressed in terms of rotation and tilt angles in the "type-Z"
coordinate system (<i>cf.</i> CfRadial documentation), These particular files
have been set up for this test case to produce sweeps in the plane
perpendicular to the direction of flight.

## CR-SIM configuration table

The components of CR-SIM used by AOSPRE use a configuration file to
specify characteristics of the simulated radar.  This configuration file is:

> <code>CONFIG_crsim</code>

Details of the CR-SIM configuration may be found in the CR-SIM
documentation.  Note that only certain entries in this configuration
setup are used by AOSPRE; the entries related to radar position, lidar,
Micropulse lidar (MPL) and postprocessing are not used by AOSPRE.

## AOSPRE namelist files

AOSPRE processing is controlled by one or more Fortran namelist files:

> <code>namelist.BOT</code><br>
> <code>namelist.LHS</code><br>
> <code>namelist.RHS</code><br>
> <code>namelist.TOP</code>

Where the bottom, left-hand side, right-hand side, and top panels each have their own namelist files.

The namelist files specify the flight path and radar scanning strategy
for a flight, and the input and output files for the AOSPRE processing..
See the [namelist documentation](namelist.md) for details of the
namelist entries.  Several of the critical namelist entries are described below:

### wrf_glob_pattern

> The <code>wrf_glob_pattern</code> entry defines how AOSPRE finds the
  atmospheric model output files.  Any file matching the glob pattern
  (using standard wildcard matching where '?' matches any single
  character and '*' matches any number of characters) will be
  interpreted as a model output file.  In this example, the
  <code>wrf_glob_pattern</code> has been set up to match the names of
  both WRF output files provided for the test.

### output_filename_format_string

> The <code>output_filename_format_string</code> determines how the
  CfRadial-formatted files produced by AOSPRE will be named.  This entry
  must specify a Fortran format string with exactly two 'A' format
  specifiers, used to build output file names.  The 'A' format
  specifiers will resolve to the beginning and ending time of the
  sweeps in any particular output file.  With each panel of the APAR
  platform producing its own set of output volumes, the file name
  construction should be unique for each panel.  Compare the
  <code><nobr>output_filename_format_string</nobr></code> for each
  namelist in the test directory.

### leg_initial_time

> The <code>leg_initial_time</code> entry specifies the time (in
  seconds from model initialization) to start a data-collection flight
  leg.  Since the example WRF output files contain only 10 seconds of
  data (model output times at 4650 to 4660 seconds past
  initialization), we've set up the namelists for an 8-second flight
  leg starting at time 4651 seconds past initialization.

### leg_time_seconds

> The <code>leg_time_seconds</code> entry specifies the duration in
  seconds of the data-collection fligh leg.  As mentioned above, only
  a minimal set of data is provided for this small test.  An 8-second
  flight will fit within the data provided.

### flight_waypoints_x, y, and z

> The entries <code><nobr>flight_waypoints_x</nobr></code>,
  <code><nobr>flight_waypoints_y</nobr></code>, and
  <code><nobr>flight_waypoints_z</nobr></code> specify the flight path
  taken through the model output files.  The X and Y coordinates are
  specified in terms of grid-point indices in the model grid (x and y
  both equal 1 at the lower left grid point).  The Z coordinate may be
  specified in meters above MSL, or as pressure in hPa, depending on
  the setting of <code><nobr>flight_level_coordinate</nobr></code>.

> The test dataset (in the wrfout files) consists of three vertical
  slices, in the south-north direction, from a larger model grid.  The
  dimensions of this subset grid (X x Y) are 3 x 400.  Our initial
  flight is set up to move from west to east at y=105, south of the
  main body of the storm, at 2000 meters above sea level.

> When running AOSPRE with for multiple APAR panels, the namelists for
  each panel must be set up with the same flight path and flight time
  information.

### CRSIM_Config

> The <code>CRSIM_Config</code> entry must specify the file name of
  the CR-SIM configuration file, in this case, "CONFIG_crsim".

### scanning_table

> The <scanning_table> entry specifies the name of the file that
  contains the details of the radar beam directions for the panel
  represented by the namelist.  Compare the
  <code><nobr>scanning_table</nobr></code> entries for the four
  namelists provided.

## Run the AOSPRE program

The AOSPRE executable takes as command-line arguments the names of one or
more namelist files.  As a first run of AOSPRE, invoke a flight using the
original namelist settings, with only the left-hand side panel
active:

> <code>../code/embed_crsim/aospre namelist.LHS</code>

If you get an immediate segmentation fault, you may need to unlimit
the stacksize in your terminal session:

> <code>ulimit -s unlimited</code> for bash shell,

> <code>unlimit stacksize</code> for csh shell,

and try again.

A successful run will result in several CfRadial-formatted files from
the left-hand side panel of the APAR platform.  E.g.:

<code>LHS_rhi_20010101_011731.092_to_20010101_011732.149.nc</code><br>
<code>LHS_rhi_20010101_011732.244_to_20010101_011733.301.nc</code><br>
<code>LHS_rhi_20010101_011733.396_to_20010101_011734.453.nc</code><br>
<code>LHS_rhi_20010101_011734.548_to_20010101_011735.605.nc</code><br>
<code>LHS_rhi_20010101_011735.700_to_20010101_011736.757.nc</code><br>
<code>LHS_rhi_20010101_011736.852_to_20010101_011737.909.nc</code><br>

These contents of these files may be examined with the NetCDF utility
program <code>ncdump</code>, and may be explored with the
<code>ncview</code> utility program.  Since these files are in radar
coordinates, it is often difficult to make sense of the quick-look
views provided by the <code>ncview</code> program.  Other graphics
capabilities that work with the radar coordinates, such as those in
the LROSE package, are beyond the scope of this mini-tutorial.

Note also the file <code><nobr>flightpath.ascii</nobr></code>.  This
file preserves details of the aircraft's flight through the
simulation, with locations, directions, and wind values at the
aircraft location, every tenth of a second.

## Modify the flightpath to run again

To explore the multi-panel capability, edit all the namelists to move
the <code><nobr>flight_waypoints_y</nobr></code> entries to somewhere
near the middle of the supercell, at Y index about 150, where the
various panels may all see interesting weather.

Rerun the flight with the left panel:

> <code>../code/embed_crsim/aospre namelist.LHS</code>

And do a repeat flight with the right panel:

> <code>../code/embed_crsim/aospre namelist.RHS</code>

Note the additional RHS CfRadial files produced by this second run.

You can also run multiple APAR panels in the same run, by specifying
multiple namelists on as command-line arguments to the AOSPRE executable.

E.g., try:

> <code>../code/embed_crsim/aospre namelist.RHS namelist.LHS</code>

or

> <code>../code/embed_crsim/aospre namelist.BOT namelist.TOP</code>

or

> <code>../code/embed_crsim/aospre namelist.BOT namelist.TOP namelist.LHS namelist.RHS</code>

and note the output generated.

## Further exploration

As noted before, the minimal WRF output dataset allows only limited
AOSPRE runs.  But there are a variety of runs you can do to further
explore the data or familiarize yourself with AOSPRE configuration.  E.g.:

<ul>
<li>Try setting up flights north of the storm.
<li>Try setting up flights from west to east.
<li> See how results change with flights at different levels.
<li> Change the direction of the flight and the corresponding orientation of the beams in the scanning tables and still get sweeps within the slices of the WRF output data.
</ul>
