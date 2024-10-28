# APAR Observing Simulation, Processing, and Research Environment (AOSPRE)

The <b>APAR Observing Simulation, Processing, and Research Environment (AOSPRE)</b> 
is intended to create data sets that simulate results from flights of the 
<b>Airborne Phased Array Radar (APAR)</b> observing platform.  AOSPRE relies 
on high-resolution simulations of real or idealized weather events as provided 
by, e.g., the <b>Weather Research and Forecasting Model</b> (WRF; Skamarock
et. al., 2019) or the <b>Cloud Model version 1</b> (CM1; Bryan and
Morrison, 2012).  AOSPRE takes a user-specified flight plan through such
atmospheric simulations, together with radar scanning instructions
(e.g., series of rotation and tilt angles), to determine which
portions of the atmosphere APAR pulses might "see" during such a
flight.  Details of the thermodynamical and microphysical fields at
these locations are extracted from the model output, and passed to
routines from the <b>Cloud-resolving Radar Simulator</b> (CR-SIM; Oue
et al., 2020).  The CR-SIM routines compute various radar moments
(reflectivity, radial velocity, etc.) from this model data, and AOSPRE
writes these results in the <b>CfRadial File Format</b> (Dixon and
Lee, 2016).

Some main goals of AOSPRE include:
<ul>
<li>providing researchers with realistic
(though simulated) data sets as an aid to planning APAR flight
strategies, scanning strategies, and analysis techniques.</li>
<li>creating datasets
useful for testing and developing APAR software processing tools.</li>
<li>creating or updating various radar algorithms</i>
</ul>

The atmospheric models used as input to AOSPRE should ideally be run at
resolutions and output time frequencies that are comparable to the
spatial and temporal scales of radar observations.  In practice,
simulations with grid spacing of 100 to 1000 m, and time output of a
few seconds to a minute, are probably the resolutions at which AOSPRE is
expected to be applied.

## <h3>APAR</h3>

APAR is the Airborne Phased Array Radar, an observing platorm currently under development at NSF NCAR. Learn more about APAR at [www.eol.ucar.edu/airborne-phased-array-radar-apar-0](https://www.eol.ucar.edu/airborne-phased-array-radar-apar-0).

## <h3>CR-SIM</h3>

The Cloud-resolving Radar Simulator (CR-SIM) derives radar moments
from high resolution atmospheric model output.  Learn more about
CR-SIM at <a
href=https://you.stonybrook.edu/radar/research/radar-simulators>you.stonybrook.edu/radar/research/radar-simulators</a>.

<!--
The AOSPRE uses a radar simulator known as the Cloud-resolving Radar
Simulator (CR-SIM, Oue et al. 2020), which was developed by
researchers at one of NSF NCAR’s partner institutes, SUNY Stony Brook.
-->

## <h3>WRF</h3>

The Weather Research and Forecasting model (WRF) is a numerical
weather prediction model used extensively in weather research and
real-time forecasting operations.  Learn more about WRF at <a
href=https://www.mmm.ucar.edu/models/wrf>www.mmm.ucar.edu/models/wrf</a>.

## <h3>CM1</h3>

The Cloud Model version 1 (CM1) is designed for detailed idealized
simulations of weather systems at resolutions of hundreds of meters to
a few kilometers.  Learn more about CM1 at <a
href=https://www2.mmm.ucar.edu/people/bryan/cm1/>www2.mmm.ucar.edu/people/bryan/cm1</a>.

## <h3>CfRadial File Format</h3>

The CfRadial File format is a CF-Compliant NetCDF format for radar and
lidar moments in radial coordinates.  Learn more about CfRadial at <a
href=https://www.eol.ucar.edu/sites/default/files/files_live/private/CfRadialDoc.v1.4.20160801.pdf>www.eol.ucar.edu/sites/default/files/files_live/private/CfRadialDoc.v1.4.20160801.pdf</a>.



<!--
# Welcome to MkDocs

For full documentation visit [mkdocs.org](https://www.mkdocs.org).

## Commands

* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs build` - Build the documentation site.
* `mkdocs -h` - Print help message and exit.

## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages, images and other files.
-->
