# Acquire the CR-SIM lookup tables

AOSPRE depends on a large set of lookup tables provided with the CR-SIM
software download.  The tables require more than 20GB of disk space,
and are not included in the AOSPRE download.  The CR-SIM package can be
downloaded from <a href="https://www.bnl.gov/cmas/cr-sim.php">https://www.bnl.gov/cmas/cr-sim.php</a>.
Download the full software package.

## Download the CR-SIM software package

The recommendation is to move to the top level
<code></nobr>APAR-Observing-Simulator/</nobr></code> directory, and download the CR-SIM package
into this directory.  The download file is a compressed tar file more
than 12GB in size.  The file will be named according to the version of CR-SIM downloaded:  <code><nobr>crsim-<<i>version</i>>.tar.gz</nobr></code>.

## Unpack CR-SIM tar file

Once downloaded untar the crsim tar file:

> <code>tar -xvf crsim-<<i>version</i>>.tar.gz</code>

or you may choose to unpack only the tables needed for AOSPRE:

> <code>tar -xvf crsim-<<i>version</i>>.tar.gz crsim-<<i>version</i>>/share/crsim/aux</code>

Either way, contents of the tar file will be unpacked into a directory named <code><nobr>crsim-<<i>version</i>>/</nobr></code>.

Once unpacked, the downloaded crsim tar file may be removed.

## Create a link to the CR-SIM tables

The CR-SIM tables are located in directory
<code><nobr>crsim-<<i>version</i>>/share/crsim/aux/</nobr></code>.
When AOSPRE runs, it assumes it will find tables under a directory called <code>aux/</code> immediately
above the working directory.  Further examples in this documentation
assume that execution of AOSPRE will occur in a working directory under
the top-level <code><nobr>APAR-Observing-Simulator/</nobr></code> directory.  For this reason,
it is recommended to create a soft link in this top-level directory to
the <code>aux/</code> tables directory:

> <code>ln -s crsim<<i>version</i>>/share/crsim/aux .</code>

If you set up to run AOSPRE in other locations, remember to link the
CR-SIM <code>aux/</code> directory with a soft link <code>aux</code>
just above your AOSPRE run directory.
