# APAR Observing Simulation, Processing, and Research Environment (AOSPRE)

The APAR Observing Simulation, Processing, and Research Environment (AOSPRE) is intended create data sets that simulate results from flights of the Airborne Phased Array Radar (APAR) observing platform.  It simulates flights and radar returns from weather events as simulated by WRF and CM1 numerical weather prediction models.

## Preliminaries

Before attempting to build or run AOSPRE, please see the ["Preliminaries and Prerequisites"](docs/preliminaries.md) documentation page.

## Building

For instructions on acquiring and building the AOSPRE code, see the ["Build the AOSPRE executable"](docs/building.md) documentation page.

## Testing

For instructions on testing the AOSPRE executable with some sample data, see the ["Run the test case"](docs/testcase.md) documentation page.

## Run your own experiments with your own model output

- Acquire NWP model output files
- Design flight plan
- Design scanning strategy (scanning tables)
- Configure AOSPRE namelist
- Configure CR-SIM configuration file
- Link or copy CR-SIM lookup tables
- Unlimit stacksize
- run AOSPRE
