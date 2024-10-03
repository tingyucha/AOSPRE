#!/bin/bash

# RUN_AOS.sh
#
# Brad Klotz, NCAR Project Scientist
# Created: 06/16/2021

# This program  will initiate the APAR Observing Simulation (AOS)
# system code. You will need to specify a few items in order to complete the
# code requirements. You will need to configure your namelist
# file prior to running this script for now. That functionality
# may change in the future.

# The script will read information that you sepcify at the top of
# this script, which include: a base directory, operating
# directory, scan type, namelist file, experiment_directory,
# test directory, and output directory.

# Recommended script execution - creates a log file and runs in the background
# ./RUN_AOS.sh > out_aos.log &

# Base Directory: This should be your main AOS directory you've created
#                 example: /export/wind1/bradklotz/git/APAR-Observing-Simulator
#
# Operating Directory: This is the directory containing your namelist and config files
#                 example: /export/wind1/bradklotz/AOS/tests/WRF_supercell_100m
#
# Scan Mode: This is simply the type of scan - surveillance or rhi
#
# Namelist file: This is the file that controls the functionality of the AOS code
#                 example: namelist.surveillance
#
# Experiment Directory: If you are performing several runs as part of a larger
#                       experiment, enter the experiment directory here. Otherwise leave
#                       as an empty string.
#                       example: experiments/test_01
#
# Test Directory: This is the name of the test-specific directory you
#                 want to create. It will typically have information about the
#		  particular AOS run you are executing.
#		  example: supercell_100m_10s_Cband_1kmAGL
#
# Output Directory: This is the combined strings of the Base Directory, the
#                   Experiment Directory, and the Test Directory.
#                   No user input is needed to create this directory.
#

##### USER-DEFINED VARIABLES #####

# Base directory
base_dir=/export/wind1/bradklotz/APAR-Observing-Simulator

# Operating directory
operating_dir=/export/wind1/bradklotz/AOS/run_tests/WRF_supercell_100m

# Output base directory
output_base_dir=/export/wind1/bradklotz/AOS

# Scan mode: rhi or surveillance
scan_mode="surveillance"

# Namelist file
namelist_file=namelist.surveillance

# Experiment directory - if none, leave as empty string
experiment_dir="output_tests"
#experiment_dir="experiments/distance_tests"

# Test directory
# This will typically provide some specific information
# about your run of the AOS - such as the storm type, resolution
# aircraft initial height, microwave band, etc.
test_dir="supercell_100m_10s_1kmAGL_Cband_2.0BW_x215_y315_4490-5090s_FVEL"

#### No user input needed beyond this point ####

# Output directory
# First check that the experiment directory exists (if it is not an
# empty string)

if [ ! -z "$experiment_dir" ]; then
   
   if [ ! -d ${output_base_dir}/${experiment_dir} ]; then
      mkdir -p ${base_dir}/${experiment_dir}
   fi

   output_dir=${output_base_dir}/${experiment_dir}/${test_dir}
else
   output_dir=${output_base_dir}/${test_dir}
fi

# Get your current working directory
current_dir=$(pwd)

echo "Running the AOS code..."
echo "Your current working directory is: " ${current_dir}
echo "Your operating directory is: " ${operating_dir}
echo "Your scan mode is: " ${scan_mode}
echo "Your namelist file is: " ${namelist_file}
echo "Your output directory is: " ${output_dir}

##### END USER VARIABLES #####

##### EXECUTION OF AOS CODE #####

# Verifying the base and operating directories  exist
if [[ ! -d ${base_dir} || -z "$base_dir" ]]; then
	echo "Your base directory is not correct...please fix and re-run."
	echo "Exiting..."
	exit
fi

if [[ ! -d ${operating_dir} || -z "$operating_dir" ]]; then
	echo "Your operating directory is not correct...please fix and re-run."
	echo "Exiting..."
        exit
fi

# Create the output directory if it doesn't exist
echo ""
echo "All specified directories are verified...Proceeding with AOS execution..."

if [ ! -d ${output_dir} ]; then
	   mkdir ${output_dir}
fi

# Go to the operating directory
echo ""
echo "The current working directory is now the operating directory: " ${operating_dir}
cd ${operating_dir}

# Execute the code
echo ""
echo "Initiating AOS code..."
echo ""
echo "And off we go..."
echo ""
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "%%        AOS SPECIFIC OUTPUT BEGINS        %%"
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

./${base_dir}/code/embed-crsim/a.out ${namelist_file}

echo ""
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "%%         AOS SPECIFIC OUTPUT ENDS         %%"
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

# Once complete, move created files to output directory
echo ""
echo "Processing complete...moving files to output directory..."
if [ ${scan_mode} == 'rhi' ]; then

   # Creating a directory within your output for the particular scan mode
   if [ ! -d ${output_dir}/${scan_mode} ]; then
       mkdir ${output_dir}/${scan_mode}
   fi

   mv *RHI.*.nc ${output_dir}/${scan_mode}
   cp 'flightpath.ascii' ${output_dir}/${scan_mode}
   cp ${namelist_file} ${output_dir}/${scan_mode}
else

   # Creating a directory within your output for the particular scan mode
   #newstr1="_hifreq"
   newstr1=""
   if [ ! -d ${output_dir}/${scan_mode}${newstr1} ]; then
       mkdir ${output_dir}/${scan_mode}${newstr1}
   fi

   mv surveillance_*.nc ${output_dir}/${scan_mode}${newstr1}
   cp 'flightpath.ascii' ${output_dir}/${scan_mode}${newstr1}
   cp ${namelist_file} ${output_dir}/${scan_mode}${newstr1}
fi

# Return to the original directory when you started
echo ""
echo "Returning to your original working directory: " ${current_dir}
cd ${current_dir}

echo ""
echo "AOS code execution is complete...Done!"


