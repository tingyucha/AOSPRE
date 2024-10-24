#!/bin/bash

# RUN_AOSPRE_test.sh
#
# Brad Klotz, NCAR Project Scientist
# Created: 06/16/2021
# Modified: 11/15/2022 - for testing with the parallel computing mode

# This program  will initiate the APAR Observing Simulation,
# Processing, and Research Environment (AOSPRE)
# system code. You will need to specify a few items in order to complete the
# code requirements. You will need to configure your namelist
# file prior to running this script for now. That functionality
# may change in the future.

# The script will read information that you sepcify at the top of
# this script, which include: a base directory, operating
# directory, scan type, namelist file, experiment_directory,
# test directory, and output directory.

# Recommended script execution - creates a log file and runs in the background
# ./RUN_AOSPRE_test.sh fpar ftime > out_aospre.log &
#
# Even though the <fpar> and <ftime> variables are listed elsewhere, they need to be included
# in the command as well. So <fpar> indicates whether you want to run parallel or not and <ftime>
# is the flight time in seconds.
#
# fpar: Enter 1 for parallel or 0 for serial
# ftime: Enter the length of the flight in seconds
#
# Other entries of interest
#
# Base Directory: This should be your main AOSPRE directory you've created
#                 example: /export/wind1/bradklotz/AOSPRE
#
# Operating Directory: This is the directory containing your namelist and config files
#                 example: ${basedir}/code/WRF_supercell_100m
#
# Scan Mode: This is simply the type of scan - surveillance or rhi
#
# Namelist file: This is the file that controls the functionality of the AOSPRE code
#                 example: namelist.surveillance
#
# Experiment Directory: If you are performing several runs as part of a larger
#                       experiment, enter the experiment directory here. Otherwise leave
#                       as an empty string.
#                       example: experiments/test_01
#
# Test Directory: This is the name of the test-specific directory you
#                 want to create. It will typically have information about the
#		  particular AOSPRE run you are executing.
#		  example: supercell_100m_10s_Cband_1kmAGL
#
# Output Directory: This is the combined strings of the Base Directory, the
#                   Experiment Directory, and the Test Directory.
#                   No user input is needed to create this directory.
#

##### USER-DEFINED VARIABLES #####

# Base directory
base_dir="/export/wind1/bradklotz/AOSPRE"

# Operating directory
operating_dir=${base_dir}/code/WRF_supercell_100m
#operating_dir=${base_dir}/code/WRF_supercell_500m
#operating_dir=${base_dir}/code/WRF_squall_line_500m
#operating_dir=${base_dir}/code/WRF_squall_line_100m
#operating_dir=${base_dir}/code/WRF_TC_500m
#operating_dir=${base_dir}/code/WRF_HNR1_1km
#operating_dir=${base_dir}/code/WRF_SL_PECAN_15JUL

# Scan mode: rhi or surveillance
scan_mode="rhi"

# Namelist file
#namelist_file=namelist_0.5_AK
namelist_file=namelist.rhi
#namelist_file=namelist_rhi_supercell_100m_0.1kmAGL_WestEast_x215_y655_singleleg_10minute

# Experiment directory - if none, leave as empty string
experiment_dir="individual_tests"
#experiment_dir="experiments/distance_tests"

# Experiment resolution (100m or 500m)
exp_res=500

# Test directory
# This will typically provide some specific information
# about your run of the AOSPRE - such as the storm type, resolution
# aircraft initial height, microwave band, etc.
#test_dir="supercell_100m_10s_1kmAGL_Cband_2.0BW_x215_y605_4490-5090s_FVEL_Aft360_test_v1_TDR"
#test_dir="supercell_100m_10s_0.05kmAGL_Cband_2.0BW_x215_y375_4490-4580s_FVEL_3.0deg"
#test_dir="supercell_100m_10s_0.1kmAGL_Cband_2.0BW_x215_y655_4490-5090s_FVEL"
#test_dir="supercell_100m_10s_0.1kmAGL_Cband_2.0BW_x215_y315_4490-5090s_FVEL"
test_dir="supercell_100m_10s_0.5kmAGL_Cband_2.0BW_x215_y375_4490-5090s_FVEL_test_v2"
#test_dir="supercell_500m_10s_0.5kmAGL_Cband_2.0BW_x066_y100_4490-5390s_FVEL"
#test_dir="supercell_500m_10s_3.0kmAGL_Cband_2.0BW_x066_y150_4490-5090s_FVEL_Aft360_test_v1_TDR"
#test_dir="squall_line_500m_10s_1.0kmAGL_Cband_2.0BW_x475_y175_15990-17250s_FVEL_v2"
#test_dir="squall_line_100m_10s_1.0kmAGL_Cband_2.0BW_x1502_y1020_15096-15156s_DVtest"
#test_dir="squall_line_100m_10s_1.0kmAGL_Cband_2.0BW_x1452_y960_15030-15210s_FVEL_NewPHIdp"
#test_dir="tropical_cyclone_500m_10s_3kmAGL_x415_y610_216990-217890s_FVEL_Aft360_test_v4"
#test_dir="tropical_cyclone_500m_10s_3kmAGL_x775_y676_h250deg_216990-218620s_FVEL_Fore360_test_v1"
#test_dir="tropical_cyclone_500m_10s_3kmAGL_x775_y610_h270deg_216990-218620s_FVEL_Aft360_test_v1"
#test_dir="tropical_cyclone_500m_10s_3kmAGL_x415_y610_216990-218620s_FVEL_Fore360_test_v1_TDR"
#test_dir="tropical_cyclone_500m_10s_3kmAGL_x415_y610_216990-218620s_FVEL_TOP_updated_bw3.0deg_v2"
#test_dir="tropical_cyclone_500m_10s_3kmAGL_x415_y610_216990-218620s_FVEL_PRT_updated_v3"
#test_dir="tropical_cyclone_500m_10s_3kmAGL_x600_y815_216990-218620s_FVEL_BOT_updated_v1"
#test_dir="tropical_cyclone_500m_10s_3kmAGL_x458_y610_217170-217770s_FVEL_2panel_sur_test"
#test_dir="tropical_cyclone_500m_10s_3kmAGL_x430_y610_217050-218490s_FVEL_2panel_sur_test"
#test_dir="tropical_cyclone_500m_10s_3kmAGL_x600_y690_216990-218790s_FVEL_circumnav_TOP_updated"
#test_dir="hnr1_tropical_cyclone_1km_360s_3kmAGL_x45_y240_497160-505960s_v6_STB"
#test_dir="pecan_squall_line_3000m_600s_3.0kmAGL_Cband_2.0BW_x227_y138_135720-139045s_FVEL_MP08"

#test_dir="AK_supercell_100m_10s_1kmAGL_Cband_x425_y435_4680-4920s"
#test_dir="AK_supercell_100m_10s_1kmAGL_Cband_x425_y375_4680-4920s"
#test_dir="AK_squall_line_500M_10s_1kmAGL_Cband_x475_y175_15990-17250s"
#test_dir="AK_tropical_cyclone_500M_10s_3kmAGL_Cband_x530_y610_217470-217710s"

#### No user input needed beyond this point ####

# Output directory
# First check that the experiment directory exists (if it is not an
# empty string)

if [ ! -z "$experiment_dir" ]; then
   
   if [ ! -d ${base_dir}/${experiment_dir} ]; then
      mkdir -p ${base_dir}/${experiment_dir}
   fi

   output_dir=${base_dir}/${experiment_dir}/${test_dir}
else
   output_dir=${base_dir}/${test_dir}
fi

# Get your current working directory
current_dir=$(pwd)

echo "Running the AOSPRE code..."
echo "Your current working directory is: " ${current_dir}
echo "Your operating directory is: " ${operating_dir}
echo "Your scan mode is: " ${scan_mode}
echo "Your namelist file is: " ${namelist_file}
echo "Your output directory is: " ${output_dir}

##### END USER VARIABLES #####

##### EXECUTION OF AOSPRE CODE #####

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
echo "All specified directories are verified...Proceeding with AOSPRE execution..."

if [ ! -d ${output_dir} ]; then
	   mkdir ${output_dir}
fi

# Go to the operating directory
echo ""
echo "The current working directory is now the operating directory: " ${operating_dir}
cd ${operating_dir}

echo ""
echo "Command line entries: (par1, flighttime)   $1 $2"

# Get the starting date and time of the AOSPRE run
start_time=`date +%s.%N`

# Pull some information out of the namelist file?
# We need to know whether we will need to operate in parallel and the time designated
# for the flight

#fp_line=`grep -nw fpar ${namelist_file}`
#echo ${fp_line}
numprocs=`nproc`
nphalf=$(( $numprocs / 2 ))

if [ ${exp_res} -eq 500 ]; then
    max_procs=19;
else
    max_procs=7;
fi

#echo "Do you need to run the AOSPRE in parallel? (1=yes,0=no): "
par1=$1

if [ $par1 -eq 0 ]; then
   npuse=1;
   echo ""
   echo "You chose to run in serial mode, so only 1 processor will be used..."
   echo ""

else
   echo ""
   echo "You chose to run in parallel mode...determining the suggested number of processors..."
   #echo "Enter the length of your flight in seconds: "
   flighttime=$2

   echo ""
   echo "The current breakdown time window is set to 30 seconds for each processor..."
   echo "You have ${numprocs} CPUs on this machine...but the recommendation is to use no more than half of them..."
   echo "This means you should use no more than ${nphalf} processors..."
   echo "The maximum allowed for this resolution is ${max_procs}...this will be the number of processors allowed..."

   # Compute the expected number of processors needed based on time
   nwin=30   # seconds for each flight window
   npuse=$(( $flighttime / $nwin ))
   
   if [ $npuse -gt $nphalf ]; then
      #echo "The maximum number of processors was exceeded...using the 60 second time window..."
      #nwin=$(( $nwin * 2 ))
      #npuse=$(( $flighttime / $nwin ))
      
      if [ $nphalf -gt ${max_procs} ]; then
         nwin=$(( $flighttime / ${max_procs} ))
         echo "The maximum number of processors was exceeded...using ${nwin} time window..."
         npuse=${max_procs}
      fi
         
   fi
   
   echo "Given your flight time specified, you should use ${npuse} processors with $nwin second windows"
   echo ""
   echo "Moving forward with ${npuse} processors and $nwin second windows"
   
fi
echo ""
#echo "Sleeping..."
#sleep 20


# Check to make sure mpirun exists, if it doesn't, then you have to run using the traditional commands
# Assign the check variable to 0, if the command isn't found, set it to 1
runver=0
if ! command -v mpirun &> /dev/null
then
   echo "mpirun does not exist on this computer...running AOSPRE in the traditional manner..."
   echo ""
   runver=1
fi

# Execute the code
echo ""
echo "Initiating AOSPRE code..."
echo ""
echo "And off we go..."
echo ""
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "%%      AOSPRE SPECIFIC OUTPUT BEGINS       %%"
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

if [ ${runver} -eq 1 ]; then
#   ./../embed-crsim/a.out ${namelist_file}
   ./../embed-crsim_standard/a.out ${namelist_file}
else
   mpirun --mca bil openid,self,vader -np ${npuse} -report-child-jobs-separately -tag-output ../embed-crsim_openmpi/atest.out ${namelist_file}
fi

echo ""
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "%%       AOSPRE SPECIFIC OUTPUT ENDS        %%"
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

# Get the end date and time of the AOSPRE run
end_time=`date +%s.%N`

runtime=$( echo "$end_time - $start_time" | bc -l )

echo ""
echo "AOSPRE processing time = ${runtime} seconds"

# Create a new flightpath file from the multiple thread file version
numfp=`ls -1 flightpath_* | wc -l`
numfp2=$(expr ${numfp} - 1)
echo $numfp2

for (( i=0 ; i<=$numfp2 ; i++ ));
do
    echo $i
    if [ $i -eq 0 ]; then
       if [ -f "flightpath.ascii" ]; then
          rm flightpath.ascii
       fi
       cat flightpath_${i}.ascii > flightpath.ascii
    else
       cat flightpath_${i}.ascii >> flightpath.ascii
    fi
    
    rm flightpath_${i}.ascii
done

# Once complete, move created files to output directory
echo ""
echo "Processing complete...moving files to output directory..."
if [ ${scan_mode} == 'rhi' ]; then

   # Creating a directory within your output for the particular scan mode
   if [ ! -d ${output_dir}/${scan_mode} ]; then
       mkdir ${output_dir}/${scan_mode}
   fi

   mv RHI.*.nc ${output_dir}/${scan_mode}
   cp 'flightpath.ascii' ${output_dir}/${scan_mode}
   cp ${namelist_file} ${output_dir}/${scan_mode}
else

   # Creating a directory within your output for the particular scan mode
   if [ ! -d ${output_dir}/${scan_mode} ]; then
       mkdir ${output_dir}/${scan_mode}
   fi

   mv surveillance_*.nc ${output_dir}/${scan_mode}
   cp 'flightpath.ascii' ${output_dir}/${scan_mode}
   cp ${namelist_file} ${output_dir}/${scan_mode}
fi

# Return to the original directory when you started
echo ""
echo "Returning to your original working directory: " ${current_dir}
cd ${current_dir}

echo ""
echo "AOSPRE code execution is complete...Done!"


