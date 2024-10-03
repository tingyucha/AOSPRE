#!/bin/bash

# run_flight_planner.sh
# This script is the parent script to flight_planner_auto.sh
# It asks the user to input some information regarding the
# AOS flight track, computes the track information, and passes
# the information to the namelist file of interest.
#
# This script assumes that you have the following directory
# structure set up:
# main directory --> model output directory
# main directory --> AOS directory
# main directory --> AOS directory --> code directory
# main directory --> AOS directory --> code directory --> AOS processing directory

echo "***********************************************************"
echo "** You are running run_flight_planner.sh. It is designed **"
echo "** to help you plan your flight path and copy it into a  **"
echo "** namelist file to be used with the AOS. This script    **"
echo "** does not run the AOS (yet) as its intent is to make   **"
echo "** flight track planning a smoother process.             **"
echo "**                                                       **"
echo "** IMPORTANT: Make sure you know which directory you are **"
echo "** initially storing your AOS output and create a file   **"
echo "** template that can be filled. Examples are located in  **"
echo "** the main code directory but may require changes based **"
echo "** on your model output naming and temporal output.      **"
echo "***********************************************************"
echo ""

sleep 3

echo "Let the flight planning begin..."
echo ""

# Set up your directory structure
maindir="/export/wind1/${USER}"
aosdir="${maindir}/AOS"
codedir="${aosdir}/code"

cd ${maindir}

echo "Now working from ${maindir}..."

# First, let's ask the user which storm they are interested in sampling
# Ask which directory they want to search in
# List the directories available
ls -1 -d */ > moddirs.txt
ii=0
declare -a modarr=()
while read -r line1; do
   ii=$(( $ii + 1 ))
   echo "${ii}) ${line1}";
   modarr+=(${line1})
done < moddirs.txt

echo "Please enter the directory number where your model ouput are stored: "
echo "If the directory isn't listed, please enter 999..."
read modeldir
echo ""
if (( ${modeldir} > 900 )); then
   echo "Please enter the directory name for your model output: "
   read modeldir2
else
   midx=$(( ${modeldir} - 1 ))
   modeldir2=${modarr[$midx]}
fi
full_modeldir=${maindir}/${modeldir2}
echo "Model directory set: ${modeldir2}"

# Next, ask the user where their namelist file and AOS output will be store

# Change to this directory
cd ${codedir}

echo "Now operating from ${codedir}..."

# List the directories available
ls -1 -d */ > dirlist.txt
ii=0
declare -a dirarr=()
while read -r line; do
   ii=$(( $ii + 1 ))
   echo "${ii}) ${line}";
   dirarr+=(${line})
done < dirlist.txt

echo "Please enter the directory number where your AOS namelist is stored: "
read processdir
didx=$(( ${processdir} - 1 ))
full_processdir=${codedir}/${dirarr[${didx}]}
echo "AOS processing directory set: ${dirarr[${didx}]}"

# Ask if the user is creating a new flight track or wants to use an existing one
echo "Are you creating a new flight track and namelist file (y/n)?: "
read newfile

if [ ${newfile} == "y" ]; then

# Now we need to get the track parameters of interest
# Determine the bounds of the available model grids
# The way this will work for now, is to use ncview to open
# the first, middle, and last model output file

   cd ${full_modeldir}

   echo "Now operating from ${full_modeldir}..."
   echo "The program ncview will now be opened for you to examine the model output"
   echo "and plan your flight. If you do not have ncview installed, exit this program"
   echo "and install it prior to running again."
   echo ""
   echo "If you have ncview installed, use it now to determine your starting location"
   echo "leg headings, and flight and leg durations. If you already have this information,"
   echo "you have the option to skip this step."
   echo ""

# Check to see if ncview is installed
   if ! command -v ncview &> /dev/null
   then
      echo "ncview is not installed on this machine. Please install to continue running...exiting program..."
      exit
   else
      echo "ncview is installed...moving forward with program..."
   fi

# Ask if the user wants to open ncview for flight planning
   echo "Do you need to open ncview for flight planning? (y/n) "
   read tf_ncv

   if [ ${tf_ncv} == "y" ]; then
      echo "You decided to open ncview. Once you are finished, please close and this"
      echo "program will continue."
      echo ""
   
   # Get a filelist in the modeldir and number of output files
      ls -1 wrfout* > filelist.txt
      numfiles=$( ls wrfout* | wc -l )
      halffiles=$(( $numfiles / 2 ))
      file1=$( sed -n 1p filelist.txt )
      file2=$( sed -n ${halffiles}p filelist.txt )
      file3=$( sed -n ${numfiles}p filelist.txt )
      echo $file1 $file2 $file3
   
   # Opening ncview with the 3 files listed
      ncview $file1 $file2 $file3
      echo ""
   else
      echo "You must already have your flight track information...moving on..."
      echo ""
   fi

   echo "Redirecting to the code directory"
   cd ${codedir}
   pwd

# Now, let's get the flight information you just gathered
   echo "We are now passing the flight track information to the flight_planner script."
   echo "You will need to provide flight path information in this order:"
   echo "Total flight time in seconds, start time in seconds, Model dx in meters, starting x, starting y, starting alt in meters,"
   echo "leg headings in degrees, leg times in seconds, leg altitudes in meters, and aircraft speed for each leg..."
   echo "You will be asked to specify each of these items in an upcoming set of commands..."
   echo ""
   echo "Note that the total time should be about 1 minute longer than you intend as a"
   echo "way to deal with a path error that can sometimes occur."
   echo ""

   echo "Please enter the total time, model dx, and starting position of x, y and alt:"
   echo "Example: 300 4490 500 275 315 1000"
   read startvals

   echo ""
   echo "Please enter the leg headings - the first entry is a dummy value,"
   echo "typically the same as your first leg heading."
   echo "Example: 90 90 45"
   read flhead

   echo ""
   echo "Please enter the leg time lengths in seconds - the first entry is a dummy value,"
   echo "and should be left at 0."
   echo "Example: 0 150 150"
   read fltime

   echo ""
   echo "Please enter the leg altitude in meters - the first entry is a dummy value,"
   echo "and should be set to your initial altitude."
   echo "Example: 1000 1000 1000"
   read flalt

   echo ""
   echo "Please enter the leg aircraft speeds in m/s - the first entry is a dummy value,"
   echo "and should be left at 120. Note that the C-130 should fly at speeds between 110-135 m/s."
   echo "Example: 120 120 120"
   read flspeeds

   . ./flight_planner_auto.sh "startvals[@]" "flhead[@]" "fltime[@]" "flalt[@]" "flspeeds[@]"

   #echo "Testing to see if variable was passed..."
   #echo "XI_USE = ${XI_USE[@]}"
   #echo "YI_USE = ${YI_USE[@]}"
   #echo "ZI_USE = ${ZI_USE[@]}"
   echo ""

# Change to processing directory
   echo "Changing to the AOS processing directory:"
   echo "${full_processdir}"

   cd ${full_processdir}

# Check to see if the template namelists exist
   nlf1="namelist_template_auto.rhi"
   nlf2="namelist_template_auto.surveillance"
   if [[ ! -f "$nlf1" || ! -f "$nlf2" ]]; then
      echo "The template namelist files do not exist...create them and try again..."
      echo "Exiting..."
      exit
   else
      echo "The template namelist files exist...continuing"
   fi
   echo ""

# Ask the user if they are running in parallel or serial mode
   echo "Will you process in parallel mode? (1=y, 0=n): "
   read parproc

# Ask the user whether this is surveillance or rhi scan mode
   echo "Are you running surveillance or RHI mode (Enter 'S' or 'R') ?"
   read scanmode

# Ask which panel the user wants to use
   if [ ${scanmode} == "R" ]; then
      echo "Which panel are you planning to use (LHS, RHS, TOP, BOT)?"
      read panelname
   
   # Copy the pertinent variables to the namelist file
      sed -e "s|FULLMD|${full_modeldir}|" \
      -e "s/PAR_PROC/${parproc}/" \
      -e "s/STIME/${START_TIME}/" \
      -e "s/TOTALTIME/${TOTAL_TIME}/" \
      -e "s/XIND/${XIU[*]}/" \
      -e "s/YIND/${YIU[*]}/" \
      -e "s/ZIND/${ZIU[*]}/" \
      -e "s/ACSPD/${AC_SPEED_USE[0]}/" \
      -e "s/PANEL/${panelname}/" \
      namelist_template_auto.rhi > namelist_filled_auto.rhi
   
   else
   # Copy the pertinent variables to the namelist file
      sed -e "s|FULLMD|${full_modeldir}|" \
      -e "s/PAR_PROC/${parproc}/" \
      -e "s/STIME/${START_TIME}/" \
      -e "s/TOTALTIME/${TOTAL_TIME}/" \
      -e "s/XIND/${XIU[*]}/" \
      -e "s/YIND/${YIU[*]}/" \
      -e "s/ZIND/${ZIU[*]}/" \
      -e "s/ACSPD/${AC_SPEED_USE[0]}/" \
      namelist_template_auto.surveillance > namelist_filled_auto.surveillance
   fi
   echo ""
   echo "Namelist file successfully populated..."
   echo ""
   
   echo "Would you like to save this namelist file for later use (y/n)?: "
   read savenlf
   echo ""
   
   if [ ${savenlf} == "y" ]; then
      echo "Enter the name of the namelist file you'd like to save: "
      echo "Note: this will be saved in your processing directory..."
      echo "Example: namelist_rhi_supercell_500m_WestEast_singleleg_10minute"
      read filesv
      
      if [ ${scanmode} == "R" ]; then
         cp namelist_filled_auto.rhi ${filesv}
      else
         cp namelist_filled_auto.surveillance ${filesv}
      fi
      echo ""
      echo "Namelist file saved..."
      echo ""
   fi
else
   echo "You chose to use an existing namelist file..."
   echo ""
   echo "Let's load the saved file and copy it to the generic auto file..."
   echo ""
   
   echo "Redirecting to the AOS processing directory..."
   cd ${full_processdir}
   
   # List the available saved files
   ls -1 namelist* > nlflist.txt
   kk=0
   declare -a nlfarr=()
   while read -r line; do
      kk=$(( $kk + 1 ))
      echo "${kk}) ${line}";
      nlfarr+=(${line})
   done < nlflist.txt

   echo "Please enter the file number you wish to use: "
   read nlfnum
   nidx=$(( ${nlfnum} - 1 ))
   full_nlname=${nlfarr[${nidx}]}
   
   if grep -q "rhi" <<< "${full_nlname}"; then
      scanmode="R"
      cp ${full_nlname} namelist_filled_auto.rhi
   else
      scanmode="S"
      cp ${full_nlname} namelist_filled_auto.surveillance
   fi
fi
echo "Returning to code directory..."
cd ${codedir}

echo "Now you can run the AOS with your namelist file"

echo "Exporting variables to use with RUN_AOS_auto.sh"
echo ""

echo CODEDIR=${codedir} >> ${aosdir}/AOS_env.txt
echo FULL_PD=${full_processdir} >> ${aosdir}/AOS_env.txt
echo SCANMODE=${scanmode} >> ${aosdir}/AOS_env.txt

echo "Finished running the flight planning script..."
echo "Exiting..."


