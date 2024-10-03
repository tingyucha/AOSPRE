#!/bin/bash

# flight_planner.sh
# This script allows the user to assign waypoints for flight planning in 
# the NCAR/EOL AOS software. Please use flight_config to set up your desired
# leg specifications. Instructions for different assignments are in that config
# file.

# Load the information from the config file
source flight_config

# Get number of elements in the pre-defined arrays - this corresonds
# to AC_HEAD
numel=${#AC_HEAD[@]}
numel2=${#LEG_TIME[@]}
numel3=${#AC_ALT[@]}
#echo $numel

# Check to see if all the parameters in the config file are valid.
# Invoke an error if something is incorrect

if [[ -z ${TOTAL_TIME} || -z ${MODEL_DX} || -z ${START_X} || \
     -z ${START_Y} || -z ${START_ALT} || -z {AC_HEAD} || -z {LEG_TIME} || \
     -z ${AC_ALT} || -z {AC_SPEED} ]]; then
     
     echo "One or more of the variables in flight_config is not set...Please check."
     echo "Exiting script..."
     exit
fi

if [[ ${numel} -ne ${numel2} || ${numel} -ne ${numel3} ]]; then
    echo "Your array sizes for the AC_HEAD, LEG_TIMES, or AC_ALT are not equal...Please check."
    echo "Exiting script..."
    exit
fi

# Declare the arrays to store the flight information
declare -a XI_USE
declare -a YI_USE
declare -a ZI_USE
PI=3.141592654

for idx in "${!AC_HEAD[@]}"; do
    lnum=$(($idx+1))
    echo "Getting waypoints for leg $lnum..."

    echo "Heading= ${AC_HEAD[$idx]} deg., Duration= ${LEG_TIME[$idx]} s"
    echo ""

    if [ $lnum -eq 1 ]; then
	   XI_USE[$idx]=$START_X
	   YI_USE[$idx]=$START_Y
       ZI_USE[$idx]=$START_ALT
	
       dist1=$(( ${LEG_TIME[$idx]} * ${AC_SPEED} ))
	   echo "Leg ${lnum} distance = ${dist1} m..."

	   sin1=$( echo "scale=6;s(${AC_HEAD[$idx]}*${PI}/180)" | bc -l )
	   cos1=$( echo "scale=6;c(${AC_HEAD[$idx]}*${PI}/180)" | bc -l )

	   xdist=$( echo "${sin1} * ${dist1}" | bc -l )
	   ydist=$( echo "${cos1} * ${dist1}" | bc -l )
 
       xdindex=$( echo "${xdist} / ${MODEL_DX}" | bc -l )
       ydindex=$( echo "${ydist} / ${MODEL_DX}" | bc -l )

       #echo $xdist $ydist $xdindex $ydindex

	   XI_INT=$( echo "($xdindex+0.5)/1" | bc )
	   YI_INT=$( echo "($ydindex+0.5)/1" | bc )
	   XI_ADD=$(( ${XI_USE[$idx]} + $XI_INT ))
	   YI_ADD=$(( ${YI_USE[$idx]} + $YI_INT ))

	   #echo $XI_INT $YI_INT $XI_ADD $YI_ADD
 
       XI_USE[$lnum]=$XI_ADD
       YI_USE[$lnum]=$YI_ADD
       ZI_USE[$lnum]=${AC_ALT[$idx]}
       
       #echo "X array: ${XI_USE[*]}"
       #echo "Y array: ${YI_USE[*]}"
       #echo "Z array: ${ZI_USE[*]}"
       
       echo "Finished assigning waypoints for leg ${lnum}..."
       echo ""
       
    else
       dist2=$(( ${LEG_TIME[$idx]} * ${AC_SPEED} ))
       echo "Leg ${lnum} distance = ${dist2} m..."

       sin2=$( echo "scale=6;s(${AC_HEAD[$idx]}*${PI}/180)" | bc -l )
       cos2=$( echo "scale=6;c(${AC_HEAD[$idx]}*${PI}/180)" | bc -l )

       xdist=$( echo "${sin2} * ${dist2}" | bc -l )
       ydist=$( echo "${cos2} * ${dist2}" | bc -l )
 
       xdindex=$( echo "${xdist} / ${MODEL_DX}" | bc -l )
       ydindex=$( echo "${ydist} / ${MODEL_DX}" | bc -l )

       #echo $xdist $ydist $xdindex $ydindex

       XI_INT=$( echo "($xdindex+0.5)/1" | bc )
       YI_INT=$( echo "($ydindex+0.5)/1" | bc )
       XI_ADD=$(( ${XI_USE[$idx]} + $XI_INT ))
       YI_ADD=$(( ${YI_USE[$idx]} + $YI_INT ))

       #echo $XI_INT $YI_INT $XI_ADD $YI_ADD
 
       XI_USE[$lnum]=$XI_ADD
       YI_USE[$lnum]=$YI_ADD
       ZI_USE[$lnum]=${AC_ALT[$idx]}
       
       # Need to add a buffer at the end of the track to ensure
       # there are no errors. I have it set to 10 km, which is
       # equivalent to about 1.3 minutes of extra flight time.
       if [ $lnum -eq $numel ]; then
          track_buffer=10000
          xtb=$( echo "${sin2} * ${track_buffer}" | bc -l )
          ytb=$( echo "${cos2} * ${track_buffer}" | bc -l )
          xbind=$( echo "${xtb} / ${MODEL_DX}" | bc -l )
          ybind=$( echo "${ytb} / ${MODEL_DX}" | bc -l )
          XBI_INT=$( echo "($xbind+0.5)/1" | bc )
          YBI_INT=$( echo "($ybind+0.5)/1" | bc )
          
          XI_USE[$lnum]=$(( ${XI_USE[$lnum]} + ${XBI_INT} ))
          YI_USE[$lnum]=$(( ${YI_USE[$lnum]} + ${YBI_INT} ))
          
          echo ""
          echo "Finished processing all flight legs..."
          echo ""
          echo "The waypoints for your flight are: "
          echo "X waypoints: ${XI_USE[*]}"
          echo "Y waypoints: ${YI_USE[*]}"
          echo "Z waypoints: ${ZI_USE[*]}"
          
          echo ""
          echo "Done."
       else
          echo "Finished assigning waypoints for leg ${lnum}..."
       fi
       
       
    fi
done

