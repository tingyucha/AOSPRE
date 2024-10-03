#!/bin/bash

# flight_planner_auto.sh
# This script allows the user to assign waypoints for flight planning in 
# the NCAR/EOL AOS software. Please use flight_config to set up your desired
# leg specifications. Instructions for different assignments are in that config
# file.

# Load the information from the config file
#source flight_config_auto $1 $2 $3 $4 $5 $6 $7 $8

STVALS=("${!1}")
STV_USE=(${STVALS[@]})

TOTAL_TIME=${STV_USE[0]}
START_TIME=${STV_USE[1]}
MODEL_DX=${STV_USE[2]}
START_X=${STV_USE[3]}
START_Y=${STV_USE[4]}
START_ALT=${STV_USE[5]}
AC_HEAD=("${!2}")
AC_HEAD_USE=(${AC_HEAD[@]})
LEG_TIME=("${!3}")
LEG_TIME_USE=(${LEG_TIME[@]})
AC_ALT=("${!4}")
AC_ALT_USE=(${AC_ALT[@]})
AC_SPEED=("${!5}")
AC_SPEED_USE=(${AC_SPEED[@]})

echo "Total Time = ${TOTAL_TIME} s"
echo "Start Time = ${START_TIME} s"
echo "Model DX = ${MODEL_DX}"
echo "Start X = ${START_X}"
echo "Start Y = ${START_Y}"
echo "Start Alt = ${START_ALT} m"
echo "Headings = ${AC_HEAD_USE[@]} deg"
echo "Leg Times = ${LEG_TIME_USE[@]} s"
echo "AC Alts = ${AC_ALT_USE[@]} m"
echo "AC Speeds = ${AC_SPEED_USE[@]} m/s"
echo ${AC_SPEED_USE[0]}

# Get number of elements in the pre-defined arrays - this corresonds
# to AC_HEAD
numel=${#AC_HEAD_USE[@]}
numel2=${#LEG_TIME_USE[@]}
numel3=${#AC_ALT_USE[@]}
numel4=${#AC_SPEED_USE[@]}
echo "Number of elements per array: ${numel}, ${numel2}, ${numel3}, ${numel4}"

# Check to see if all the parameters in the config file are valid.
# Invoke an error if something is incorrect

if [[ -z ${TOTAL_TIME} || -z ${MODEL_DX} || -z ${START_X} || \
     -z ${START_Y} || -z ${START_ALT} || -z {AC_HEAD_USE} || -z {LEG_TIME_USE} || \
     -z ${AC_ALT_USE} || -z {AC_SPEED_USE} ]]; then
     
     echo "One or more of the variables in flight_config is not set...Please check."
     echo "Exiting script..."
     exit
fi

if [[ ${numel} -ne ${numel2} || ${numel} -ne ${numel3} || ${numel} -ne ${numel4} ]]; then
    echo "Your array sizes for the AC_HEAD, LEG_TIMES, AC_SPEED, or AC_ALT are not equal...Please check."
    echo "Exiting script..."
    exit
fi

# Declare the arrays to store the flight information
declare -a XI_USE
declare -a YI_USE
declare -a ZI_USE
PI=3.141592654

for idx in "${!AC_HEAD_USE[@]}"; do
    lnum=$(($idx+1))
    echo "Getting waypoints for leg $lnum..."

    echo "Heading= ${AC_HEAD_USE[$idx]} deg., Duration= ${LEG_TIME_USE[$idx]} s"
    echo ""

    if [ $lnum -eq 1 ]; then
	   XI_USE[$idx]=$START_X
	   YI_USE[$idx]=$START_Y
       ZI_USE[$idx]=$START_ALT
	
       dist1=$(( ${LEG_TIME_USE[$idx]} * ${AC_SPEED_USE[$idx]} ))
	   echo "Leg ${lnum} distance = ${dist1} m..."

	   sin1=$( echo "scale=6;s(${AC_HEAD_USE[$idx]}*${PI}/180)" | bc -l )
	   cos1=$( echo "scale=6;c(${AC_HEAD_USE[$idx]}*${PI}/180)" | bc -l )

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
       ZI_USE[$lnum]=${AC_ALT_USE[$idx]}
       
       #echo "X array: ${XI_USE[*]}"
       #echo "Y array: ${YI_USE[*]}"
       #echo "Z array: ${ZI_USE[*]}"
       
       echo "Finished assigning waypoints for leg ${lnum}..."
       echo ""
       
    else
       dist2=$(( ${LEG_TIME_USE[$idx]} * ${AC_SPEED_USE[$idx]} ))
       echo "Leg ${lnum} distance = ${dist2} m..."

       sin2=$( echo "scale=6;s(${AC_HEAD_USE[$idx]}*${PI}/180)" | bc -l )
       cos2=$( echo "scale=6;c(${AC_HEAD_USE[$idx]}*${PI}/180)" | bc -l )

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
       ZI_USE[$lnum]=${AC_ALT_USE[$idx]}
       
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
          echo "Done processing waypoint planner...moving on..."
       else
          echo "Finished assigning waypoints for leg ${lnum}..."
       fi
  
    fi
    
    for ii in "${!XI_USE[@]}"; do
       if [ $ii -gt 0 ]; then
          newi=$(( $ii - 1 ))
          XIU[$newi]="${XI_USE[$ii]}, "
          YIU[$newi]="${YI_USE[$ii]}, "
          ZIU[$newi]="${ZI_USE[$ii]}, "
       fi
    done
done

