for file in surveillance*.nc
do
	echo $file
	flen=${#file}
	echo $flen
	newlen=$(expr $flen - 3)
	echo $newlen
	newstr=`echo $file | cut -c 14-$newlen`
	echo $newstr
	newstr2=cfrad.${newstr}_APAR_simulated_by_WRF_AIR.nc
	echo $newstr2
	cp -r $file $newstr2
	#mv $file $newstr

done
