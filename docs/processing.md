#Processing
The AOS code execution has several phases, including the initial
determination of model output file availability, radar volume
allocation based on scan parameters provided, model file opening and
reading, looping through the scan beam and gate information to fill
the volume for the specified variables, and outputting the data into
the common radar output format (CfRadial, Dixon and Lee 2016). This
process is completed over the number of model files that exist within
the time frame configured for the flight. In most cases, AOS will
interpolate between model output times, given that the aircraft
currently used in the simulation (the NCAR C-130) has an airspeed of
120 m s<sup>-1</sup> and can produce output about every ~2.0-2.3 s.
