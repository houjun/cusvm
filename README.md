# CENTRAL UNITED STATES VELOCITY MODEL (VERSION 1.3) ![Model overview. The regional structure is constructed by extending the physical properties of most of the well-defined units to less well-constrained areas in the Mississippi Embayment and the CUS: Paleozoic, Crystalline Upper Crust, Lower Crust, Modified Lower Crust and Upper Mantle. The Mesozoic to Cenozoic materials are represented by several units within the Mississippi Embayment and a Vs30 based layer outside the embayment.](IMAGES/model.jpg)
This set of utilities is used to extract velocity profiles from the central U.S. velocity model (Ramirez-Guzman et al., 2012).

## ABSTRACT

We have developed a new three-dimensional seismic velocity model of the central United States (CUSVM) that includes the New Madrid Seismic Zone (NMSZ) and covers parts of Arkansas, Mississippi, Alabama, Illinois, Missouri, Kentucky, and Tennessee (Ramírez‐Guzmán et al, 2012). The model represents a compilation of decades of crustal research consisting of seismic, aeromagnetic, and gravity profiles; geologic mapping; geophysical and geological borehole logs; and inversions of the regional seismic properties. The density and P- and S-wave velocities are synthesized in a stand-alone spatial database that can be queried to generate the required input for numerical seismic-wave propagation simulations. The velocity model has been tested and calibrated by simulating ground motions of the 18 April 2008 Mw 5.4 Mt. Carmel, Illinois, earthquake and comparing the results with observed records within the model area (see associated publication). The selected stations in the comparisons reflect different geological site conditions and cover distances ranging from 10 to 430 km from the epicenter. The results, based on a qualitative and quantitative goodness-of-fit (GOF) characterization, indicate that both within and outside the Mississippi Embayment the CUSVM reasonably reproduces: (1) the body and surface- wave arrival times; and (2) the observed regional variations in ground-motion amplitude, cumulative energy, duration, and frequency content up to a frequency of 1.0 Hz.

## SUMMARY AND GENERAL REQUIREMENTS

This README explains how to query the CUS velocity model version 1.3 from the command line. Prior to running the code, you must download the [database](https://www.sciencebase.gov/catalog/item/60880867d34e0fbbe3ddba0a) (Boyd and Ramirez-Guzman, 2021) and place it in the DATABASE directory. A matlab gui is also provided (see [GUI/README.md](../GUI/README.md)) but is not required to query the model. 

The extension of the velocity model is

          LON_min         = -92.9
          LON_max         = -85.1
          LAT_min         =  33.1
          LAT_max         =  39.9
          elevation_max   =  777 meters
          depth_min       =  0
          depth_max       =  60e3 meters
	
If you request a value outside this region it will return the closest value within the region.

You will need the following:

           DATABASE :   INCLUDES ALL THE UNCONFORMITIES, HORIZONS AND SURFACES 
                        OF THE MODEL AS WELL AS THE VELOCITY DISTRIBUTIONS. THE
                        DATABASE CAN BE DOWNLOADED FROM 
                        https://doi.org/10.5066/P995PCQY 
                    
           GEODATAQ :   INCLUDES THE SOURCE AND MAKE FILES TO BUILD THE PROGRAM 
                        GEODATAQUERY.

           EXAMPLE  :   INCLUDES ONE EXAMPLE ON HOW TO QUERY THE DATABASE.

           GUI      :   A MATLAB GUI THAT CALLS THE SYSTEM AND PLOTS CROSS-
                        SECTIONS, MAP VIEWS AND POINT PROPERTY PROFILES.

    GEOGRAPHIC_INFO :   GEOGRAPHIC REFERENCES USED BY THE GUI.

       compile_copy :   EXECUTABLE USED TO BUILD geodataquery

## INSTALLATION

1. Clone the repository using git in a unix terminal or download the project 
   zip file and unzip it to GDDIR (directory in which resides the velocity model 
   utilities and database, for example, /home/user/cusvm). Note that you will
   also have to download the database files from [ScienceBase](https://www.sciencebase.gov/catalog/item/60880867d34e0fbbe3ddba0a) 
   and place in the GDDIR/DATABASE directory.
2. Within a unix terminal and the directory GDDIR, make the file compile_copy 
   executable and run it.

   chmod +x compile_copy
   ./compile_copy

## EXAMPLE, INPUT AND OUTPUT FILE FORMAT

There are several ways you can query the velocity model (lines, grids,
request depth to surfaces). The most useful for wave propagation
simulation grids is to request points. In order to inspect the model we
recommend using the [Matlab GUI](../GUI/README.md), but note that MATLAB is not required
to query the model. Elevations and depths are given in meters. The 
following examples may also be used to verify the code was installed
correctly.

        The command line for geodataquery is:

        Usage: geodataquery typeofquery velorveldepthorbdrortopo db in rout 
                  typeofquery:  (0) points or (1) lines or planes
                                    Points is only applicable to velorveldepthorbdrortopo 0, 1, 
                                    and 2
     velorveldepthorbdrortopo:  (0) velocity & density (elevation above sea level) 
                                (1) velocity & density (depth from surface) 
                                (2) velocity & density (depth from surface, 
                                      topography squeezed to mean sea level) 
                                      Squeezing scales the depth by:
                                      (50e3-topo[lon_point][lat_point])/(50e3-elevation_max)
                                (3) elevation topo (relative to datum, positive down)
                                (4) elevation bedrock (relative to datum, positive down)
                                (5) depth to bedrock
                           db:  path to a database
                           in:  input file, query specifications
                         rout:  output path (directory), 
                                output file named output.out

        INPUT:
        The typeofquery for points is (0). The input file format is

                NumberOfStations
                lon_point_1 lat_point_1 elevation_or_depth_point_1
                lon_point_2 lat_point_2 elevation_or_depth_point_2

        The second parameter, velorveldepthorbdrortopo, specifies, for the first three
        cases (0,1,2), whether the location of interest is given in elevation (0) or
        depth (1,2). For example, the query of two points at an elevation of 10 and 100 
        meters above sea level would have an input file containing:

                2
                -90.5 35.5 10
                -88.3 37.4 100

        The type of query for lines and planes is (1). In this case, one of four types 
        of input file is needed. For velorveldepthorbdrortopo (0,1,2) and a horizontal 
        plane, the input file has the form:

                hororvert   = 0
                originLon   = -91
                originLat   = 35
                originDepth = 0
                lengthLon   = 1
                lengthLat   = 1
                numLon      = 100
                numLat      = 100

        a vertical plane:

                hororvert   = 1
                originLon   = -91
                originLat   = 35
                originDepth = 0
                finalLon    = -90
                finalLat    = 35
                finalDepth  = 1e4
                numAlong    = 100
                numD        = 100

        velorveldepthorbdrortopo (3,4,5) and a horizontal plane:

                hororvert = 0
                originLon = -91
                originLat = 35
                lengthLon = 1
                lengthLat = 1
                numLon    = 100
                numLat    = 100

        a vertical plane:

                hororvert = 1
                originLon = -91
                originLat = 35
                finalLon  = 1
                finalLat  = 1
                numAlong  = 100

        In the EXAMPLE directory you should find several example files. 
       
        For example, the command line to run input_Points_example.in if geodataquery 
        or geodataquery.exe is in GDDIR

         ./geodataquery 0 0 GDDIR/DATABASE/ GDDIR/EXAMPLE/input_Points_example.in GDDIR/

        OUTPUT:
        A file named output.out (in GDDIR) is created. The output has the same
        information as the input and the values of vp, vs and rho in m/s and kg/m^3. 
        The last column indicates the unit where the point is located. The output.out 
        file for input_Points_example.in is 

        -90.500000 35.500000 10.000000 1550.292660 382.796149 1667.000000 3
        -88.300000 37.400000 100.000000 2657.534211 409.354537 2139.454134 7

        Key for units in output file:
                0  - Air
                1  - St Louis 3D model
                2  - Quaternary Embayment
                3  - Upper Tertiary
                4  - Claiborne
                5  - Paleocene
                6  - Cretaceous
                7  - Meso/cenozoic
                8  - Paleozoic
                9  - upper crust
                10 - lower crust
                11 - modified lower crust
                12 - upper mantle

        The damping (Qs/Qp) is left as a parameter to be chosen by the user. Brocher
        (2008) provides an example that's used in the San Francisco Bay Area velocity 
        model.

## REFERENCES
*Boyd, O.S., Ramirez-Guzman, L. (2021). Database for the Central United States Velocity Model, v1.3: U.S. Geological Survey data release, https://doi.org/10.5066/P995PCQY.*

*Leonardo Ramírez‐Guzmán, Oliver S. Boyd, Stephen Hartzell, Robert A. Williams (2012). Seismic Velocity Model of the Central United States (Version 1): Description and Simulation of the 18 April 2008 Mt. Carmel, Illinois, Earthquake: Bulletin of the Seismological Society of America, 102 (6), 2622–2645. doi: https://doi.org/10.1785/0120110303.*

*Thomas M. Brocher (2008). Compressional and Shear-Wave Velocity versus Depth Relations for Common Rock Types in Northern California: Bulletin of the Seismological Society of America, 98 (2), 950–968, doi: https://doi.org/10.1785/0120060403.*
