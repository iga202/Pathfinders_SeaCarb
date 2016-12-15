# Pathfinders_SeaCarb
Code to implement SeaCarb using Python
Imports data from specified netCDF to be run through SeaCarb R carb function. 
   Tested on SeaCarb version 3.0.8, https://cran.r-project.org/web/packages/seacarb/index.html
   Inputs: filename for input data (must exist). Filename for output data (does not need to exist)
   The variable names for at least 2 of the carbonate system parameters to be used. Non 
   entries should be left blank, []. Other suitable inputs can be defined using a variable
   name and will be exctracted from the netCDF file (e.g. temperature). 
   Those properties not in the input netCDF file can be defined by a numerical input or 
   no input to revert to defaults.
   A full description of Seacarb, the inputs and default values can be found here: 
   https://cran.r-project.org/web/packages/seacarb/seacarb.pdf
   
   Ian Ashton, Dec 2015

The code should identify which inputs are provided and use these to calculate the remaining carbonate parameters. Outputs are in the form of a NEtCDF file that includes a flag value in the variable name. This flag value reflects those in the SeaCarb documentation and refers to the carbonate parameters used to calculate each output.
This has not been robustly tested and errors are likely to arise. Use with caution.
