# History
# date, decsription, author, email, institution
# v1 29/07/2015  First functional program with output, Ian Ashton, I.G>C>Ashton@exeter.ac.uk, University of Exeter.
# v2 29/07/2015  Now only opens one netCDF file, Ian Ashton, I.G>C>Ashton@exeter.ac.uk, University of Exeter.
# v3 10/09/2015  Lots of changes! Peter Land, peland@pml.ac.uk
# v4 18/11/2015  Varies SeaCarb parameters k1k2, kf, ks and b, Peter Land
# v2.1 - 15/12/2016 Revert to simple run and include mask option and algorithm naming for output (for Pathfinders poject)

from __future__ import print_function
import rpy2.robjects as ro
import numpy as np,csv
import scipy as sp
import itertools, math, random
from netCDF4 import Dataset
from datetime import datetime as dt
import shutil
import traceback

# Some definitions
defaultN = 50
defaultRmses = {'ALK':2.5, 'DIC':2.5, 'pH':.005, 'pCO2insitu':2.} # nominal 'state-of-art' errors
defaultRmses2 = defaultRmses.copy()

writeNames = {'pH':'pH', 'ALK':'Total alkalinity', 'DIC':'Dissolved inorganic carbon', 'pCO2':'pCO2','OmegaAragonite':'OmegaAragonite','OmegaCalcite':'OmegaCalcite'}
inNames = ['pH', 'ALK', 'DIC', 'pCO2insitu', 'S', 'T', 'Patm', 'P', 'Pt',
   'Sit', 'k1k2', 'kf', 'ks', 'pHscale', 'b', 'gas'] # Names required for the SeaCarb program
defaults = {'S':35., 'T':25., 'Patm':1, 'P':0, 'Pt':0, 'Sit':0,
   'k1k2':'x', 'kf':'x', 'ks':'d', 'pHscale':'SWS', 'b':'u74',
   'gas':'potential'}
seacarbAlgorithms = {'k1k2':['m10', 'm06', 'l', 'r'], # ensemble info
   'kf':['dg', 'pf'], 'ks':['d', 'k'], 'b':['l10', 'u74']}
defaultAlgorithms = {'kf':'dg', 'ks':'d'}
Slimits = {'m10':[1.,50.], 'm06':[.1, 50.], 'l':[19., 43.], 'r':[5., 45.],
   'dg':[0., 45.], 'pf':[10., 40.], 'd':[5., 45.], 'k':[20., 45.]}
Tlimits = {'m10':[0.,50.], 'm06':[1., 25.], 'l':[2., 35.], 'r':[0., 45.],
   'dg':[0., 45.], 'pf':[9., 33.], 'd':[0., 45.], 'k':[5., 40.]}
carbNames = inNames[:4]
statsNames = carbNames + ['OmegaAragonite', 'OmegaCalcite']
# The seacarb documentation has a full list of outputs that can be added
# to this list to enable them as outputs. Each entry in the list refers to
# the corresponding parameter in inNames. This formulation allows more
# attributes to be specified if necessary, by adding more fields to the
# dictionary at the relevant location in the list.    
outNames = ['fCO2', 'CO2', 'ALK', 'P', 'fCO2pot', 'pCO2insitu', 'flag',
   'OmegaAragonite', 'fCO2insitu', 'S', 'pCO2', 'T', 'OmegaCalcite',
   'CO3', 'pH', 'DIC', 'Patm', 'HCO3', 'pCO2pot']
outputAttributes = [{'units':'uatm'}, {'units':'mol/kg'},
   {'units':'mol/kg'}, {'units':'bar'}, {'units':'uatm'},
   {'units':'uatm'}, {'units':'1'}, {'units':'1'},
   {'units':'uatm'}, {'units':'1'}, {'units':'uatm'},
   {'units':'K'}, {'units':'1'}, {'units':'mol/kg'},
   {'units':'1'}, {'units':'mol/kg'}, {'units':'bar'},
   {'units':'mol/kg'}, {'units':'uatm'}]
attributes = {'missing_value':-999}
#And now define the flag for each pair
pairlist = [[0,1], [0,2], [1,2], [3,0], [3,1], [3,2]] # Pairs of carbonate parameters accepted in SeaCarb
flaglist = [8, 9, 15, 21, 24, 25] # Flags associated with each pair
pairFlagList = zip(pairlist, flaglist)
npairs = len(pairlist)
idx = []
missing_value=-999

def SeaCarb_IGA(infilename, outfilename, finalMask, pH, ALK, DIC, PCO2, S=35., T=25.,
   Patm=1., P=0., Pt=0., Sit=0., k1k2='x', kf='x', ks='d', pHscale='T',
   b='u74', gas='potential'):
   '''Imports data from specified netCDF to be run through SeaCarb R carb function. 
   Tested on SeaCarb version 3.0.8, https://cran.r-project.org/web/packages/seacarb/index.html
   Inputs: filename for input data (must exist). Filename for output data (does not need to exist)
   The variable names for at least 2 of the carbonate system parameters to be used. Non 
   entries should be left blank, []. Other suitable inputs can be defined using a variable
   name and will be exctracted from the netCDF file (e.g. temperature). 
   Those properties not in the input netCDF file can be defined by a numerical input or 
   no input to revert to defaults.
   A full description of Seacarb, the inputs and default values can be found here: 
   https://cran.r-project.org/web/packages/seacarb/seacarb.pdf
   Ian Ashton, 22/07/2015'''

   Uinputs = [pH,ALK,DIC,PCO2,S,T,Patm, P, Pt, Sit, k1k2, kf, ks, pHscale, b, gas]
   #outNames = {}
   outData = {}
   units = {}
   #N and RMSE values here. Note that in runSeaCarbEnsemble, the RMSE value is later calculated as a proportion of the mean value.
   N = 50
   RMSE = 0.1
   for mn in range(1,13):
     month = "%02d" % mn
     print(month)
     (indata,in_list,nc_dims,formask,dtype,unit_change,lat,lon) = SeaCarb_inputs(infilename,outfilename,Uinputs,inNames,mn-1)
     if len(formask[~formask.mask])>0:
       #get all possible pairs from optional inputs
       pairs = SeaCarb_getpairs(in_list)
       for pair in pairs:
          tflag = 0 # match the current pair with the relevant flag value
          for i,flag in enumerate(flaglist):
              if pair==pairlist[i]:
                  tflag = flag
                  stflag = str(tflag)
                  break
          if tflag:#If this combination is one accepted by SeaCarb
              out,fails = runSeaCarb(indata,tflag, inNames[pair[0]], inNames[pair[1]],verb = True)#Run the R function
       ## Run the R program and allocate the output in NetCDF file (based on the original)
      #return(out)
          for outvarind,name in enumerate(out):#For each output parameter
            if mn == 1:
               if name in writeNames.keys():
                  writeNames[name] = name+'_SeaCarb_flagvalue_'+stflag#Provide a new name dependent on the flag value
               outData[name] = np.empty((12,180,360))
            tdata = np.copy(formask)#Copied to make a masked array of same size as inputs
            tdata[~formask.mask] = out[name]#As this is the same mask used for the inputs, the outputs should always fit into the non-masked array
            tdata[tdata<=-999] = missing_value#This masks any -9999's that were set in runSeaCarb when it failed
            if outputAttributes[outvarind]['units'] == 'mol/kg' and unit_change == 'umol/kg':
                units[name] = 'umol/kg'
                #p1data = p1data*1000000
                tdata[np.where(tdata!=missing_value)] = tdata[np.where(tdata!=missing_value)]*1000000
            else:
                units[name] = outputAttributes[outvarind]['units']
            outData[name][mn-1,:,:] = np.reshape(tdata,(1,180,360))#Assign data to new variable in netCDF file
     else:
       for outvarind,name in enumerate(out):
        outData[name][mn-1,:,:] = np.zeros((1,180,360))-999
        #outData[name][mn-1,:,:].mask = True

   with Dataset(outfilename,'w') as ncout:#Open netCDF file for reading
     ncout.createDimension('longitude',360)
     ncout.createDimension('latitude',180)
     ncout.createDimension('time',12)

     longOut = ncout.createVariable('longitude', 'f4', ('longitude',))
     longOut.units = 'degrees_east'
     longOut.axis = 'X'
     longOut.long_name = 'Longitude'
     longOut[:] = lon
     latOut = ncout.createVariable('latitude', 'f4', ('latitude',))
     latOut.units = 'degrees_north'
     latOut.axis = 'Y'
     latOut.long_name = 'Latitude'
     latOut[:] = lat

     for name in writeNames:#
       print('Writing ',writeNames[name])
       output = ncout.createVariable(writeNames[name], dtype,('time', 'latitude', 'longitude'), fill_value=-999.)#Create output variable for the new estimates
       output.missing_value = attributes['missing_value']
       output.units = units[name]
       toutData = outData[name]
       toutData[toutData<=0.05] = attributes['missing_value']
       toutData[finalMask] = attributes['missing_value']
       output[:] = toutData
       output.long_name = writeNames[name]
     setattr(ncout, 'Conventions', 'CF-1.6') 
     setattr(ncout, 'Data source','Created using SeaCarb R function, https://cran.r-project.org/web/packages/seacarb/seacarb.pdf')
     setattr(ncout, 'Institution', 'European Space Agency (ESA) Pathfinders Ocean Acidification - University of Exeter; Plymouth Marine Laboratory (PML) and the Institut Francais de Recherches pour l\'Exploitation de la Mer (IFREMER)') 
     setattr(ncout, 'contact', 'email: J.D.Shutler@exeter.ac.uk') 
     setattr(ncout, 'sst_product',T)
     setattr(ncout, 'salinity_product',S)
     setattr(ncout, 'Version','SeaCarb IGA v3.1')
  
   return(0)

def SeaCarb_inputs(infilename,outfilename,Uinputs,inNames,idx):
   #Arranges input data for SeaCarb
   #Outputs the data to be analysed as a dictionary, dimensions, mask and datatype for the netCDF variables (to be used when writing)
   print(idx)
   unit_change = []#Switch in case we need to change the units
   shutil.copy2(infilename,outfilename)#Create output file as a copy of input file
   indata = dict()
   netCDFvars = []
   with Dataset(outfilename) as ncout:#Open netCDF file for reading
        lat = ncout.variables['latitude'][:]
        lon = ncout.variables['longitude'][:]
        ## Get NetCDF inputs
        for ind,Uinput in enumerate(Uinputs[0:9]):
            if type(Uinput) is str:
                try:
                    var = ncout.variables[Uinput]
                    tdata = var[idx,:,:]#If the variable name is given, extract from netcdf
                    # if idx:
                    #     tdata = var[idx,:,:]#If the variable name is given, extract from netcdf
                    #     print('idx recognised')
                    # else:
                    #     tdata = var[:]
                    #     print('idx not recognised')
                    nc_dims = tdata.shape
                    #print(Uinput,'data is',type(tdata),'With size',tdata.shape)
                except:
                    raise NameError(Uinput, 'not found in', outfilename, 'Or issue with input index',idx)
                #Units check
                try: 
                  if var.units == 'umol/kg':
                    tdata /= 1000000
                    unit_change = 'umol/kg'
                except:
                    print('\nNo units\n')
                    if inNames[ind] is 'DIC' or inNames[ind] is 'ALK':
                      tdata /= 1000000
                      unit_change = 'umol/kg'#Units are assumed to need changing

                indata[inNames[ind]] = tdata
                if netCDFvars:#If this is not the first time we have had a netCDF we add to existing which effectively combines the mask of all netCDF inputs
                    formask = formask+tdata#Adding will mask any values masked in either data set, the values are not important. (Could be done with just the mask)
                else:
                    formask = tdata
                    dtype = ncout.variables[Uinput].datatype
                netCDFvars.append(ind)#creating a list of netCDF inputs
                
            else:
                indata[inNames[ind]] = Uinput#Otherwise, not netCDF so take the input (or default) value

   if len(formask[~formask.mask])==0:
      print('THERE ARE NO UNMASKED VALUES - NO DATA FOR ONE OR MORE INPUT, OR NO CO-LOCATED DATA ACROSS INPUTS')
      #NEEDS ERROR BREAK

   for v in netCDFvars:
      tdata = indata[inNames[v]][~formask.mask].data#For each netcdf variable, make them the values that are not masked (by any input variable)   
      #print('\nv',inNames[v],'And the data are',tdata) 
      indata[inNames[v]] = tdata.tolist()

   ## The non-NetCDF variables

   for i,Uinput in enumerate(Uinputs[9:]):
      indata[inNames[i+9]] = Uinput#take defaults or entered values for switches
 
   ## Establishing the chosen parameters (at least 2 out of 4)     
   in_list = [ind for ind,Uinput in enumerate(Uinputs[0:4]) if Uinput]
   if len(in_list)<2:#If there are less than 2 of the first 4 required inputs
      print('Requires two input variables')
      #NEEDS ERROR BREAK

   return indata,in_list,nc_dims,formask,dtype,unit_change,lat,lon

def SeaCarb_getpairs(bob):
   #Returns index values for every possible pair from selected inputs (out of 4)
   #i.e. [1,2,3] returns [1,2],[1,3],[2,3],[3,1],[3,2],[2,1]

   pairs = [list(x) for x in itertools.combinations(bob, 2)]#Gets every possible pair of the Carbonate parameters input to see if they are an accepted pair in SeaCarb
   pairs2 = pairs[:]
   for pair in pairs2:
      pairs.append([pair[1],pair[0]])#And their reverse, in case that is available
   return pairs

def runSeaCarb(indata, tflag, var1name, var2name, verb = False):
   '''Creates an R function to arrange data and run the carb function within the SeaCarb package.
   It assumes data are in a python dictionary, with the optional inputs defined by variable name
   At present this only takes the pHscale switch
   Ian Ashton, 22/07/2015
   Added k1k2, kf, ks and b switches
   Peter Land 18/11/15'''

   ro.r('library(seacarb)')
   #Perturb the data here (or make a sister function that perturbs the data to be called if necessary.
   ro.r('''
      f<-function(flags,val1,val2,S,T,Patm,P,Pt,Sit,k1k2,kf,ks,b,pHs){
         flags = as.numeric(matrix(data = flags, nrow = length(flags), ncol = 1))
         val1 = as.numeric(matrix(data = val1, nrow = length(val1), ncol = 1))
         val2 = as.numeric(matrix(data = val2, nrow = length(val2), ncol = 1))
         S = as.numeric(matrix(data = S, nrow = length(S), ncol = 1))
         T = as.numeric(matrix(data = T, nrow = length(T), ncol = 1))
         out = carb(flags,val1,val2,S,T,Patm,P,Pt,Sit,k1k2=k1k2,kf=kf,ks=ks,b=b,pHscale=pHs)
         return(out)
      }
      ''')

   rCarb = ro.r['f']#define rCarb as above function
   nd = len(indata[var1name])
   out = {}
   fails = []
   if 1:#try:
      res = rCarb([tflag]*nd, indata[var1name], indata[var2name],
         indata['S'], indata['T'], indata['Patm'], indata['P'],
         indata['Pt'], indata['Sit'], indata['k1k2'], indata['kf'],
         indata['ks'], indata['b'], indata['pHscale']) # Call R function
      for l,v in res.items(): # For each parameter in the output, convert into Python dictionary, out
         out[l] = []
         for j in range(len(v)):
            if not type(v[j]) in [float, int]:
               raise ValueError(res.items(),l,v,j,v[j])
            out[l].append(v[j])
   if 0:#except:
      out = {}
      nfails = 0
      for index in xrange(nd):
         try:
            res = rCarb(tflag, indata[var1name][index],
               indata[var2name][index], indata['S'][index],
               indata['T'][index], indata['Patm'], indata['P'],
               indata['Pt'], indata['Sit'], indata['k1k2'], indata['kf'],
               indata['ks'], indata['b'], indata['pHscale']) # Call R function
            for l,v in res.items(): # For each parameter in the output, convert into Python dictionary, out
               if l in out:
                  out[l].append(v[0])
               else:
                  out[l] = [v[0]]
               if not type(v[0]) in [float, int]:
                  raise ValueError(res.items(),l,v)
         except:
            #traceback.print_exc()
            nfails+=1
            fails.append(index)
            for key in out:
                out[key].append(-99999)

            print(nfails, 'fails out of', nd)
            # print('GOOD ->',indata[var1name][index-10],indata[var2name][index-10],indata['S'][index-10],indata['T'][index-10], indata['Patm'], indata['P'],indata['Pt'], indata['Sit'], indata['k1k2'], indata['kf'],indata['ks'], indata['b'], indata['pHscale'])
            # print('BAD ->',indata[var1name][index],indata[var2name][index],indata['S'][index],indata['T'][index], indata['Patm'], indata['P'],indata['Pt'], indata['Sit'], indata['k1k2'], indata['kf'],indata['ks'], indata['b'], indata['pHscale'])
         # if verb or nfails == nd:
         #    print(nfails, 'fails out of', nd)
   for l in out:
      out[l] = np.array(out[l])

   return(out,fails)

def isArray(item):
   '''checks whether item is a numpy array or masked array'''
   return type(item) in [np.ndarray, np.ma.MaskedArray, np.float32]

def inLimits(algorithm, S, T):
   return (S >= Slimits[algorithm][0] and S <= Slimits[algorithm][1] and
      T >= Tlimits[algorithm][0] and T <= Tlimits[algorithm][1])

