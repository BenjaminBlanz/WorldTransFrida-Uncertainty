import numpy as np
from scipy.stats import qmc
import pandas as pd
import os, sys, re

from policyCreationFunctions import *


allowedPolicies = ['SixParm', 'ThreeParm', 'OneParm', 'DivBy5Parm', 'Ones', 'Zeros']


############### Change inputs here #################
discountRate=0.04

# Actually need ~1000 times more priors than selections to get a good spread
# 2**16 ~ 65,000, which is enough to get around 100 evenly distributed samples,
# this will take around 10-15 min on Levante interactive node.
# For testing set this to  something like 1024 and 8
# (the SobolSequence sampler complains if the number is not a power of 2)
Nprior=int(2**16) 
Nselect=128

###################################################


######## Reading and ordering the policy lines from the policy info sheet #################
df = pd.read_csv('PolicyAnalysisVariables.csv', skip_blank_lines=True)
df = df.dropna(how='all')

# Check that every value in the 7th column is one of the allowed categories
if not df.iloc[:, 6].isin(allowedPolicies).all():
    # List rows with invalid category values for debugging
    invalid_rows = df.loc[~df.iloc[:, 6].isin(allowedPolicies)]
    sys.exit(f"Invalid policy type found in :\n{invalid_rows}")


policyList=df.values[:,0].tolist()
varList=df.values[:,2].tolist()
policyInputList=np.asarray(df.values[:,3:6])
policyTypes=df.values[:,6].tolist()
fullVarList=df.values[:,7].tolist()


### Ordering the input list by policyType
polTypeCounts = [policyTypes.count(s) for s in allowedPolicies]
writtenCounts = np.zeros(len(allowedPolicies)).astype(int).tolist()
newIndices = np.zeros(len(policyTypes)).astype(int).tolist()

for i, polType in enumerate(policyTypes):
    pidx = allowedPolicies.index(polType)
    newIdx = writtenCounts[pidx] + np.sum(polTypeCounts[:pidx]).astype(int)
    writtenCounts[pidx] += 1
    newIndices[newIdx] = i
    

policyList = [policyList[i] for i in newIndices]
varList = [varList[i] for i in newIndices]
fullVarList = [fullVarList[i] for i in newIndices]
policyInputList = [policyInputList[i] for i in newIndices]
policyTypes = [policyTypes[i] for i in newIndices]
############################################################################################




############### Create the new files #####################################

# Delete existing Folder?
if True:
    directory = "./PolicyFolder"
    if not os.path.exists(directory): os.makedirs(directory)
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        # Check if it's a file (and not a directory)
        if os.path.isfile(file_path):
            os.remove(file_path)

columns=['polID', 'policyString']
time=np.linspace(2020,2150,131)
timeFile=np.linspace(2020,2150,27)


print('Going through the policies')
for i, policy in enumerate(policyList):
                
    varname=fullVarList[i]
    
    outfile='PolicyFolder/'+policy+'.csv'
    if os.path.exists(outfile): fileExists=True
    else: fileExists=False

    polType=policyTypes[i]
    startvalue=policyInputList[i][0]
    minvalue=float(policyInputList[i][1].replace(",", ""))
    maxvalue=float(policyInputList[i][2].replace(",", ""))
    

    if polType in ['ThreeParm', 'SixParm']:
        sortedSamples, sortedGraphicals, sorted_idx, time = createSamplingOrderedByDPS(Nprior,
                                                                                       startvalue,
                                                                                       minvalue,
                                                                                       maxvalue,
                                                                                       policyType=polType,
                                                                                       time=time,
                                                                                       discountRate=discountRate)
        
        DPS=sortedSamples[-1]
        selected_DPS, selectedSamples, selectedGraphicals, indices = selectFromOrderedSampling(DPS, 
                                                                                      sortedSamples,
                                                                                      sortedGraphicals,
                                                                                      Nselect=Nselect)

        policyString_Xs = []
        policyString_Ys = []
        polIDs = []
        # File creation
        for iPol, sample in enumerate(np.transpose(selectedSamples)):
            polID=iPol+1


            graphical=selectedGraphicals[:,iPol]
            policyString_X = varname+':x'
            policyString_Y = varname+':y'
            for iy, year in enumerate(timeFile):
                policyString_X += ','+str(int(year))
                policyString_Y += ','+np.array2string(graphical[iy*5], formatter={'float_kind': lambda x: f'{x:.4e}'})

            policyString_Xs.append(policyString_X)
            policyString_Ys.append(policyString_Y)
            polIDs.append(polID)
            
        if fileExists:
            # Exit here because the sampling would be independent from the already existing sampling!
            sys.exit('Watch out: writing ThreeParm or SixParm PolType, but file already exists!')

            addSamplesToExistingFile(outfile, polIDs, policyString_Xs)
            addSamplesToExistingFile(outfile, polIDs, policyString_Ys)       
            
        else:
            df = pd.DataFrame(zip(polIDs, policyString_Xs), columns=columns)
            df.to_csv(outfile, index=False, header=True)
                        
            addSamplesToExistingFile(outfile, polIDs, policyString_Ys)


    elif polType in ['OneParm', 'DivBy5Parm']:
    
        if fileExists:
            # Exit here because the sampling would be independent from the already existing sampling!
            pass #sys.exit('Error: writing OneParm PolType, but file already exists: '+outfile)
                   
            # check how many variables are already in that file
            df = pd.read_csv(outfile)
            nvars=0
            while df.values[nvars,0]==1: nvars+=1
                
            bundlePolTypes = []
            bundleVarnames = []
            bundleMinValues = []
            bundleMaxValues = []
            # get the information on the variables that are already in there
            for i in range(nvars):                
                idx = fullVarList.index(re.split(r'[,:]', df.values[i,1])[0])
                bundlePolTypes.append(policyTypes[idx])
                bundleVarnames.append(fullVarList[idx])
                bundleMinValues.append(float(policyInputList[idx][1].replace(",", "")))
                bundleMaxValues.append(float(policyInputList[idx][2].replace(",", "")))
            
            bundlePolTypes.append(polType)
            bundleVarnames.append(varname)
            bundleMinValues.append(minvalue)
            bundleMaxValues.append(maxvalue)
            
            nVPlus=nvars+1
            
            # create the sampling
            sobol_engine = qmc.Sobol(d=nVPlus, scramble=True)
            samples_cont = sobol_engine.random(n=Nselect)
            
            for ii, polType in enumerate(bundlePolTypes):
                varname=bundleVarnames[ii]
                minvalue=bundleMinValues[ii]
                maxvalue=bundleMaxValues[ii]
                if polType=='OneParm':
                    selectedSamples = samples_cont[:,ii]*(maxvalue-minvalue) + minvalue
                elif polType=='DivBy5Parm':
                    choices = np.linspace(minvalue,maxvalue,int((maxvalue-minvalue)/5)+1)
                    indices = (samples_cont[:,ii] * len(choices)).astype(int)
                    selectedSamples = choices[indices.flatten()]
                    
                policyStrings = []
                polIDs = []
                for iPol, sample in enumerate(selectedSamples):
                    polID=iPol+1
                    if polType=='OneParm': policyString = varname+','+np.array2string(selectedSamples[iPol], formatter={'float_kind': lambda x: f'{x:.4e}'})
                    elif polType=='DivBy5Parm': policyString = varname+','+str(int(selectedSamples[iPol]))
                    else: sys.exit('Something went wrong!')
                    
                    policyStrings.append(policyString)
                    polIDs.append(polID)
            
                if ii == 0:
                    df = pd.DataFrame(zip(polIDs, policyStrings), columns=columns)
                    df.to_csv(outfile, index=False, header=True)
                    
                else:
                    addSamplesToExistingFile(outfile, polIDs, policyStrings)
                     
        else:
            if polType=='OneParm': selectedSamples = createSampling1(Nselect,minvalue,maxvalue)
            elif polType=='DivBy5Parm': selectedSamples = createSamplingYears(minvalue,maxvalue,step=5)
            policyStrings = []
            polIDs = []
            for iPol, sample in enumerate(selectedSamples):
                polID=iPol+1            
                if polType=='OneParm': policyString = varname+','+np.array2string(selectedSamples[iPol], formatter={'float_kind': lambda x: f'{x:.4e}'})
                elif polType=='DivBy5Parm': policyString = varname+','+str(int(selectedSamples[iPol]))
                    
                policyStrings.append(policyString)
                polIDs.append(polID)
            df = pd.DataFrame(zip(polIDs, policyStrings), columns=columns)
            df.to_csv(outfile, index=False, header=True)
                      
    elif polType in ['Ones', 'Zeros']:

        if fileExists:
            df = pd.read_csv(outfile)
            length=df.values[-1,0]
                
            policyStrings = []
            polIDs = []
            for iPol in range(length):
                polID=iPol+1            
                if polType=='Ones': policyString = varname+',1'
                elif polType=='Zeros': policyString = varname+',0'       
                    
                policyStrings.append(policyString)
                polIDs.append(polID)
                
            addSamplesToExistingFile(outfile, polIDs, policyStrings)
                     
        else:
            # Exit here because this should not happen, but theoretically works
            sys.exit('Error: writing Ones/Zeros PolType, but file does not yet exist!')
            
            if polType=='Ones': policyString = varname+',1'
            elif polType=='Zeros': policyString = varname+',0'
            
            df = pd.DataFrame(zip([1], [policyString]), columns=columns)
            df.to_csv(outfile, index=False, header=True)
        
    else: sys.exit('PolType not found: '+polType+' Variable = '+varname)
    print('... '+str(i+1)+' out of '+str(len(policyList))+' done!')

