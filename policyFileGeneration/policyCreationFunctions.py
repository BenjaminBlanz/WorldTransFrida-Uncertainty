import pandas as pd
import numpy as np
from scipy.stats import qmc
import sys

################################################################################################
#######################           Sampling                ######################################
################################################################################################

def createSampling1(n_samples, minvalue, maxvalue):
    
    #sobol_engine = qmc.Sobol(d=1, scramble=True)
    #samples_cont = sobol_engine.random(n=n_samples)
    #values=samples_cont[:,0]*(maxvalue-minvalue) + minvalue

    return np.linspace(minvalue, maxvalue,n_samples)

def createSamplingYears(minvalue, maxvalue, step=5):
    
    #sobol_engine = qmc.Sobol(d=1, scramble=True)
    #samples_cont = sobol_engine.random(n=n_samples)
    
    choices = np.linspace(minvalue,maxvalue,int((maxvalue-minvalue)/step)+1)
    #indices = (samples_cont[:,0] * len(choices)).astype(int)
    #years = choices[indices.flatten()]
    
    return choices


def createSampling3(n_samples, startvalue, minvalue, maxvalue):
        
    sobol_engine = qmc.Sobol(d=3, scramble=True)
    samples_cont = sobol_engine.random(n=n_samples)
    
    p1vs=samples_cont[:,0]*(maxvalue-minvalue) + minvalue
    startYears = np.linspace(2025,2100,16)    
    indices = (samples_cont[:,1] * len(startYears)).astype(int)
    sYs = startYears[indices.flatten()]

    s1ds=np.zeros(n_samples)
    p1ds=np.zeros(n_samples)

    for i in range(n_samples):
        
        maxDuration=2150-sYs[i]-10
        choices = np.linspace(5, maxDuration, int((maxDuration-5)/5)+1)
        idx = (samples_cont[i,2] * len(choices)).astype(int)
        s1ds[i] = choices[idx]
            
        p1ds[i] = 2150 - sYs[i] - s1ds[i] 

    return np.zeros(n_samples)+startvalue, sYs, s1ds, p1ds, p1vs
    
    
def createSampling6(n_samples, startvalue, minvalue, maxvalue):

    sobol_engine = qmc.Sobol(d=6, scramble=True)
    samples_cont = sobol_engine.random(n=n_samples)
    

    p1vs=samples_cont[:,0]*(maxvalue-minvalue) + minvalue
    p2vs=samples_cont[:,1]*(maxvalue-minvalue) + minvalue
    
    startYears = np.linspace(2025,2100,16)    
    indices = (samples_cont[:,2] * len(startYears)).astype(int)
    sYs = startYears[indices.flatten()]
        
    s2ds=np.zeros(n_samples)
    p2ds=np.zeros(n_samples)
    s1ds=np.zeros(n_samples)
    p1ds=np.zeros(n_samples)

    
    for i in range(n_samples):
        
        maxDuration=2150-sYs[i]-25
        choices = np.linspace(5, maxDuration, int((maxDuration-5)/5)+1)
        idx = (samples_cont[i,3] * len(choices)).astype(int)
        s1ds[i] = choices[idx]
        
        maxDuration=2150-sYs[i]-s1ds[i]-15
        choices = np.linspace(10, maxDuration, int((maxDuration-10)/5)+1)
        idx = (samples_cont[i,4] * len(choices)).astype(int)
        p1ds[i] = choices[idx]
        
        maxDuration=2150-sYs[i]-s1ds[i]-p1ds[i] - 10
        choices = np.linspace(5, maxDuration, int((maxDuration-5)/5)+1)
        idx = (samples_cont[i,5] * len(choices)).astype(int)
        s2ds[i] = choices[idx]
        
        p2ds[i] = 2150 - sYs[i] - s1ds[i] - p1ds[i] - s2ds[i]  
    
    return np.zeros(n_samples)+startvalue, sYs, s1ds, p1ds, s2ds, p2ds, p1vs, p2vs

################################################################################################
##################      Graphicals creation               ######################################
################################################################################################


def createGraphicalFromSampling3(samples):
    sv=samples[0]
    v1=samples[-1]
    
    time = np.linspace(2020,2150, 131)
    
    init=np.linspace(sv, sv, int(samples[1])-2020)
    slope1=np.linspace(sv, v1, int(samples[2])+1)
    policy1=np.linspace(v1, v1,int(samples[3]))

    graphical=np.append(init, np.append(slope1, policy1))
    return time, graphical

def createGraphicalFromSampling6(samples):
    sv=samples[0]
    v1=samples[-2]
    v2=samples[-1]
    
    time = np.linspace(2020,2150, 131)
    
    init=np.linspace(sv, sv, int(samples[1])-2020)
    slope1=np.linspace(sv, v1, int(samples[2])+1)
    policy1=np.linspace(v1, v1,int(samples[3])-1)
    slope2=np.linspace(v1, v2, int(samples[4])+1)
    policy2=np.linspace(v2, v2, int(samples[5]))
    
    graphical=np.append(init, np.append(slope1, np.append(policy1, np.append(slope2, policy2))))
    return time, graphical


def calcDPSfromGraphical(startvalue, graphical, time=np.linspace(2020,2150,131), discountRate=0.04):
    DPS=0.0
    for i, year in enumerate(time): DPS+=(graphical[i]-startvalue)*(1.0-discountRate)**int(max(0,year-2025))
    return DPS


################################################################################################
################################################################################################
################################################################################################



def createSamplingOrderedByDPS(N, startvalue, minvalue, maxvalue, policyType='SixParm',
                               time=np.linspace(2020,2150,131), discountRate=0.04):
        
    if policyType=='ThreeParm':
        Samples = np.asarray(createSampling3(N, startvalue, minvalue, maxvalue))
    elif policyType=='SixParm':
        Samples = np.asarray(createSampling6(N, startvalue, minvalue, maxvalue))
    else:
        sys.exit('policyTypes other than ThreeParm and SixParm are not allowed:', policyType)

        
    SamplesPlus = np.append(Samples, np.zeros((1,N)), axis=0)
    Graphicals = np.zeros((len(time),N))
    
    for i in range(N):
        sample=Samples[:,i]
        if policyType=='ThreeParm': time, graphical = createGraphicalFromSampling3(sample)
        elif policyType=='SixParm': time, graphical = createGraphicalFromSampling6(sample)
        SamplesPlus[-1,i] = calcDPSfromGraphical(sample[0],graphical, discountRate=discountRate)
        Graphicals[:,i] = np.copy(graphical)
        
    sorted_idx = np.argsort(SamplesPlus[-1, :])
    sortedSamples = SamplesPlus[:, sorted_idx]
    sortedGraphicals = Graphicals[:,sorted_idx]
    
    return sortedSamples, sortedGraphicals, sorted_idx, time

    
def selectFromOrderedSampling(DPS, sortedSamples, sortedGraphicals, Nselect=100):
    
    target_values = np.linspace(DPS[0], DPS[-1], Nselect)

    indices = []
    for target in target_values:
        
        # This index is the first index where arr[index] >= target_value.
        search_idx = np.searchsorted(DPS, target)
        # To choose the closest actual array value, we check the neighbor if possible.        
        if search_idx > 0 and (search_idx == len(DPS) or np.abs(target - DPS[search_idx-1]) <= np.abs(target - DPS[search_idx])): 
            search_idx -= 1
            
        # Now check if this was already selected!
        if search_idx in indices:
            # Find upper and lower index that was not selected yet:
            i_up = search_idx
            i_lo = search_idx
            while i_up in indices: i_up+=1
            while i_lo in indices: i_lo+=-1
                
            if i_up >= len(DPS) or np.abs(target - DPS[i_up]) > np.abs(target - DPS[i_lo]):
                search_idx = i_lo
            else:
                search_idx = i_up
            
        indices.append(search_idx)        
    
    
    return DPS[indices], sortedSamples[:,indices], sortedGraphicals[:,indices], indices


def addSamplesToExistingFile(outfile, polIDs, polStrings, columns=['polID', 'policyString']):
    df = pd.read_csv(outfile)
                
    oldPolIDs = df.values[:,0].tolist()
    oldPolStrings = df.values[:,1].tolist()
        
    if len(list(dict.fromkeys(oldPolIDs))) != len(polIDs):
        sys.exit('Error: You want to write more or less policies than are currently in the file.\n'+
                 ' Outfile: '+outfile+' '+str(len(list(dict.fromkeys(oldPolIDs))))+' '+str(len(polIDs)))
    
    NPs=oldPolIDs.count(1)

    
    newPolIDs = []
    newPolStrings = []
    for i in range(len(polIDs)):
        for oldPolID in oldPolIDs[i*NPs:(i+1)*NPs]: newPolIDs.append(oldPolID)
        newPolIDs.append(polIDs[i])

        for oldPolString in oldPolStrings[i*NPs:(i+1)*NPs]: newPolStrings.append(oldPolString)
        newPolStrings.append(polStrings[i])    
    
    pd.DataFrame(zip(newPolIDs, newPolStrings)).to_csv(outfile, index=False, header=columns)
    return
