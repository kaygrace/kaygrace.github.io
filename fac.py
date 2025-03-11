import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
import lightkurve as lk
from Pyriod import Pyriod

import warnings
warnings.filterwarnings('ignore')

def false_prob_prewhitening(lc, pyriod, cycles):
    """Determines the false amplitude cutoff of a given periodogram when searching for pulsation periods.
    
    lc is the lightcurve you created earlier in your notebook
    pyriod is the pyriod object you created from the lightcurve
    cycles is the number of cycles you want the program to run for - 1000 cycles means a 1 in 1000 chance of something below the cutoff being a real frequency"""

    np.savetxt('tempresiduals.dat', np.vstack((pyriod.lc.time.value,np.array(pyriod.lc.resid.value))).T) # take the time and residuals columns from your lightcurve and save them in a table

    resid_array = pd.read_table("tempresiduals.dat", sep=r'\s+', header=None) # read in the table of time and residual flux created in the last line

    resid_array[1] = resid_array[1] + 1 # set the baseline for flux to center around 1 instead of 0

    max_amp_array = np.zeros(cycles) # an array to store each calcualted fac in
    
    time = np.zeros(len(resid_array[0])) # taking the time column from your residuals array
    flux = np.zeros(len(resid_array[1])) # taking the flux column from your residuals array

    for j in range(len(resid_array[0])): # cleaning up the arrays
        time[j] = str(resid_array[0][j]) # making time into a string
        element = str(resid_array[1][j]) # making residual flux into a string
        tempflux = element.split(' ') # getting rid of any units tacked onto residual flux

        # this piece may not be needed anymore but i haven't checked yet
        if tempflux[0] == '———': # deal with the --- that appears in TESS data sometimes
            tempflux[0] = np.nan # replace with nan

        flux[j] = tempflux[0] # write the residual flux value to our flux array after dealing with issues like units and ---

    for i in range(cycles): # this is where the monte carlo simulation of it all happens
        
        randlist = random.sample(range(0, len(time)), len(time)) # randomly sample the time list, and assign each time index to a new random index
        
        newflux = np.zeros(len(time)) # an array to store our new flux's in

        for k in range(len(time)): # assign each flux in the residual lightcurve a new index based on the randlist we created above
            newflux[k] = flux[randlist[k]]
            
        lc2 = lk.LightCurve(time = time, flux = newflux) # create a new lc with the randomized indices

        pg = lc2.to_periodogram(freq_unit = 'microHertz') # create a new periodogram from the new lc with the randomized indices
        
        max_amp_array[i] = pg.max_power*1000 # get our units correct for amplitude (make it into mma)
        
    return np.mean(max_amp_array) # return the false amplitude probability cutoff, which is the mean of the peak amplitude for however many cycles you asked for