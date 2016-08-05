import math
import os,sys
import numpy as np
import SterileNeutrinosUtils.DecayWidths as DW

# Physics parameters
muMass = 105.658 # MeV
piMass = 139.570 # MeV
kMass = 493.677 # MeV

#Detector parameters
volumeSBND = 55.0 # m^3 [fiducial]
volumeUBOONE = 46.2 # m^3 [fiducial]
volumeICARUS = 259.6 # m^3 [fiducial]
distanceSBND = 110 # m
distanceUBOONE = 470 # m
distanceICARUS = 600 # m
fluxCorrFactorSBND = pow(distanceUBOONE,2)/float(pow(distanceSBND,2))
fluxCorrFactorUBOONE = pow(distanceUBOONE,2)/float(pow(distanceUBOONE,2))
fluxCorrFactorICARUS = pow(distanceUBOONE,2)/float(pow(distanceICARUS,2))


# Calculate the number of heavy sterile neutrino events in LArTPCs
def EventRate_Nu2MuPi(nuMass,nuTheta2,POT):
    # Store current working directory
    pwd = os.path.dirname(__file__)
    flux = np.genfromtxt(pwd+"/BnbNuFlux/flux_numu.dat",\
                         delimiter=" ",skip_header=1,names=True,dtype=None)

    for i in range(len(flux['Energy'])):
        if flux['Energy'][i]==0: flux['Flux'][i] = 0
        else: flux['Flux'][i] = flux['Flux'][i]

    rate = [[],[]] # Heavy sterile neutrinos rate spectrum as a function of energy
    for i in range(len(flux['Energy'])):
        # Decay length gets broken if it gets unphysical values (e.g. total energy less than mass)
        if nuMass<=piMass+muMass or nuMass+muMass>=kMass or flux['Energy'][i]<=nuMass: calcRate = 0
        else: calcRate = POT*nuTheta2*flux['Flux'][i]/DW.Length_Nu2MuPi(nuMass,nuTheta2,flux['Energy'][i])
        rate[0].append(flux['Energy'][i])
        rate[1].append(calcRate)

    totalRate = 0
    for i in range(len(rate[0])):
        # Get the bin width for integration, usually same size but you never know
        if i!=len(flux)-1: deltaE = flux['Energy'][i+1]-flux['Energy'][i]
        else: deltaE = flux['Energy'][i]-flux['Energy'][i-1]
        totalRate += rate[1][i]*deltaE

    totalRateSBND = totalRate*volumeSBND*fluxCorrFactorSBND
    totalRateUBOONE = totalRate*volumeUBOONE*fluxCorrFactorUBOONE
    totalRateICARUS = totalRate*volumeICARUS*fluxCorrFactorICARUS

    return totalRateSBND, totalRateUBOONE, totalRateICARUS
