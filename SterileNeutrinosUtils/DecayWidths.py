import math

gFermi = 1.16637e-11 # 1/MeV^2
hbar = 6.582e-22 # MeV*seconds
c = 2.998e8 # m/seconds

eMass = 0.511 # MeV
muMass = 105.658 # MeV
piMass = 139.570 # MeV
piV = 0.97425
piF = 130 # MeV

def Length_Nu2MuPi(nuMass,lTheta2,nuEnergy):
    hMass = piMass
    lMass = muMass
    hV = piV
    hF = piF

    nuMomentum = math.sqrt(pow(nuEnergy,2)-pow(nuMass,2))
    gamma = nuEnergy/nuMass
    beta = nuMomentum/nuEnergy

    term1 = lTheta2*pow(gFermi,2)*pow(hV,2)*pow(hF,2)*pow(nuMass,3)
    term2 = pow((1 - (pow(lMass,2)/pow(nuMass,2))),2) - (pow(hMass,2)/pow(nuMass,2))*(1 + (pow(lMass,2)/pow(nuMass,2)))
    term3 = math.sqrt((1 - (pow(hMass-lMass,2)/pow(nuMass,2)))*(1 - (pow(hMass+lMass,2)/pow(nuMass,2))))
    width = term1*term2*term3

    lifetime = hbar/width
    length = beta*gamma*c*lifetime

    return length # m
