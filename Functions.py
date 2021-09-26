
from typing import SupportsComplex
import numpy as np
import scipy.constants as cst
import pandas as pd
import isa 
import glob
### Plank function in nu/cm^(-1)
def plank(wavelength, temperature, Flux = False, units = 'cm'):
    if units == 'm': 
        wavelength = wavelength
        k = 1
    elif units == 'cm':
        wavelength = wavelength/100
        k = 1/100
    elif units == 'um':
        k = 10**-6
        wavelength = wavelength*10**-6
    
    if Flux == True: ### Returns Flux in w/(m^2 cm^-1)
        pifac = cst.pi
    if Flux == False:### Returns Radiance in w/(m^2 cm^-1 ster)
        pifac = 1 
    c1 = 2*cst.h*cst.c**2*k
    c2 = cst.h*cst.c/cst.k

    return c1/wavelength**5*(1/(np.exp(c2/(wavelength*temperature))-1))

        
def plank_nu(nu,temperature,Flux = False, units = 'cm'):
    ###The factors of 100 come from the conversion to cm because the constants are given in si units
    if units == 'm': ####These are inverse units 
        nu = nu
        k = 1
    elif units == 'cm':
        nu = nu*100
        k = 100
    elif units == 'um':
        k = 10*6
        nu = nu*10**6

    c1 = 2*cst.h*cst.c**2/k**2
    c2 = cst.h*cst.c/cst.k
    
    if Flux == True: ### Returns Flux in w/(m^2 cm^-1)
        pifac = cst.pi
    if Flux == False:### Returns Radiance in w/(m^2 cm^-1 ster)
        pifac = 1 
    
    return pifac*c1*nu**3/(np.exp(c2*nu/temperature)-1)*10**6



### This Class Takes HITRAN Gas Data and packs it nicely, I might add to it in future The path is to the HITRAN DATA file 
### The Hitran Data is in a non default output

### Description of the output, 
### Rename The output file to 'somegas'.csv,eg. CO2.csv

### The field Separator is a comma 
### Line Endings  Unix / Linux / Mac OS X (LF)
### Parameters in order:
### molecule ID
### nu
### S
### A
### gamma air 
### n air 
### gamma self 

class gas():
    ###load Data 
    ### The Path is to the 'somegas'.csv file 
    def __init__(self,path,mr,ppm = 415.1):##Mr is relative atomic mass in g/mol ##ppm ref is that of CO2
        HITRAN_data = pd.read_csv(path)
        self.nu = HITRAN_data['nu']
        self.nu_unit = '$cm^-1$'
        self.absorbtion_coeff = HITRAN_data['sw'] ###
        self.absorbtion_coeff_unit = r'$cm^-1/molecule \cdot cm^-1' 
        self.mr = mr
        self.mr_unit = '$g/mol'
        self.data = HITRAN_data
        self.gamma_air = HITRAN_data['gamma_air']
        self.gamma_self = HITRAN_data['gamma_self']
        self.n_air = HITRAN_data['n_air']
        self.ppm = ppm ##Parts per million        
        self.pv =  ppm/10**6 ###Fractional content per volume 
        
        ### This Function returns the tau broadening for the Tau Broadening going through layers of atmosphere\ 
        ### This is quite computationally intensive  
    def Tau_for_whole_atmosphere(self,alt_1, alt_2,steps = 10):
        altitudes = np.linspace(alt_1,alt_2,steps)
        delta_alt = abs(altitudes[0]-altitudes[1])*100
        Tau = np.zeros_like(self.nu)
        for i , _ in enumerate(altitudes):
            u = isa.getDensity(altitudes[i])/1000*self.pv/(self.mr)*cst.N_A*delta_alt
            Tau += lorentzain_fit(self,altitudes[i],Dataframe= False)*u
        
        return Tau
### Specific integral for plank or similar function
def func_int(start,stop,f_x,step = 10):
    xs = np.linspace(start,stop,step) ###heights are the z co-ord in m
    integral = 0 
    for i in range(len(xs)):
        integral += f_x(xs[i])*(xs[1])
    return integral

class HITRAN_CONST(): ##Class to store the HITRAN Database reference values 
    def __init__(self):
        self.tref = 296
        self.pref = 1*10**5



###Lorentzian centred on x0 
def lorenzian(x,x0,gamma): ### x and gamma are 1d arrays or pd series of the same length 
    return(1/cst.pi*(gamma/((x-x0)**2+(gamma)**2)))


def main_peaks(gas,threshold = 0.99,start = 0, stop = None,Quantile = True): ##Function that \
    #Returns a 1D array of the  top percentile of peak indexes/
    # Use this to index the absorption coeff 
    from scipy.signal import find_peaks
    absorbtion_coeff = np.array(gas.absorbtion_coeff)[start:stop]
    # Section functionality enables selection of a small window of wavelengths
    peaks , _  = find_peaks(absorbtion_coeff)
    #The Scipy function detects turning points 
    if Quantile ==True :
        mean_height = np.quantile(absorbtion_coeff[peaks],threshold)
    else:
        mean_height = 0 
    #take top 1% of peaks
    peaks = np.array((absorbtion_coeff[peaks]>mean_height)*peaks)
    peaks = peaks[peaks != 0]
    return peaks

def gamma(gas,altitude): ### gamma broadening \ ### gas is gas class altitude is altitude in m 
    ### For both self and air broadening 
    T_ref = 296 ###K
    P_ref = 1 ##atm
    P = isa.getPressure(altitude,atm = True)
    T = isa.getTemperature(altitude)
    n_air = gas.n_air
    P_self = P*gas.pv*isa.getDensity(altitude)/isa.getDensity(0) ##This Line assumes that the gas concentration falls off 
    # at the same rate that air density does it calculates 
    tfrac = (T_ref/T)**n_air
    return tfrac*(gas.gamma_air*(P-P_self)+gas.gamma_self*P_self)

def lorentzain_fit(gas,altitude,threshold = 0.99 ,start = None,stop = None,Quantile = True, Dataframe = True): 
    peaks = main_peaks(gas,start= start,stop= stop,threshold= threshold,Quantile=True)#find peaks for gas 
    n_peak = np.shape(peaks)[0]#number of peaks
    if Dataframe == True: 
        spec_copy = np.array(gas.absorbtion_coeff.iloc[start:stop].copy())
    else:
        spec_copy = np.array(gas.absorbtion_coeff[start:stop].copy())
    spec_copy[peaks] = 0
    ###For loop applies a lorentzian fit 
    ### Applies the lorenzian profile of a central peak between n_peak_coverage on either side 
    ### Larger gamma values require larger values For n_peak_coverage, small gamma means that the profile is deltafunction like and doesn't effect
    ### further out peaks. 
    n_peak_coverage = 1 ###Th
    for i in range(0,(n_peak)):
        if i < n_peak_coverage :
            sl_ = 0
            slp = 0
        
        if i >= n_peak -n_peak_coverage:
            sl_ = int(peaks[i-n_peak_coverage])
            slp = int(peaks[-1])
        if i>= n_peak_coverage and i< n_peak-n_peak_coverage:
            sl_ =int(peaks[i-n_peak_coverage])
            slp = int(peaks[i+n_peak_coverage])
        if Dataframe == True:
            spec_copy[sl_:slp] += (gas.absorbtion_coeff.iloc[start:stop]).iloc[peaks[i]]\
                *lorenzian((gas.nu.iloc[start:stop]).iloc[sl_:slp]\
                    ,(gas.nu.iloc[start:stop]).iloc[peaks[i]],(gamma(gas,altitude).iloc[peaks[i]]))
        else:
            spec_copy[sl_:slp] += (gas.absorbtion_coeff[start:stop])[peaks[i]]\
            *lorenzian((gas.nu[start:stop])[sl_:slp]\
                ,(gas.nu[start:stop])[peaks[i]],(gamma(gas,altitude)[peaks[i]]))
    
    return spec_copy

### Fast integral for a 1d flux array with wavenumber or wavelength
def Flux_Int(flux,nu):
    spacing  = nu-np.roll(nu,1)
    spacing[0] = spacing[1]
    integral = np.sum(spacing*flux)
    return integral



### Returns the flux through a layer of atmosphere 1d array 
def delta_I(gas,incident_flux,alt_1,alt_2,steps, Ground  = True):# the radiation goes from alt 1 to  2 takes a list of gas classes
    tau = gas.Tau_for_whole_atmosphere(alt_1,alt_2,steps)
    if Ground == True:
        flux = np.array(np.exp(-tau)*incident_flux) + plank_nu(gas.nu,isa.getTemperature(alt_2),Flux = True)
    else :
        flux = np.array(np.exp(-tau)*incident_flux) 
    return flux

### Returns 3 lists of arrays; of outgoing flux between the two layers , wavenumber for each gas, Incoming fluxes 
def multigas(gases,alt_1,alt_2,steps = 10, Ground  = True ):
    Incident_fluxes = []
    nuspec = []
    Outgoing_flux = []
    for i in range(len(gases)):    
        Incident_fluxes.append(plank_nu(gases[i].nu,isa.getTemperature(alt_1),Flux = True))
        nuspec.append(gases[i].nu)
        Outgoing_flux.append(delta_I(gases[i],Incident_fluxes[i],alt_1,alt_2, steps  = steps, Ground  = Ground ))
    return Outgoing_flux , nuspec , Incident_fluxes





###
###
###
###
###
###
###
### The Gas with Cross section Class is to deal with the Ozone profiles that are of A different Format Their is a short tutorial notebook on these functions
### Download all files for one gas for the different temperatures into one folder. 
### Path leads to the folder.  


class gas_with_crossection(): ### Similar to the gas class but for the diffrent data type 
    def __init__(self,path,mr,ppm):
        self.mr = mr ### relative atomic mass g/mol
        self.path = path
        self.entries  = len(crossection(glob.glob(path+'*')[0])) ### number of points in array 
        self.nu = HITRAN_crossection_nu(glob.glob(path+'*')[0], self.entries) ###wavenumber dataframe 1d for the cross section 
        self.temps = list(Crossections_temp_dataframe(self.path).columns) ### This returns a list of the temperatures available 
        self.ppm = ppm ### parts per million of the gas
    ### Returns Crossection dataframe for a certain Temperature
    def corossection(self,temp): ### Returns a 1d dataframe with temperature closest too the temperature input
        Crossections_df = Crossections_temp_dataframe(self.path)
        headers= list(Crossections_df.columns)
        new_head = np.zeros(6)
        for i in range(6):
            new_head[i] = int(headers[i])
        new_head = np.array(new_head)
        temp = np.argmin(abs(temp-new_head))
        return Crossections_df.iloc[:,temp]

### Function reads the dataframe asides from the first row 
def crossection(path):
    crossection = pd.read_csv(path,delimiter = " ", skiprows =1, header= None )
    crossection = np.array(crossection.iloc[:,1:])
    crossection = crossection.flatten()
    return crossection[:-1]



### Extracts the maxima and minima from the string at the top of the crossection file 
def HITRAN_crossection_nu(path,entries): ###Returns 1D Array of NU values with the same number of entries as the 
    data = open(path)
    first_line = str(data.readlines(1)[0])
    i = 0
    stops = 0
    zeros = 0 
    jj = 0
    while zeros <1: 
        jj = jj+1
        zeros = first_line[:jj].count('O')
        
    while stops < 1:
        i = i+1
        stops = first_line[:i].count('.')
    
    j = i 
    while stops < 2:
        j = j+1
        stops = stops = first_line[:j].count('.')
        
    wave_start = float(first_line[i:j])
    wave_stop = float(first_line[jj:i])
    nu = np.linspace(wave_start,wave_stop,entries)
    return nu   


def Crossections_temp_dataframe(path):
    all_crossect_paths = glob.glob(path+'*')
    entries = len(crossection(glob.glob(path+'*')[0]))
    header  = []
    corossections  = np.zeros((entries,len(all_crossect_paths)))
    i = 0
    for path in all_crossect_paths:
        header.append(path[46:49])
        corossections[:,i] = crossection(path)
        i = i+1
    corossections = pd.DataFrame(corossections, columns = header)
    return(corossections)

###
###
###
###
###
### Rayleigh Scattering 
### Function implements https://doi.org/10.1175/1520-0426(1999)016<1854:ORODC>2.0.CO;2 relationship for rayleigh wavenumber 
def Reigligh_crossection(wavenumber): ###Note that this crossection is accurate for 0.25 to 1 *10^-6um#                                    ##1000-4000 cm^-1
    wavenumber = wavenumber/10000
    numerator = 1.0455996 -341.29061*(wavenumber)**2 -0.90230850*wavenumber**(-2)
    denominator = 1+0.0027059889*(wavenumber)**2- 85.968563*wavenumber**(-2)
    return numerator/denominator


def Rayligh_optical_depth(wavenumber):
    sigma = Reigligh_crossection(wavenumber)
    return 0.0021520*sigma


### Further Functions 

def Wavenumber_concantenate(gas,upper_lim):
    gas_nu = np.array(gas.nu).flatten()
    nu = np.linspace(gas_nu[-1],upper_lim,1000)
    nu = np.hstack([gas_nu,nu])
    return nu 