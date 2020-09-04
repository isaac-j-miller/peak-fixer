from peakfixer import *
import matplotlib.pyplot as plt

####PARAMETERS####
filename ='./unfixed/GA_0.6torr_72M_10X_600LPF_337.95K_02-06-19_13.58_Emission back parallel input_Si Bolometer [External Pos.5]_0.000960_80 kHz_6 micron Mylar_BEB_AVG186.csv'
min_wavelength=187.0  # cm-1
max_wavelength=215.0  # cm-1
medium='vacuum'  # other option is 'air'
gas_temperature=337.95  # K
water_partial_pressure=0.000000008 #I'm not actually sure if this is the water's partial pressure or the overall pressure, 
#but the peaks are WAY too broad if I use the overall pressure, so I assume it is the partial pressure
isotopologues = {
    1:0.9978, 
    2:0.002,
    3:0.0003,
    4:0.000149,
    5:0.0000022
}  #abundance of water isotopologues. 1,2,3,4,5 are the HITRAN isotopologue identifiers
relative_abundance = 1.0e8
####__________####

water = HapiSpectrum('water')
test = Spectrum('test')
test.select_data_from_file(filename,headers=False) 
water.simulate('H2O', min_wavelength=min_wavelength, max_wavelength=max_wavelength, medium=medium,
    wavelengths=test.get_data()['wavenumber'], 
    Tgas=gas_temperature, 
    pressure=water_partial_pressure, isotopologues=isotopologues)
test.trim(min_wavelength, max_wavelength)
water.reconcile_wavelengths(test.get_data()['wavenumber'])

water.trim(min_wavelength, max_wavelength)
water.plot()
test.plot()
pair = SpectrumPair(test, water, 'Combined Spectrum')
pair.subtract_secondary_spectrum_with_matching_index(relative_abundance,fill_negative=True)
pair.primary.plot()
plt.show()

