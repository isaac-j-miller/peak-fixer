from peakfixer import *
import matplotlib.pyplot as plt

water = HapiSpectrum('water')

# water.select_data_from_hitran('H2O',187, 215, isotopologues=('1','2','3'))


filename ='./unfixed/GA_0.6torr_72M_10X_600LPF_337.95K_02-06-19_13.58_Emission back parallel input_Si Bolometer [External Pos.5]_0.000960_80 kHz_6 micron Mylar_BEB_AVG186.csv'
test = Spectrum('test')
test.select_data_from_file(filename,headers=False)  # test sets its spectrum to the data in this path
water.simulate('H2O', min_wavelength=187.0, max_wavelength=215.0, wavelengths=test.get_data()['wavenumber'], Tgas=337.95, pressure=0.0008, isotopologues={
    1:0.9978, 
    2:0.002,
    3:0.0003,
    4:0.000149,
    5:0.0000022
})
test.trim(187, 215)
water.reconcile_wavelengths(test.get_data()['wavenumber'])

water.trim(187, 215)
water.plot()
test.plot()
pair = SpectrumPair(test, water, 'Combined Spectrum')
pair.subtractSecondarySpectrumWithMatchingIndex(1.0e5,fill_negative=True)
pair.primary.plot()
plt.show()

