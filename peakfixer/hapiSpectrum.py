from peakfixer import Spectrum
from .utilities import *
import radis
import pandas as pd
import numpy as np

class HapiSpectrum(Spectrum):
    def __init__(self, name=None):
        Spectrum.__init__(self,name)
        self.molec_name=None
        self.isotopologues=None
        self.simulated = None
    
    def simulate(self, molec_name, min_wavelength=0, max_wavelength=None, isotopologues=None, wavelengths=None, **kwargs):
        self.molec_name=molec_name
        self.isotopologues = isotopologues
        if wavelengths is not None:
            minWavelengthStep = abs(wavelengths.diff().min())*0.25
        else:
            minWavelengthStep = 0.01
        if self.isotopologues is not None:
            isotopes = ','.join([str(k) for k in self.isotopologues.keys()])
        else:
            isotopes = getAllIsotopes(self.molec_name)
        spec = radis.calc_spectrum(wavenum_min=min_wavelength, wavenum_max=max_wavelength, molecule=molec_name, isotope=isotopes, wstep=minWavelengthStep, databank='fetch', name=f'{self.name} (simulated)', **kwargs)
        spec.apply_slit(minWavelengthStep*4, 'cm-1')
        nu, coeff = spec.get('abscoeff', wunit='cm-1')
        self.simulated = pd.DataFrame(data = {'wavenumber':nu, 'intensity':coeff})

    def reconcile_wavelengths(self, wavelengths:pd.Series):
        blank_df = pd.DataFrame(data={'wavenumber':wavelengths})
        blank_df['intensity']=np.nan
        blank_df = blank_df.set_index('wavenumber')
        sim = self.simulated.set_index('wavenumber')
        joined = pd.concat((sim,blank_df))
        joined = joined[~joined.index.duplicated()]
        joined = joined.sort_index()
        interpolated=joined
        interpolated['intensity'] = joined['intensity'].interpolate(method='index')
        interpolated = interpolated[interpolated.index.isin(blank_df.index)]
        self.set_data_from_spectrum(pd.DataFrame(data = {'wavenumber':interpolated.index, 'intensity':interpolated['intensity']}))     

        
        
    
