from .utilities import *
from .spectrum import Spectrum
import numpy as np


class SpectrumPair(object):
    def __init__(self, primary, secondary, name=None):
        if name is not None:
            self.name = name
        else:
            self.name = randomString()
        if isinstance(primary, Spectrum):
            self.primary = Spectrum(self.name+', '+primary.name)
            self.primary.set_data_from_spectrum(primary.get_data())
            self.min, self.max = primary.get_data().min(0)['wavenumber'], primary.get_data().max(0)['wavenumber']
        else:
            raise TypeError('primary must be Spectrum object')
        if isinstance(secondary, Spectrum):
            self.secondary = Spectrum(self.name + ', ' + secondary.name)
            self.secondary.set_data_from_spectrum(secondary.get_data())
            self.secondary.par_data = secondary.par_data
        else:
            raise TypeError('secondary must be Spectrum object')

        self.reconciled = False
        self.combined = None
        self.discardedSecondaryUpper = None
        self.discardedSecondaryLower = None

    def remove_secondary(self, scale = 1.0, peak_width = 0.04, lower_intensity_threshold=0.0, in_place=True, **kwargs):
        print(dt.now().time(), 'normalizing...')
        self.primary.normalize_spectrum()
        print(dt.now().time(), 'trimming...')
        self.secondary.trim(self.min, self.max)
        self.secondary.normalize_spectrum()
        print(dt.now().time(), 'finding peaks...')
        data = self.primary.get_data()
        allSecPeaks = self.secondary.get_data()
        secPeaks = allSecPeaks[(allSecPeaks.intensity > lower_intensity_threshold) ]
        secPeaks = secPeaks[['wavenumber','intensity']]
        print('peaks detected: ', len(secPeaks))
        print(dt.now().time(), 'removing peaks...')
        std = peak_width/2
        base = self.primary.determine_baseline()
        baseline = base.mean()
        basestd=base.std()
        noise = pd.Series(np.random.normal(baseline,basestd/4),data.index)
        arrPeaks = list(np.asarray(secPeaks).transpose((1, 0)))
        three_std=peak_width*1.5
        print(dt.now().time(), 'preparing to generate gaussians...')
        gauss_peaks_locations = [data.copy()[(data.wavenumber >= wavenumber - three_std) & (
                    data.wavenumber <= wavenumber + three_std)] for wavenumber, intensity in zip(*arrPeaks)]
        print(dt.now().time(), 'generating gaussians...')
        gauss_peaks_intensity = [gauss.apply(df_gauss, axis=1, args=(std, intensity * scale, wavenumber))
                                 for gauss, wavenumber, intensity in zip(gauss_peaks_locations, *arrPeaks)]
        print(dt.now().time(), 'processing gaussians...')
        total_intensity = pd.concat(gauss_peaks_intensity).sort_index().groupby(level=0).sum()
        print(dt.now().time(), 'subtracting gaussians...')
        backup = data.intensity.copy()
        data.intensity -= total_intensity
        data.intensity = data.intensity.fillna(backup)
        data.intensity[data.intensity < noise] = noise
        tempSpectrum=Spectrum()
        tempSpectrum.set_data_from_spectrum(data)
        base = tempSpectrum.determine_baseline()
        baseline = base.mean()
        minbaseline = baseline.dropna().index.min()
        data.intensity[data.intensity < (noise-0.1)] = noise
        data.intensity[data.intensity > 1] = noise
        data = data[data.index > minbaseline]
        tempSpectrum.set_data_from_spectrum(data)
        if in_place:
            self.primary.set_data_from_spectrum(tempSpectrum.get_data())
            return self
        else:
            primary = Spectrum(self.name+'(SECONDARY REMOVED)')
            primary.set_data_from_spectrum(data)
            return SpectrumPair(primary, self.secondary,self.name+'(SECONDARY REMOVED)')

    def subtract_secondary_spectrum_with_matching_index(self, relative_scale=1, fill_negative=False, in_place=True):
        self.primary.normalize_spectrum()
        self.secondary.normalize_spectrum(0,relative_scale)
        combined = self.primary.get_data()
        secondary = self.secondary.get_data()
        combined = combined.set_index('wavenumber')
        secondary = secondary.set_index('wavenumber')
        tempSpectrum1 = Spectrum('TEMP')
        tempCombined = combined
        tempCombined['wavenumber']=tempCombined.index
        tempSpectrum1.set_data_from_spectrum(tempCombined)
        base = tempSpectrum1.determine_baseline()
        combined['intensity'] = combined['intensity']-secondary['intensity']
        combined['wavenumber'] =combined.index
        combined=combined.reindex()
        
        if fill_negative:
            
            baseline = base.mean()
            basestd=base.std()
            noise = pd.Series(np.random.normal(baseline,basestd/4),combined.index)
            minbaseline = baseline.dropna().index.min()
            combined.intensity[combined.intensity < noise] = noise
            combined = combined[combined.index > minbaseline]
        tempSpectrum = Spectrum(self.name + ' (SECONDARY REMOVED)')
        tempSpectrum.set_data_from_spectrum(combined)
        if in_place:
            self.primary.set_data_from_spectrum(tempSpectrum.get_data())
            return self
        else:
            return SpectrumPair(tempSpectrum, self.secondary,self.name+'(SECONDARY REMOVED)')


    def plot(self, num_points=None, show_baseline=False):
        if self.combined is not None:
            self.combined.plot(num_points, show_baseline)
        self.primary.plot(num_points, show_baseline)
        self.secondary.plot(num_points, show_baseline)

    def subtract_and_plot(self, **kwargs):
        self.remove_secondary(**kwargs)
        self.plot(**kwargs)

    def save_spectrum(self, file_path):
        self.primary.save_spectrum(file_path)

    def trim(self, min_, max_):
        self.primary.trim(min_,max_)
        self.secondary.trim(min_,max_)

    def normalize_spectra(self, min_=0,max_=1):
        self.primary.normalize_spectrum(min_,max_)
        self.secondary.normalize_spectrum(min_, max_)

    # def estimate_scale(self, wavelength_of_peak_to_use, tolerance=0.005, first_guess=1e6, max_iterations=30, peak_width=0.003, acceptable_difference=0.05):
    #     # TODO: fix formula that governs this
    #     testingSP = SpectrumPair(self.primary, self.secondary)
    #     tolPlusWidth=tolerance+peak_width
    #     testingSP.trim(wavelength_of_peak_to_use-tolPlusWidth*500, wavelength_of_peak_to_use+tolPlusWidth)
    #     testingSP.normalize_spectra(0,1)
    #     primaryData = testingSP.primary.get_data()
    #     secondaryData = testingSP.secondary.get_data()
    #     maxPrimaryData = max(primaryData[(primaryData['wavenumber']<=wavelength_of_peak_to_use+tolPlusWidth) & (primaryData['wavenumber']>=wavelength_of_peak_to_use-tolPlusWidth)]['intensity'].values)
    #     maxPrimaryIndex=primaryData[primaryData['intensity']==maxPrimaryData].index
    #     primaryIndex=maxPrimaryIndex.tolist()
    #     primaryIndex=primaryIndex[len(primaryIndex)//2]
    #     maxSecondaryData = max(secondaryData[(secondaryData['wavenumber']<=wavelength_of_peak_to_use+tolPlusWidth) & (secondaryData['wavenumber']>=wavelength_of_peak_to_use-tolPlusWidth)]['intensity'].values)
    #     maxSecondaryIndex = secondaryData[secondaryData['intensity'] == maxSecondaryData].index
    #     secondaryIndex = maxSecondaryIndex.tolist()
    #     secondaryIndex = secondaryIndex[len(secondaryIndex) // 2]
    #     guess=first_guess
    #     differences = []
    #     for i in range(max_iterations):
    #         print(dt.now().time(), f'testing guess {i} ({guess:e})...')
    #         newGuess, difference = self._estimate_scale_repeater(testingSP, guess, peak_width, primaryIndex, secondaryIndex)
    #         print(dt.now().time(), f'guess {i} difference: {difference:e}')
    #         guess=newGuess
    #         differences.append(difference)
    #         if abs(difference) < acceptable_difference:
    #             break
    #     return guess, differences

    # def _estimate_scale_repeater(self, spectrum_pair, guess, peak_width, primaryIndex, secondaryIndex):
    #     originalIntensity = spectrum_pair.primary.get_data().loc[primaryIndex].intensity
    #     secondaryIntensity = spectrum_pair.secondary.get_data().loc[secondaryIndex].intensity
    #     newPair=spectrum_pair.remove_secondary(guess, peak_width, in_place=False)
    #     newIntensity=newPair.primary.get_data().loc[primaryIndex].intensity
    #     newGuess = newIntensity/(guess*secondaryIntensity**2)
    #     actualDifference=newIntensity
    #     return newGuess, actualDifference