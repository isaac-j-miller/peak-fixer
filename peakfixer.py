import gc
import os
import random
import string
import math
from datetime import datetime as dt
# import copy
from tkinter import Tk, filedialog, messagebox

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

spectra = []

par_template={'molec_num':{'type':'int16',
                           'start': 0,
                           'width': 2},
              'iso_num':{'type':'int16',
                           'start': 2,
                           'width': 1},
              'wavenumber':{'type':'float64',
                           'start': 3,
                           'width': 12},
              'intensity':{'type':'float64',
                           'start': 15,
                           'width': 10},
              'einstein_a':{'type':'float16',
                           'start': 25,
                           'width': 10},
              'air_broad_w':{'type':'float16',
                           'start': 35,
                           'width': 5},
              'self_broad_w':{'type':'float16',
                           'start': 40,
                           'width': 5},
              'lower_state_e':{'type':'float16',
                           'start': 45,
                           'width': 10},
              't_dep_ir_w':{'type':'float16',
                           'start': 55,
                           'width': 4},
              'p_shift':{'type':'float16',
                           'start': 59,
                           'width': 8},
              'upper_v_q':{'type':'object',
                           'start': 67,
                           'width': 15},
              'lower_v_q':{'type':'object',
                           'start': 82,
                           'width': 15},
              'upper_l_q':{'type':'object',
                           'start': 97,
                           'width': 15},
              'lower_l_q':{'type':'object',
                           'start': 112,
                           'width': 15},
              'err_codes':{'type':'object',
                           'start': 127,
                           'width':6},
              'ref_codes':{'type':'object',
                           'start': 133,
                           'width':12},
              'flag_line_mix':{'type':'object',
                           'start': 145,
                           'width': 1},
              'upper_stat_wt':{'type':'float16',
                           'start': 146,
                           'width': 7},
              'lower_stat_wt':{'type':'float16',
                           'start': 153,
                           'width': 7},
              }

txt_types = {'.csv': ',',
             '.dpt': '\t'}


def randomString(stringLength=10):
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))


def fast_interleave_overlay(df1, df2):
    print(dt.now().time(),'interleaving...')
    # divisor = 1
    # done = False
    print(dt.now().time(),'assembling xlist and dropping duplicates')
    output = pd.DataFrame(columns=['wavenumber', 'intensity'])
    output['wavenumber'] = df1['wavenumber'].append(df2['wavenumber']).drop_duplicates(inplace=False).sort_values(kind = 'quicksort')
    print(dt.now().time(),'dropped',len(df1)+len(df2)-len(output),'duplicates')
    # df2.loc[:, 'intensity'] = np.nan
    output['wavenumber'] = np.nan
    # here, output is the common freq axis. df1 is the primary dataframe, df2 is a dataframe with blank values
    print(dt.now().time(),'adding secondary dataframe intensities...')
    output['intensity'] = df2[df2['wavenumber']==output['wavenumber']]['intensity']
    del df2
    print(dt.now().time(),'adding primary dataframe intensities...')
    output['intensity'] = df1[df1['wavenumber']==output['wavenumber']]['intensity']
    del df1
    return output
    # while not done:
    #     chunks1 = chunker(df1, divisor)
    #     chunks2 = chunker(df2, divisor)
    #     if len(chunks1) == len(chunks2):
    #         print(dt.now().time(),dt.now().time(),'divisor is ', divisor)
    #         try:
    #             inputs = [chunks1[i].append(chunks2[i]).drop_duplicates(inplace=False, subset = 'wavenumber') for i in range(len(chunks1))]
    #             name = saver(inputs)
    #             done = True
    #         except:
    #             print(dt.now().time(),dt.now().time(),'failure. increasing divisor')
    #             divisor += 1
    #     else:
    #         print(dt.now().time(),dt.now().time(),'unequal size')
    #         raise('SizeError: Unequal arrays')
    # del chunks1,chunks2,df1,df2

    # print(dt.now().time(),dt.now().time(),'sorting...')
    # return retrieve_from_saver(name)


def chunker(df,divisor):
    print(dt.now().time(),'chunking...')
    length = len(df)
    size = length // divisor
    return [df[i:i+size] for i in range(0,length-1,size)]


def saver(list_of_dfs):
    print(dt.now().time(),'saving temporary file...')
    fname = randomString() + '_temp.csv'
    temp_file = open(fname, 'w')
    temp_file.close()
    for item in list_of_dfs:
        item.to_csv(fname, mode='a', index_label=False, index=False, columns=['wavenumber', 'intensity'])
    del list_of_dfs
    return fname


def retrieve_from_saver(fname):
    print(dt.now().time(),'retrieving data...')
    data = pd.read_csv(fname, delimiter=',', header=0, names=['wavenumber', 'intensity'])
    print(dt.now().time(),'deleting temp file...')
    os.remove(fname)
    return data


def within_percent_of(x,y,percent):
    if abs((x-y)/np.average([x,y]))*100 <=percent:
        return True
    else:
        return False


def gauss(x, stdev, maximum, mean):
    return maximum*math.exp((-1*(x-mean)**2)/(2*stdev**2))


def df_gauss(row, stdev, maximum, mean):
    return gauss(row['wavenumber'], stdev, maximum, mean)


def specific_gauss(x, y, stdev, maximum, mean, x_min, x_max):
    if x_min < x < x_max:
        return gauss(x, stdev, maximum, mean)
    else:
        return y


def reject_greater(x, thresh):
    if x > thresh:
        return 0
    if x < thresh*-1:
        return 0
    return x


def reject_lesser(x,thresh):
    if x < thresh:
        return thresh
    return x


class Spectrum(object):
    def __init__(self, name = None):
        if name is not None and name not in spectra:
            self.name = name
        else:
            self.name = randomString(10)
        spectra.append(self)
        self._data = None
        self.cols = ['wavenumber', 'intensity']
        self._orig = None
        self.peak_list = None
        self.peak_tuples = None
        self.peak_tuples2 = None
        self.tolerance = .1
        self.min=0
        self.max=0

    def select_data_from_file(self, file_path=None, headers=None, silenced=False, overwrite=True, id_mods=None):
        print(dt.now().time(), 'loading...')
        if not silenced:
            root.deiconify()
        if file_path is None:
            file_path=filedialog.askopenfilename(title='Select Spectrum')
        index = file_path.rfind('.')
        filetype = file_path[index:].lower()

        if filetype in txt_types.keys():

            if headers is None:
                headers = messagebox.askyesnocancel('Headers?', 'Does this file have headers?')
            delimiter = txt_types[filetype]
            try:
                dtypes = {item: par_template[item]['type'] for item in self.cols}
                if headers:
                    file = pd.read_csv(file_path, delimiter=delimiter, header=0, names=self.cols)
                    self._data = file
                    del file
                    success = True
                elif not headers:
                    file = pd.read_csv(file_path, delimiter=delimiter,header = None, names=self.cols)
                    self._data = file
                    del file
                    success = True
                else:
                    success = None

            except:
                success = False

        elif filetype == '.par':
            #try:

            file = open(file_path,'r')
            data = file.readlines()
            file.close()
            to_pd = [{item: line[par_template[item]['start']:par_template[item]['start']+par_template[item]['width']] for item in par_template} for line in data]
            self.cols.append('iso_num')
            dtypes = {item: par_template[item]['type'] for item in self.cols}
            self._data = pd.DataFrame.from_dict(to_pd)[self.cols]
            test=self._data.intensity
            if id_mods is not None:
                isonum=self._data.iso_num.astype(int)
                self._data['mods'] = isonum.map(id_mods).astype('float64')
                self._data['intensity']=self._data.mods*self._data.intensity.astype(par_template['intensity']['type'])
            test2 = self._data.intensity
            self.cols.remove('iso_num')
            del to_pd, data, file
            success = True

            #except:
                #success = False
        else:
            success = False

        if success:
            self._data = self._data.astype(dtypes)
            # del dtypes, delimiter, file_path, index
            msg = 'The file "' + file_path + '" has been loaded into '+self.name+' spectrum.'
            if not silenced:
                messagebox.showinfo('Success!', msg)
            else:
                print(dt.now().time(),msg)
        elif success is not None:
            msg = 'The file "' + file_path + '" was not loaded.'
            if not silenced:
                messagebox.showerror('Error', msg)
            else:
                print(dt.now().time(), msg)
        else:
            msg = 'The operation has been cancelled. No data has been saved to '+self.name+' spectrum.'
            if not silenced:
                messagebox.showwarning('Canceled', msg)
            else:
                print(dt.now().time(), msg)
        if overwrite:
            self._orig = self._data
        if '_temp.csv' in file_path:
            os.remove(file_path)
            print(dt.now().time(), 'deleted temporary file', file_path)
        self.min, self.max = self._data.min(0)['wavenumber'], self._data.max(0)['wavenumber']

    def get_data(self):
        gc.collect()

        return self._data.copy()

    def get_orig(self):
        gc.collect()

        return self._orig.copy()

    def set_data_from_spectrum(self, source):
        self._data = source
        self.min, self.max = self._data.min(0)['wavenumber'], self._data.max(0)['wavenumber']

    def normalize_spectrum(self, min_ = 0.0, max_ = 1.0):
        print(dt.now().time(), 'normalizing',self.name,'...')

        self._data.intensity =\
            (self._data.intensity-self._data.intensity.min())/(self._data.intensity.max() - self._data.intensity.min())\
            * (max_ - min_) + min_
        return self.get_data()

    def trim(self, min_, max_):
        print(dt.now().time(), 'trimming', self.name,'to (', min_,',',max_, ')...')
        data = self.get_data()
        self._data=data[(data.wavenumber>min_) & (data.wavenumber<max_)]
        print(self._data.wavenumber.min(),self._data.wavenumber.max())
        return self._data

    def get_length(self):
        return len(self.get_data())

    def identify_peaks(self):
        print(dt.now().time(), 'identifying', self.name, 'peaks...')
        data = self.get_data()
        data['wavenumber'] = data.wavenumber[(data.intensity.shift(1)<data.intensity) & (data.intensity.shift(-1)<data.intensity)]
        data['intensity'] = np.nan
        data.drop_duplicates(subset='wavenumber', keep=False, inplace=True)
        self.peak_list = data
        return data

    def determine_baseline(self):
        # TODO: improve fit under peaks
        print(dt.now().time(), 'determining', self.name, 'baseline...')
        intensity = self._data.intensity
        intensity = intensity.apply(reject_greater,args=[0.08])
        #intensity = intensity.interpolate(method='polynomial',order=3)
        return intensity.rolling(15000)

    def identify_greater_than_baseline(self):
        baseline = self.determine_baseline().mean()
        greaters = self.get_data()
        greaters = greaters[greaters['intensity'] > baseline]
        return greaters

    def detect_peak(self, frequency): # returns min and max frequency above threshold
        # freq: centerpoint of peak
        # tolerance: absolute tolerance of peak frequency
        min_freq, max_freq = frequency-self.tolerance, frequency+self.tolerance
        data = self.get_data()
        peak = data[(data.wavenumber >= min_freq) & (data.wavenumber <= max_freq)]
        if len(peak['wavenumber'])>0:
            if len(peak['wavenumber'])==1:
                peak = peak['wavenumber'].values[0]
            else:
                peak = min(peak['wavenumber'].values, key=lambda x: abs(x - frequency))
            # df['min'] = df.data[(df.data.shift(1) > df.data) & (df.data.shift(-1) > df.data)]
            # peak is the peak frequency
            # find greater min
            try:
                temp_frame = data[data.wavenumber > peak]
                upper_bound = min(temp_frame.wavenumber[(temp_frame.intensity.shift(1) > temp_frame.intensity) & (
                        temp_frame.intensity.shift(-1) > temp_frame.intensity)].values)
                temp_frame = data[data.wavenumber < peak]
                lower_bound = max(temp_frame.wavenumber[(temp_frame.intensity.shift(1) > temp_frame.intensity) & (
                        temp_frame.intensity.shift(-1) > temp_frame.intensity)].values)
                intensity = data[data.wavenumber == peak]['intensity'].values[0]
                return peak, lower_bound, upper_bound, intensity
            except ValueError:
                return np.nan, np.nan, np.nan, np.nan
        else:
            return np.nan, np.nan, np.nan, np.nan

    def detect_peak2(self, frequency): # returns min and max frequency above threshold
        # freq: centerpoint of peak
        # tolerance: absolute tolerance of peak frequency
        min_freq, max_freq = frequency-self.tolerance, frequency+self.tolerance
        data = self.get_data()
        peak = data[(data.wavenumber >= min_freq) & (data.wavenumber <= max_freq)]
        if len(peak['wavenumber'])>0:
            if len(peak['wavenumber'])==1:
                peak = peak['wavenumber'].values[0]
            else:
                peak = min(peak['wavenumber'].values, key=lambda x: abs(x - frequency))
            # df['min'] = df.data[(df.data.shift(1) > df.data) & (df.data.shift(-1) > df.data)]
            # peak is the peak frequency
            # find greater min
            try:
                temp_frame = data[data.wavenumber > peak]
                upper_bound = min(temp_frame.wavenumber[(temp_frame.intensity.shift(1) > temp_frame.intensity) & (
                        temp_frame.intensity.shift(-1) > temp_frame.intensity)].values)
                temp_frame = data[data.wavenumber < peak]
                lower_bound = max(temp_frame.wavenumber[(temp_frame.intensity.shift(1) > temp_frame.intensity) & (
                        temp_frame.intensity.shift(-1) > temp_frame.intensity)].values)
                intensity = data[data.wavenumber == peak]['intensity'].values[0]
                return lower_bound, upper_bound, intensity
            except ValueError:
                return np.nan, np.nan, np.nan
        else:
            return np.nan, np.nan, np.nan

    def get_peak_tuples(self):
        print(dt.now().time(), 'characterizing', self.name, 'peaks...')
        peaks = self.identify_peaks().wavenumber.values.tolist()
        print(dt.now().time(), 'processing', self.name, 'peaks...')
        # peak_data = [self.detect_peak(peak, tolerance) for peak in peaks]
        peak_data = list(map(self.detect_peak, peaks))
        peak_tuples = pd.DataFrame(data=peak_data, columns=['peak', 'lower_bound', 'upper_bound', 'intensity'])
        print(dt.now().time(), 'removing invalid data...')
        peak_tuples = peak_tuples[~peak_tuples['peak'].isna()]
        self.peak_tuples = peak_tuples
        return peak_tuples

    def get_peak_tuples2(self):
        print(dt.now().time(), 'characterizing', self.name, 'peaks...')
        peaks = self.identify_peaks().wavenumber.values.tolist()
        print(dt.now().time(), 'mapping', self.name, 'peaks...')
        # peak_data = [self.detect_peak(peak, tolerance) for peak in peaks]
        peak_data = list(map(self.detect_peak2, peaks)) # i think this is the slow step
        print(dt.now().time(), 'processing', self.name, 'peaks...')
        peak_tuples = pd.DataFrame(data=peak_data, columns=['lower_bound', 'upper_bound', 'intensity'])
        print(dt.now().time(), 'removing invalid data...')
        peak_tuples = peak_tuples[ ~peak_tuples['lower_bound'].isna() & ~peak_tuples['upper_bound'].isna()]
        self.peak_tuples2 = peak_tuples
        return peak_tuples

    def get_downsampled_instance(self, num_points):
        # TODO: make this not so shitty
        print(dt.now().time(), 'getting downsampled instance of', self.name)
        data = self.get_data()
        length = self.get_length()
        if num_points > length:
            num_points = length
        ratio = length // num_points
        if num_points == length:
            return data
        else:
            return data[data.index.values % ratio == 0]

    def calibrate(self, factor):
        self._data.wavenumber *= factor

    def plot(self, num_points=None, show_baseline=False):
        if num_points is None:
            num_points = len(self.get_data())
        smaller = self.get_downsampled_instance(num_points)

        fig = plt.figure()
        axis = fig.add_subplot(111)
        axis.set_xlabel('cm^-1')
        axis.set_ylabel('intensity (a.u.)')
        axis.set_title(self.name)
        try:
            axis.plot(smaller['wavenumber'], smaller['intensity'])
            if show_baseline:
                axis.plot(self._data['wavenumber'], self.determine_baseline().mean())

        except MemoryError:
            print(dt.now().time(), 'unable to plot. not enough memory. try spectrum.plot(n) where n is smaller than before.')
        return fig

    def save_spectrum(self, file_path=None, silenced = True):
        if file_path is None:
            if not silenced:
                root.deiconify()
            file_path = filedialog.asksaveasfilename(title='Save As')

        extension=file_path[file_path.rfind('.'):]
        if extension =='.dpt':
            self.get_data().to_csv(file_path,sep='\t',columns=self.cols, index=False, header=False)


class SpectrumPair(object):
    def __init__(self, primary, secondary, name=None):
        if name is not None:
            self.name = name
        else:
            self.name = randomString()
        if type(primary) is Spectrum:
            self.primary = Spectrum(self.name+', '+primary.name)
            self.primary.set_data_from_spectrum(primary.get_data())
            self.min, self.max = primary.get_data().min(0)['wavenumber'], primary.get_data().max(0)['wavenumber']
        else:
            raise TypeError('primary must be Spectrum object')
        if type(secondary) is Spectrum:
            self.secondary = Spectrum(self.name + ', ' + secondary.name)
            self.secondary.set_data_from_spectrum(secondary.get_data())
        else:
            raise TypeError('secondary must be Spectrum object')

        self.reconciled = False
        self.combined = None

    def remove_secondary(self, scale = 1.0, peak_width = 0.04, intensity_threshold=0.0, **kwargs):

        print(dt.now().time(), 'normalizing...')
        self.primary.normalize_spectrum()
        print(dt.now().time(), 'trimming...')
        self.secondary.trim(self.min, self.max)
        self.secondary.normalize_spectrum()
        print(dt.now().time(), 'finding peaks...')
        data = self.primary.get_data()

        secPeaks = self.secondary.get_data()
        secPeaks = secPeaks[secPeaks.intensity>intensity_threshold]
        secPeaks = secPeaks[['wavenumber','intensity']]
        print('peaks detected: ', len(secPeaks))
        print(dt.now().time(), 'removing peaks...')

        std = peak_width/2
        base = self.primary.determine_baseline()
        baseline = base.mean()
        basestd=base.std()
        noise = pd.Series(np.random.normal(baseline,basestd/4),data.index)

        arrPeaks = list(np.asarray(secPeaks).transpose((1, 0)))

        print(dt.now().time(), 'preparing to generate gaussians...')
        gauss_peaks_locations = [data.copy()[(data.wavenumber >= wavenumber - std) & (
                    data.wavenumber <= wavenumber + std)] for wavenumber, intensity in zip(*arrPeaks)]
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
        self.primary.set_data_from_spectrum(data)
        base = self.primary.determine_baseline()
        baseline = base.mean()
        minbaseline = baseline.dropna().index.min()
        #data.intensity -= baseline
        data.intensity[data.intensity < (noise-0.1)] = noise
        data.intensity[data.intensity > 1] = noise
        data = data[data.index > minbaseline]
        self.primary.set_data_from_spectrum(data)

        #self.primary.normalize_spectrum(0,1)

    def plot(self,num_points = None):
        if self.combined is not None:
            self.combined.plot(num_points)
        self.primary.plot(num_points)
        self.secondary.plot(num_points)

    def subtract_and_plot(self, **kwargs):
        self.remove_secondary(**kwargs)
        self.plot(**kwargs)

    def save_spectrum(self, file_path=None, silenced = True):
        self.primary.save_spectrum(file_path,silenced)


gc.enable()
root = Tk()
root.withdraw()
# temporary stuff
if __name__ == '__main__':
    path = './'
    test = Spectrum('test')
    isotope_mods = {
        1: 1,
        2: 1,
        3: 1,
        4: 1e4,
        5: 1e4,
        6: 1e4,
        7: 1e4
    }
    isotope_mods_hcl = {
        1: 1,
        2: 1,
        3: 1e4,
        4: 1e4,
    }
    test.select_data_from_file(path+'high-DVA_Av2_w_BKG500_cal.dpt',False,True)
    real = Spectrum('simulation')
    real.select_data_from_file(path+'syn-DVA.dpt', False,True)
    water = Spectrum('water')
    water.select_data_from_file(path+'5df03339.par', silenced=True, id_mods=isotope_mods)
    hcl = Spectrum('HCl-DCl')
    hcl.select_data_from_file(path+'hcldcl.par', silenced=True, id_mods=isotope_mods_hcl)
    hcl.trim(90, 600)
    water2 = water
    test.trim(90, 600)
    real.trim(90, 600)
    water.trim(test.min, test.max)
    hcl.normalize_spectrum()
    test.normalize_spectrum()
    water.normalize_spectrum()
    real.normalize_spectrum()
    pair = SpectrumPair(test, water, 'Combined Spectrum')

    pair.plot()
    real.plot()

    mag, width, thresh = 1e4, 0.035, 1e-10
    hclmag, hclwidth = 5e2, 0.02
    pair.remove_secondary(mag, width, thresh)
    pair2 = SpectrumPair(pair.primary, hcl, 'HCl-DCl Included')
    pair2.remove_secondary(hclmag, width, thresh)
    title = './spectrum_mag{}_width{}_thresh{}.dpt'.format(mag, width, thresh)
    title2 = './spectrum_rem_hcl_dcl_mag{}_width{}_water_mag{}_width{}_thresh{}.dpt'.format(hclmag, mag, hclwidth, width, thresh)
    print('saving spectrum as', title)
    pair.save_spectrum(title)
    pair2.save_spectrum(title2)
    pair.plot()
    pair2.plot()
