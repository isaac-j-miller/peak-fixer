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
from sklearn import preprocessing as pp

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
              'intensity':{'type':'float32',
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

    def select_data_from_file(self, file_path=None, headers=None, silenced=False, overwrite=True):
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
            try:

                file = open(file_path,'r')
                data = file.readlines()
                file.close()
                to_pd = [{item: line[par_template[item]['start']:par_template[item]['start']+par_template[item]['width']] for item in par_template} for line in data]

                dtypes = {item: par_template[item]['type'] for item in self.cols}
                self._data = pd.DataFrame.from_dict(to_pd)[self.cols]
                del to_pd, data, file
                success = True

            except:
                success = False
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
                print(dt.now().time(),msg)
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

    def set_data_from_spectrum(self,source):
        self._data = source
        self.min, self.max = self._data.min(0)['wavenumber'], self._data.max(0)['wavenumber']

    def normalize_spectrum(self, min_ = 0.0, max_ = 1.0):
        print(dt.now().time(), 'normalizing',self.name,'...')
        data = self.get_data()
        min_max_scalar = pp.MinMaxScaler(feature_range=(min_,max_))
        np_scaled = min_max_scalar.fit_transform(data)
        del min_max_scalar
        temp = pd.DataFrame(np_scaled,columns=self.cols)
        del np_scaled
        normalized = data
        del data
        normalized['intensity'] = temp['intensity']
        del temp
        self._data = normalized
        del normalized
        return self._data

    def trim(self, min_, max_):
        data = self.get_data()
        print(dt.now().time(), 'trimming', self.name, '...')
        divisor = 1
        done = False
        while not done:
            print(dt.now().time(), 'divisor is',divisor)
            try:
                chunks = chunker(data,divisor)
                inputs = [chunk[(chunk.wavenumber > min_) & (chunk.wavenumber < max_)] for chunk in chunks]
                fname = saver(inputs)
                self.select_data_from_file(fname, True, True, False)
                self.min, self.max = self._data.min(0)['wavenumber'], self._data.max(0)['wavenumber']
                return self._data
            except MemoryError:
                print(dt.now().time(), 'failed. increasing divisor')
                divisor += 1
            # self._data = data

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

    def identify_greater_than_baseline(self, baseline=.01):
        self.normalize_spectrum()
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

    def plot(self, num_points=None):
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
        except MemoryError:
            print(dt.now().time(), 'unable to plot. not enough memory. try spectrum.plot(n) where n is smaller than before.')
        return fig


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

    def interpolate_spectra(self):
        if not self.reconciled:
            self.reconcile_x()
        print(dt.now().time(), 'interpolating', self.primary.name)
        # self.primary.set_data_from_spectrum(self.primary.get_data().interpolate('index', inplace=False))
        self.primary.set_data_from_spectrum(self.primary.get_data().interpolate('polynomial', inplace=False, order=2, limit_area='inside')) # test version
        print(dt.now().time(), 'interpolating', self.secondary.name)
        # self.secondary.set_data_from_spectrum(self.secondary.get_data().interpolate('index', inplace=False))
        self.secondary.set_data_from_spectrum(self.secondary.get_data().interpolate('polynomial', inplace=False, order=2, limit_area='inside')) # test version
        self.primary.normalize_spectrum()
        self.secondary.normalize_spectrum()

    def reconcile_x(self):
        print(dt.now().time(), 'setting up for reconciliation...')
        self.secondary.trim(self.min, self.max)
        pri = self.primary.get_data()
        sec = self.secondary.get_data()
        pri['intensity'] = np.nan
        sec['intensity'] = np.nan
        print(dt.now().time(), 'appending...')
        primary_base = self.primary.get_data().append(sec)
        secondary_base = self.secondary.get_data().append(pri)
        del pri, sec

        print(dt.now().time(), 'sorting primary spectrum...')
        primary_base.sort_values(by='wavenumber', kind='quicksort', inplace=True)
        print(dt.now().time(), 'sorting secondary spectrum...')
        secondary_base.sort_values(by='wavenumber', kind='quicksort', inplace=True)
        # print(dt.now().time(), 'dropping duplicates from primary spectrum...')
        # primary_base.drop_duplicates(inplace=False, subset='wavenumber')
        # print(dt.now().time(), 'dropping duplicates from secondary spectrum...')
        # secondary_base.drop_duplicates(inplace=False, subset='wavenumber')
        print(dt.now().time(), 'setting values...')
        self.primary.set_data_from_spectrum(primary_base)
        self.secondary.set_data_from_spectrum(secondary_base)
        print(dt.now().time(), 'reconciliation complete.')
        self.reconciled = True

    def remove_secondary(self, scale=1.0):
        print(dt.now().time(), 'finding peaks...')
        self.primary.normalize_spectrum()
        self.secondary.trim(self.min,self.max)
        self.secondary.normalize_spectrum(0, scale)
        sec_peaks = self.secondary.get_peak_tuples()[:10] # remove slicer
        data = self.primary.get_data()
        print(dt.now().time(), 'calculating gaussian curves...')

        temp = [[[data.index[index], gauss(x_val, (peak[1]['upper_bound'] - peak[1]['lower_bound']) / 6, peak[1]['intensity'], peak[1]['peak'])]
                for index, x_val in zip(
                data[(data.wavenumber >= peak[1]['lower_bound']) & (data.wavenumber <= peak[1]['upper_bound'])].index,
                data[(data.wavenumber >= peak[1]['lower_bound']) & (data.wavenumber <= peak[1]['upper_bound'])].wavenumber)]
                for peak in sec_peaks.iterrows()]
        print(dt.now().time(), 'transposing...')

        temp = np.vstack(temp).transpose()
        print(dt.now().time(), 'processing gaussian data...')
        temp = pd.DataFrame(temp[1], temp[0], columns=['intensity'])
        # insert line which sums intensity of temp rows with duplicate index
        self.inspect = temp
        print(dt.now().time(), 'removing secondary peaks...')
        data['intensity'] = data['intensity'] - temp['intensity']
        data.fillna(self.primary.get_data()['intensity'])
        print(dt.now().time(), 'done subtracting.')
        self.combined = data
        return data

    def zero(self, data, row):
        data.intensity[(data.wavenumber > row.lower_bound) & (data.wavenumber < row.upper_bound)] = 0
        return None

    def simple_remove_secondary(self, scale = 1.0, peak_width = 0.08, intensity_threshold=0.001):
        peak_width_2 =peak_width/2
        print(dt.now().time(), 'normalizing...')
        self.primary.normalize_spectrum()
        print(dt.now().time(), 'trimming...')
        self.secondary.trim(self.min, self.max)
        self.secondary.normalize_spectrum(0, scale)
        print(dt.now().time(), 'finding peaks...')
        #pri_peaks = self.primary.get_peak_tuples2()
        data = self.primary.get_data()

        #data['d1'] = data.intensity.diff()
        #data['d2'] = data.d1.diff()
        zero_thresh = 1E-6
        neg_zth = zero_thresh*-1
        #minima = data[(data.d1 <= zero_thresh) & (data.d1 >= neg_zth) & (data.d2 > 0)]
        secPeaks = self.secondary.get_data()
        secPeaks = secPeaks[secPeaks.intensity>intensity_threshold]
        peaksRemoved = 0
        #print(minima.head(10), len(minima))
        print(secPeaks.head(10), len(secPeaks))
        print(dt.now().time(), 'removing peaks...')
        """
        for i in range(len(minima)-1):
            if len(secPeaks):
                low = minima.wavenumber.iloc[i]
                li = minima.index[i]
                high = minima.wavenumber.iloc[i + 1]
                hi = minima.index[i + 1]
                #  print(low, high)
                try:
                    peak = secPeaks[(secPeaks.wavenumber >= low) & (secPeaks.wavenumber <= high)].wavenumber.values[0]
                    #TODO: Look for first peak in secPeaks in criteria then remove others up to high
                    secPeaks = secPeaks[secPeaks.wavenumber > high]
                    print(low, high, peak)
                    data.intensity[(data.index >= li) & (data.index <= hi)] *= (1 - scale)
                    peaksRemoved += 1
                except IndexError:
                    break
            else:
                break
        """

        backup = data.intensity.copy()
        std = peak_width/2
        for index,row in secPeaks.iterrows():

            gauss_data = data.copy()[(data.wavenumber >= row['wavenumber']-peak_width) & (data.wavenumber <= row['wavenumber']+peak_width)]
            gauss_data['intensity']=gauss_data.apply(df_gauss, axis=1,args=(std,row['intensity']*scale, row['wavenumber']))
            print('peak:',row['wavenumber'])
            print(gauss_data)
            data.intensity-=gauss_data.intensity
            #data.intensity.subtract(gauss_data.intensity)
            data.intensity=data.intensity.fillna(backup)
            peaksRemoved+=1
        print(dt.now().time(), 'peaks removed:',peaksRemoved)
        self.primary.set_data_from_spectrum(data)

    def subtract_spectrum(self, scale=1):
        print(dt.now().time(), 'beginning subtraction...')
        # just does the bare minimum of subtracting a scaled version of the second spectrum
        if self.reconciled:
            print(dt.now().time(), 'normalizing...')
            self.primary.normalize_spectrum()
            pair.secondary.trim(self.min,self.max)
            self.secondary.normalize_spectrum(0,scale)
            sub = pd.DataFrame(columns=['wavenumber', 'intensity'])
            sub['wavenumber'] = self.primary.get_data()['wavenumber']
            print(dt.now().time(), 'subtracting...')
            # chunk it
            done = False
            divisor = 1
            while not done:
                print(dt.now().time(), 'divisor is', divisor)
                try:
                    sub_chunks = chunker(sub, divisor)
                    sec_chunks = chunker(self.secondary.get_data(), divisor)
                    for i in range(len(sub_chunks)):
                        sub_chunks[i]['intensity'] = sub_chunks[i].subtract(sec_chunks[i])['intensity']
                    fname = saver(sub_chunks)
                    done = True
                except MemoryError:
                    print(dt.now().time(), 'failed. increasing divisor...')
                    divisor += 1
            del sub
            sub = retrieve_from_saver(fname)
            self.combined = Spectrum(self.primary.name+' - '+self.secondary.name)
            self.combined.set_data_from_spectrum(sub)
            return self.combined.get_data()

        else:
            print(dt.now().time(), 'not reconciled. reconciling...')
            self.reconcile_x()
            self.subtract_spectrum() # currently broken

    def plot(self,num_points = None):
        if self.combined is not None:
            self.combined.plot(num_points)
        self.primary.plot(num_points)
        self.secondary.plot(num_points)

    def easy_subtract_and_plot(self, num_points = None):
        # TODO: plot both spectra
        primary_plot = self.primary.plot(num_points)
        self.secondary.plot(num_points)
        # TODO: pick which peak to use to normalize spectrum to be subtracted

        # TODO: plot combined spectrum

    def save_spectrum(self):
        pass
        # TODO: select which to save
        # TODO: select filetype from list
        # TODO: input filename
        # TODO: save spectrum


gc.enable()
root = Tk()
root.withdraw()
# temporary stuff
path = 'F:/Isaac/peak_fixer/'
test = Spectrum('test')
#test.select_data_from_file(path+'GA_0.031torr_72M_10X_600LPF_337.95K_02-09-19_11.12_Emission back parallel input_Si Bolometer [External Pos.5]_0.000960_80 kHz_6 micron Mylar_BEB_AVG201.dpt',False,True)
test.select_data_from_file(path+'CUTDVA_Av2_w_BKG500.csv',False,True)
water = Spectrum('water')
water.select_data_from_file(path+'5d262e45.par', silenced=True)
test.trim(100, 585)
water.trim(test.min, test.max)
test.normalize_spectrum()
water.normalize_spectrum()
pair = SpectrumPair(test, water, 'Combined Spectrum')


# pair.remove_secondary(2.33)
# water.plot()
# print(dt.now().time(), water.detect_peak(202.9, .1))
# pair.interpolate_spectra()
# pair.subtract_spectrum()
pair.plot()
pair.simple_remove_secondary(3.0,0.08, 0.001)
pair.plot()
# print(dt.now().time(), done)
# print(dt.now().time(), 'duplicates:', len(test.get_data())+len(water.get_data())-len(unique))
