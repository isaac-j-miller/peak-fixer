from .utilities import *
import numpy as np
import matplotlib.pyplot as plt


class Spectrum(object):
    spectra = []

    def __init__(self, name = None):
        if name is not None and name not in Spectrum.spectra:
            self.name = name
        else:
            self.name = randomString(10)
        Spectrum.spectra.append(self)
        self._data = None
        self.cols = ['wavenumber', 'intensity']
        self._orig = None
        self.min=0.0
        self.max=0.0
        self.par_data = None
        self.tempfilename= f'./hitran/{replaceIllegalPathChars(self.name)}.par'
    
    def __del__(self):
        if os.path.exists(self.tempfilename):
            os.remove(self.tempfilename)
        if os.path.exists(self.tempfilename.replace('.par','.h5')):
            os.remove(self.tempfilename.replace('.par','.h5'))

    def load_from_par_buffer(self, data, id_mods=None):
        to_pd = [
            {item: line[par_template[item]['start']:par_template[item]['start'] + par_template[item]['width']] for item
             in par_template} for line in data]
        self.cols.append('iso_num')
        # dtypes = {item: par_template[item]['type'] for item in self.cols}
        self._data = pd.DataFrame.from_dict(to_pd)[self.cols]
        if id_mods is not None:
            isonum = self._data.iso_num.astype(int)
            self._data['mods'] = isonum.map(id_mods).astype(float)
            self._data['intensity'] = self._data.mods * self._data.intensity.astype(float)
            self._data['wavenumber'] = self._data.wavenumber.astype(float)
        self.cols.remove('iso_num')
        self.par_data = to_pd
        with open(self.tempfilename, 'w') as f:
            f.write('\n'.join(data))
        

    def select_data_from_hitran(self, molec_name, min_wavelength=0, max_wavelength=None, isotopologues=None, id_mods=None, timeout=None):
        print(dt.now().time(), 'getting data...')
        getter = SpectrumGetter(molec_name, min_wavelength, max_wavelength, isotopologues)
        data=getter.fetch()
        self._orig = data
        self.load_from_par_buffer(data, id_mods)

    def select_data_from_file(self, file_path, headers=True, overwrite=True, id_mods=None):
        print(dt.now().time(), 'loading...')
        index = file_path.rfind('.')
        filetype = file_path[index:].lower()
        if filetype in txt_types.keys():
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
            with open(file_path, 'r') as f:
                data = f.readlines()
            self.load_from_par_buffer(data, id_mods)
            dtypes = {item: par_template[item]['type'] for item in self.cols}
            success = True
        else:
            success = False

        if success:
            self._data = self._data.astype(dtypes)
            msg = 'The file "' + file_path + '" has been loaded into '+self.name+' spectrum.'
            print(dt.now().time(),msg)
        elif success is not None:
            msg = 'The file "' + file_path + '" was not loaded.'
            print(dt.now().time(), msg)
        else:
            msg = 'The operation has been cancelled. No data has been saved to '+self.name+' spectrum.'
            print(dt.now().time(), msg)
        if overwrite:
            self._orig = self._data
        if '_temp.csv' in file_path:
            os.remove(file_path)
            print(dt.now().time(), 'deleted temporary file', file_path)
        self.min, self.max = self._data.min(0)['wavenumber'], self._data.max(0)['wavenumber']

    def get_data(self):
        return self._data.copy()

    def get_orig(self):
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

    def determine_baseline(self):
        # TODO: improve fit under peaks
        print(dt.now().time(), 'determining', self.name, 'baseline...')
        intensity = self._data.intensity
        intensity = intensity.apply(reject_greater,args=[0.08])
        return intensity.rolling(15000)

    def identify_greater_than_baseline(self):
        baseline = self.determine_baseline().mean()
        greaters = self.get_data()
        greaters = greaters[greaters['intensity'] > baseline]
        return greaters

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

    def save_spectrum(self, file_path):
        extension=file_path[file_path.rfind('.'):]
        if extension == '.dpt':
            self.get_data().to_csv(file_path,sep='\t',columns=self.cols, index=False, header=False)
