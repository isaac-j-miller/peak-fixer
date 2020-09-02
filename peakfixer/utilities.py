from peakfixer import hapi
import requests
import random
import string
import math
import pandas as pd
from datetime import datetime as dt
import os

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

illegal_pathchars = '\'"*/\\[]:;|,<>'


def getAllIsotopes(molec_name:str):
    allIsos = [str(k[1]) for k, v in hapi.ISO.items() if v[4]==molec_name]
    return ','.join(allIsos)


def replaceIllegalPathChars(string):
    d = {k:'_' for k in illegal_pathchars}
    return replaceMany(string, d)

def replaceMany(string, replacementDict):
    s = string
    for k, v in replacementDict.items():
        s= s.replace(k, v)
    return s


def randomString(stringLength=10):
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))


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


def gauss(x, stdev, maximum, mean):
    return maximum*math.exp((-1*(x-mean)**2)/(2*stdev**2))


def df_gauss(row, stdev, maximum, mean):
    return gauss(row['wavenumber'], stdev, maximum, mean)


def reject_greater(x, thresh):
    if x > thresh:
        return 0
    if x < thresh*-1:
        return 0
    return x


class SpectrumGetter:
    def __init__(self, molec_name, min_wavelength=0, max_wavelength=None, isotopologues=None, timeout=None):
        self.molec_name=molec_name
        self.min_wavelength=min_wavelength
        self.max_wavelength=max_wavelength
        self.data=None
        self.molec_ids=[]
        self.timeout=timeout
        self._get_molec_ids(isotopologues)

    def _get_molec_ids(self, isotopologues):
        valids = {k:v[1] for k,v in hapi.ISO_ID.items() if v[5] == self.molec_name}

        if isotopologues is None:
            self.molec_ids=[str(k) for k in valids.keys()]
        else:
            self.molec_ids=[str(k) for k, v in valids.values() if v in isotopologues]

    def getURL(self):
        max_bit = f'&numax={self.max_wavelength}' if self.max_wavelength is not None else ''
        url=f'https://hitran.org/lbl/api?iso_ids_list={",".join(self.molec_ids)}&numin={self.min_wavelength}{max_bit}'
        req=requests.get(url,timeout=30 if self.timeout is None else self.timeout)
        self.data=req.text[:-1].split('\n')

    def fetch(self):
        self.getURL()
        return self.data


class ParWriter:
    def __init__(self, original_par:list, lines_to_keep:list, filename):
        self.original_par = original_par
        self.lines_to_keep = lines_to_keep
        self.filtered_par = None
        self.filename = filename
        self.output_str = ""

    def filter_par(self):
        par_data = self.original_par
        str_lines_to_keep = ['{:12.6f}'.format(line) for line in self.lines_to_keep]
        self.filtered_par = [line for line in par_data if line['wavenumber'] in str_lines_to_keep]

    def write_par_file(self):
        lines = [''.join(list(line.values())) for line in self.filtered_par]
        self.output_str = '\n'.join(lines)
        with open(self.filename, 'w') as f:
            f.write(self.output_str)
        return self.output_str


