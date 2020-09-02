#!/usr/bin/env python
# coding: utf-8

# In[1]:


from peakfixer import *


# Import the peakfixer module

# In[2]:


d_abundance_rel = 1  # Empirically determined abundance of D relative to natural abundance. This number can be changed
path = './'  # just means that all the files are in the same folder as this file
test = Spectrum('test')  # instantiate a new spectrum named "test"
# relative abundances of species in the water (numbers correspond to HITRAN file)
isotope_mods = {
    1: 1,
    2: 1,
    3: 1,
    4: d_abundance_rel,
    5: d_abundance_rel,
    6: d_abundance_rel,
    7: d_abundance_rel
}
filename ='unfixed/GA_0.6torr_72M_10X_600LPF_337.95K_02-06-19_13.58_Emission back parallel input_Si Bolometer [External Pos.5]_0.000960_80 kHz_6 micron Mylar_BEB_AVG186.csv'


# set the relative abundances of the different water isotopologues. the numbers correspond to HITRAN

# In[3]:


test.select_data_from_file(path+filename,headers=False)  # test sets its spectrum to the data in this path


# load the file of interest (NOTE: if the file is comma-separated, you must rename it so that it has a .csv extension)

# In[4]:


water = Spectrum('water')  # make a new spectrum for the water
water.select_data_from_file('./peakfixer/water.par', id_mods=isotope_mods)  # select the HITRAN file for water


# load the hitran linelists. the "silenced" argument just means that it doesn't try to do a file dialogue

# In[5]:


test.plot()


# plot the test spectrum to see what it looks like. this one has a really weird baseline where the frequency < 100 cm^-1 and where the frequency > 600 cm-1

# In[6]:


test.trim(188, 214)
water.trim(188, 214)
# normalize the spectra to be between 0 and 1
test.normalize_spectrum()
water.normalize_spectrum()
pair = SpectrumPair(test, water, 'Combined Spectrum')  # make a spectrum pair to compare the exp. spectrum and water


# normalize the spectra and trim the experimental spectrum to exclude the weird part because the peakfixer program can't quite handle it

# In[7]:


pair.plot()


# In[8]:


mag, width, thresh_lower= 1e9, 0.003, 0


# set the parameters for peak removal. mag refers to water, hclmag refers to hcl. thresh and width refer to both.
# mag (hclmag) is the relative abundance of the water (hcl) to the species of interest, width (hclwidth) is the assumed peak width of water (hcl), and thresh is the minimum intensity for a peak to be removed.

# In[9]:


pair.remove_secondary(mag, width, thresh_lower)


# remove the secondary spectrum (water) from the primary (experimental) in a SpectrumPair (see above, where pair was declared)
# You can use whatever parameters for relative magnitude, assumed peak width, and the peak threshold you want for whatever species you want.
# They don't have to all be the same

# In[10]:


pair.plot()


# plot the spectra again after the secondary spectrum has been removed

# In[11]:


title = path +'fixed/0.12torr_spectrum_mag{}_width{}_thresh_lower{}_.dpt'.format(mag, width, thresh_lower)


# come up with a filename for the new spectrum

# In[12]:


print('saving spectrum as', title)
pair.save_spectrum(title)  # save the spectrum with just water removed as a *.dpt


# save the new spectrum
