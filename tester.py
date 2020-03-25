from peakfixer import *

if __name__ == '__main__':
    d_abundance_rel = 2e4  # Empirically determined abundance of D relative to natural abundance. This number can be changed
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
    # relative abundances of species in the water (numbers correspond to HITRAN file)
    isotope_mods_hcl = {
        1: 1,
        2: 1,
        3: d_abundance_rel,
        4: d_abundance_rel,
    }
    test.select_data_from_file(path+'high-DVA_Av2_w_BKG500_cal.dpt',False,True)  # test sets its spectrum to the data in this path
    simulated = Spectrum('simulation')  # make another spectrum called "simulation". This is not necessary, i just did this to compare the experimental spectrum to the simulated one
    simulated.select_data_from_file(path + 'syn-DVA.dpt', False, True)  # set the data for the simulated spectrum
    water = Spectrum('water')  # make a new spectrum for the water
    water.select_data_from_file(path+'5df03339.par', silenced=True, id_mods=isotope_mods)  # select the HITRAN file for water
    hcl = Spectrum('HCl-DCl')  # make a new spectrum for the hcl/dcl
    hcl.select_data_from_file(path+'hcldcl.par', silenced=True, id_mods=isotope_mods_hcl)  # select the HITRAN file for hcl-dcl
    hcl.trim(90, 600)  # trim it to a wavelength range of interest (wavenumbers)
    test.trim(90, 600)
    simulated.trim(90, 600)
    water.trim(test.min, test.max)  # another way of trimming, test.min and test.max are the min and max wavenumbers of "test"
    # normalize the spectra to be between 0 and 1
    hcl.normalize_spectrum()
    test.normalize_spectrum()
    water.normalize_spectrum()
    simulated.normalize_spectrum()
    pair = SpectrumPair(test, water, 'Combined Spectrum')  # make a spectrum pair to compare the exp. spectrum and water

    # plot the spectra
    pair.plot()
    simulated.plot()

    # set the parameters for peak removal. mag refers to water, hclmag refers to hcl. thresh and width refer to both
    # mag (hclmag) is the relative abundance of the water (hcl) to the species of interest, width (hclwidth) is the
    # assumed peak width of water (hcl), and thresh is the minimum intensity for a peak to be removed.
    mag, width, thresh = 5e6, 0.01, 0
    hclmag= 5e2
    pair.remove_secondary(mag, width, thresh)  # remove the secondary spectrum (water) from the primary (experimental) in a SpectrumPair (see above, where pair was declared)
    pair2 = SpectrumPair(pair.primary, hcl, 'HCl-DCl Included') # make a new pair from the experimental spectrum with water removed (pair.primary) and hcl
    pair2.remove_secondary(hclmag, width, thresh) # do the same with hcl as you did with water.
    # You can use whatever parameters for relative magnitude, assumed peak width, and the peak threshold you want for whatever species you want.
    # They don't have to all be the same

    #title and title2 are the file names of the experimental spectrum with just water removed and with water and hcl/dcl removed, respectively.
    title = './spectrum_mag{}_width{}_thresh{}.dpt'.format(mag, width, thresh)
    title2 = './spectrum_rem_hcl_dcl_mag{}_width{}_water_mag{}_width{}_thresh{}.dpt'.format(hclmag, mag, width, width, thresh)
    print('saving spectrum as', title)
    pair.save_spectrum(title)  # save the spectrum with just water removed as a *.dpt
    print('saving spectrum as', title2)
    pair2.save_spectrum(title2) # save the spectrum with water and hcl/dcl removed as a *.dpt
    # plot the spectrum pairs
    pair.plot()
    pair2.plot()
