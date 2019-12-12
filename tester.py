from peakfixer import *

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
