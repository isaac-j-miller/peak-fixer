from .spectrum import Spectrum
from .spectrumPair import SpectrumPair
from .hapiSpectrum import HapiSpectrum
import os
if not os.path.isdir('./hitran'):
    os.mkdir('./hitran')