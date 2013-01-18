# astero __init__.py
"""
Object-oriented package to deal with asteroseismology frequencies.
It loads frequency sets from some standard file layouts and provides methods
to do asteroseismic analysis.
Currently these are implemented:
    Calculate large separations
    Calculate small separations
    Calculate second differences
    Approximate mean large separation by fitting radial modes
    Echelle plot
    -> other calculations are still buggy
Package depends on 'matplotlib', 'numpy', 'scipy' and 'uncertainties'
"""

from freq_class import _astero
load = _astero
echelle = _astero.echelle
plotls = _astero.plotls
plotss = _astero.plotss
plotr = _astero.plotr
mesa_create_inlist = _astero.mesa_create_inlist
