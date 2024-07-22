# allam/__init__.py

# Define a package-level variable
__version__ = '1.0.0'
__date__ = '20.07.2024'
__all__ = ['allam', 'combustion', 'recuperator', 'phasediagrCO2', 'property_sCO2_cp']

# from .allam import *
# from .combustion import *
# from .recuperator import *
from .phasediagrCO2 import *
from .spHvol import *
from .property_sCO2_cp import *

