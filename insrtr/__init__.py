"""Top-level package for INSRTR."""

__author__ = """Ajasja Ljubetic"""
__email__ = 'ajasja.ljubetic@gmail.com'
__version__ = '0.1.0'


from .ccs import COILED_COILS, LINKERS
from .main import *
from .analysis import * # is this best practice? A: You should import only the things you need
from .model import predict_positions
from .models import *
MODEL_NAME = 'gbt_classifier_v1.pkl'