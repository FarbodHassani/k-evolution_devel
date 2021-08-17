
import numpy as np
from pygadgetreader import *
import argparse
import sys
sys.path.append('./../../Cosmology_library/')
from cosmology_module import read_gadget


read_file = read_gadget.ReadGadget("every_N");
kind, file_location, file_base, num_ini, num_fin, output_name, param = read_file.read_params()
read_file.choose_every_N(kind, file_location, file_base, num_ini, num_fin, output_name, param);
