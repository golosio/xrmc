import numpy as np

def read_intensityscreen_image(filename):
    with open(filename, 'rb') as f:
        return np.fromfile(f, dtype=np.double)
