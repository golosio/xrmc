import numpy as np
import math

def Gauss(x, x0, s):
    A = 1.0/math.sqrt(2. * math.pi * s * s)
    x1 = x - x0
    return A * math.exp(-x1 * x1/(2.0 * s * s))

Nx = 100
Ny = 100
dx = 0.12
dy = 0.12
sx = 4.0
sy = 2.0

Image = np.empty((Nx,Ny), dtype=np.float64, order='F')

y0 = (-0.5 * Ny + 0.5) * dy
x0 = (-0.5 * Nx + 0.5) * dx

for ix in range(Nx):
    x = x0 + dx * ix
    for iy in range(Ny):
        y = y0 + dy * iy
        Image[ix, iy] = Gauss(x, 0, sx)*Gauss(y, 0, sy);
        
with open('intensityscreen_image_py.dat', 'wb') as outfile:
    outfile.write(Image.tobytes(order='Any'))
