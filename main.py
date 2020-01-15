from create_moments import *
from gal_params import parameters
from image_moments import create_images
from matplotlib import pyplot as plt
import numpy as np

# Basic info

galaxy = 'NGC4713'

path = '/home/nikki/Documents/Data/VERTICO/'
savepath = path + '/' + galaxy + '/'

# Have a first look at the cube to figure out some parameters

centre_x, centre_y, size, vrange, vrange_2, pbcor, cliplevel, stokes, start, stop, dv = parameters(galaxy)

#cube = moment_maps(galaxy, path, pbcor, cliplevel=0, dv=10).readfits()
#plt.imshow(np.sum(cube.data, axis=0))
#plt.figure()
#plt.plot(np.sum(cube.data, axis=(1, 2)))

#mom0, mom1, mom2 = moment_maps(galaxy, path, dv=10, pbcor=pbcor, cliplevel=cliplevel, stokes=stokes, tosave=True).calc_moms()


create_images(savepath, centre_y, centre_x, size, savepath=savepath, tosave=True).moment_zero()
create_images(savepath, centre_y, centre_x, size, savepath=savepath, tosave=True).moment_1_2(vrange=vrange)
create_images(savepath, centre_y, centre_x, size, savepath=savepath, tosave=True).moment_1_2(moment=2, vrange_2=vrange_2)