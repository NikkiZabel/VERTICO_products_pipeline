from image_moments import create_images
from create_moments import *

# Basic info
galaxy = 'NGC4713'
path = '/home/nikki/Documents/Data/VERTICO/' + galaxy + '/'

# Have a first look at the cube to figure out some parameters
#cube = moment_maps(galaxy, path, pbcor, cliplevel=0, dv=10).readfits()
#plt.imshow(np.sum(cube.data, axis=0))
#plt.figure()
#plt.plot(np.sum(cube.data, axis=(1, 2)))

# Call this to creates fits files of the moment maps, or use the "refresh" and "overwrite" options below
#clipped_cube, mom0, mom1, mom2 = moment_maps(galaxy, path, dv=10, pbcor=pbcor, cliplevel=cliplevel, stokes=stokes, tosave=True).calc_moms()

# Call these to create images of the moment maps
create_images(galaxy, path, savepath=path, refresh=True, overwrite=True, make_cutout=False, tosave=True).moment_zero()
create_images(galaxy, path, savepath=path, refresh=True, overwrite=True, make_cutout=False, tosave=True).moment_1_2()
create_images(galaxy, path, savepath=path, refresh=True, overwrite=True, make_cutout=False, tosave=True).moment_1_2(moment=2)
create_images(galaxy, path, savepath=path, refresh=True, overwrite=True, make_cutout=False, tosave=True).PVD(findcentre=False, find_velcentre=False, full_width=False)
create_images(galaxy, path, savepath=path, refresh=True, overwrite=True, make_cutout=False, tosave=True).spectrum()