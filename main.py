from image_moments import create_images
from create_moments import *

galaxy = 'NGC4713'; galtemp = 'ngc4713'
path = '/home/nikki/Documents/Data/VERTICO/' + galtemp + '/'

file_pbcorr = path + galtemp + '_7m+tp_co21_pbcorr_round_k.fits'
file_uncorr = path + galtemp + '_7m+tp_co21_flat_round_k.fits'

# Set some parameters to apply to all images below
refresh = True
overwrite = True
make_cutout = False
sun = True
tosave = True

# Have a first look at the cube to figure out some parameters
#cube = moment_maps(galaxy, path, pbcor, cliplevel=0).readfits()
#plt.imshow(np.sum(cube.data, axis=0))
#plt.figure()
#plt.plot(np.sum(cube.data, axis=(1, 2)))

# Call this to creates fits files of the moment maps, or use the "refresh" and "overwrite" options below
#clipped_cube, mom0, mom1, mom2 = moment_maps(galaxy, path, dv=10, pbcor=pbcor, cliplevel=cliplevel, stokes=stokes, sun=True, tosave=True).calc_moms()

# Call these to create images of the moment maps
#create_images(galaxy, file_pbcorr, file_uncorr, savepath=path, refresh=refresh, overwrite=overwrite,
#              make_cutout=make_cutout, sun=sun, tosave=tosave).moment_zero()
#create_images(galaxy, file_pbcorr, file_uncorr, savepath=path, refresh=refresh, overwrite=overwrite,
#              make_cutout=make_cutout, sun=sun, tosave=tosave).moment_1_2()
#create_images(galaxy, file_pbcorr, file_uncorr, savepath=path, refresh=refresh, overwrite=overwrite,
#              make_cutout=make_cutout, sun=sun, tosave=tosave).moment_1_2(moment=2)
create_images(galaxy, file_pbcorr, file_uncorr, savepath=path, refresh=refresh, overwrite=overwrite,
              make_cutout=make_cutout, sun=sun, tosave=tosave).PVD(findcentre=False, find_velcentre=False, full_width=False)
#create_images(galaxy, file_pbcorr, file_uncorr, savepath=path, refresh=refresh, overwrite=overwrite,
#              make_cutout=make_cutout, sun=sun, tosave=tosave).spectrum()
