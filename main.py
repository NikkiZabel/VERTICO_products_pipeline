from image_moments import CreateImages
from create_moments import *
import warnings; warnings.filterwarnings("ignore")
from matplotlib import pyplot as plt

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
#cube_corr, cube_uncorr = ClipCube(galaxy, file_pbcorr, file_uncorr).readfits()
#plt.imshow(np.sum(cube_corr.data, axis=0))
#plt.figure()
#plt.plot(np.sum(cube_corr.data, axis=(1, 2)))

# Call this to creates fits files of the moment maps, or use the "refresh" and "overwrite" options below
#clipped_cube, mom0, mom1, mom2 = moment_maps(galaxy, path, dv=10, pbcor=pbcor, cliplevel=cliplevel, stokes=stokes, sun=True, tosave=True).calc_moms(units='M_Sun/pc^2', alpha_co=6.25)

# Call th#ese to create images of the moment maps
#CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=path, refresh=refresh, overwrite=overwrite,
#              make_cutout=make_cutout, sun=sun, tosave=tosave).moment_zero(units='K km/s')
#CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=path, refresh=refresh, overwrite=overwrite,
#             make_cutout=make_cutout, sun=sun, tosave=tosave).moment_1_2()
#CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=path, refresh=refresh, overwrite=overwrite,
#              make_cutout=make_cutout, sun=sun, tosave=tosave).moment_1_2(moment=2)
#CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=path, refresh=refresh, overwrite=overwrite,
#              make_cutout=make_cutout, sun=sun, tosave=tosave).\
#    PVD(axis='minor', full_width=False)
#CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=path, refresh=refresh, overwrite=overwrite,
#              make_cutout=make_cutout, sun=sun, tosave=tosave).spectrum(x_axis='velocity')
#CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=path, refresh=refresh, overwrite=overwrite,
#            make_cutout=make_cutout, sun=sun, tosave=tosave).radial_profile(units='kpc',
#                                    alpha_co=6.25, table_path='/home/nikki/Documents/Data/VERTICO/VERTICO_master.fits',
#                                                                        check_aperture=True)
