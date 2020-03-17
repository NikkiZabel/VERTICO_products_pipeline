from image_moments import CreateImages
from create_moments import *
import warnings; warnings.filterwarnings("ignore")
import os
from matplotlib import pyplot as plt

# Set some parameters to apply to all images below
refresh = True
overwrite = True
sun = True
tosave = True

galaxies = ['IC3392', 'NGC4064', 'NGC4189', 'NGC4192', 'NGC4216', 'NGC4222', 'NGC4294', 'NGC4299', 'NGC4302',
            'NGC4330', 'NGC4351', 'NGC4380', 'NGC4383', 'NGC4388', 'NGC4394', 'NGC4405', 'NGC4419', 'NGC4522',
            'NGC4532', 'NGC4533', 'NGC4568', 'NGC4606', 'NGC4607', 'NGC4651', 'NGC4713', 'NGC4808', 'NGC4396',
            'NGC4567', 'NGC4772']

#galaxies = ['NGC4380']

for i in range(len(galaxies)):

    galaxy = galaxies[i]

    path = '/home/nikki/Documents/Data/VERTICO/' + galaxy + '/'

    if not os.path.exists(path + galaxy + '_Products'):
        os.mkdir(path + galaxy + '_Products/')

    if sun:
        if not os.path.exists(path + galaxy + '_Products/Sun_method/'):
            os.mkdir(path + galaxy + '_Products/Sun_method/')
        savepath = path + galaxy + '_Products/Sun_method/'
    else:
        if not os.path.exists(path + galaxy + '_Products/Smooth_method/'):
            os.mkdir(path + galaxy + '_Products/Smooth_method/')
        savepath = path + galaxy + '_Products/Smooth_method/'

    if galaxy == 'NGC4606':
        import matplotlib
        matplotlib.rcParams['text.usetex'] = False

    try:
        file_pbcorr = path + galaxy + '_7m+tp_co21_pbcorr_round_k.fits'
        file_uncorr = path + galaxy + '_7m+tp_co21_flat_round_k.fits'
        cube_corr, cube_uncorr = ClipCube(galaxy, file_pbcorr, file_uncorr).readfits()
        savepath = savepath + galaxy + '_7m+tp_co21_pbcorr_round_k_'
    except:
        file_pbcorr = path + galaxy + '_7m_co21_pbcorr_round_k.fits'
        file_uncorr = path + galaxy + '_7m_co21_flat_round_k.fits'
        cube_corr, cube_uncorr = ClipCube(galaxy, file_pbcorr, file_uncorr).readfits()
        savepath = savepath + galaxy + '_7m_co21_pbcorr_round_k_'

    # Have a first look at the cube to figure out some parameters
    #plt.imshow(np.sum(cube_corr.data, axis=0))
    #plt.figure()
    #plt.plot(np.sum(cube_corr.data, axis=(1, 2)))

    # Call this to creates fits files of the moment maps, or use the "refresh" and "overwrite" options below
    #clipped_cube, mom0, mom1, mom2 = moment_maps(galaxy, path, dv=10, pbcor=pbcor, cliplevel=cliplevel, stokes=stokes, sun=True, tosave=True).calc_moms(units='M_Sun/pc^2', alpha_co=6.25)

    # Call th#ese to create images of the moment maps
    #CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
    #              sun=sun, tosave=tosave).moment_zero(units='K km/s')
    #CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
    #              sun=sun, tosave=tosave).moment_zero()
    #CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
    #             sun=sun, tosave=tosave).moment_1_2()
    #CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
    #              sun=sun, tosave=tosave).moment_1_2(moment=2)
    #CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
    #             sun=sun, tosave=tosave).\
    #   PVD(axis='major', find_angle=False, check_slit=True)
    #CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
    #              sun=sun, tosave=tosave).\
    #    PVD(axis='minor', check_slit=False)
    #CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
    #              sun=sun, tosave=tosave).spectrum(x_axis='vel_offset')
    #CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
    #              sun=sun, tosave=tosave).spectrum(x_axis='velocity')
    #CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
    #             sun=sun, tosave=tosave).spectrum(x_axis='frequency')
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                sun=sun, tosave=tosave).radial_profile(units='arcsec',
                                    alpha_co=6.25, table_path='/home/nikki/Documents/Data/VERTICO/VERTICO_master.fits',
                                                                            check_aperture=False)
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                sun=sun, tosave=tosave).radial_profile(units='kpc',
                                    alpha_co=6.25, table_path='/home/nikki/Documents/Data/VERTICO/VERTICO_master.fits',
                                                                            check_aperture=False)

