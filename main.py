from image_moments import CreateImages
from create_moments import *
import warnings; warnings.filterwarnings("ignore")
import os
from matplotlib import pyplot as plt

# Set some parameters to apply to all images below
redo_clip = False
refresh = True
overwrite = True
sun = True
tosave = True
pbcor = True
resolution = 0
version = '1_1'
sample = 'vertico'

galaxies = ['IC3392', 'NGC4064', 'NGC4189', 'NGC4192', 'NGC4216', 'NGC4222', 'NGC4294', 'NGC4299', 'NGC4302',
            'NGC4330', 'NGC4351', 'NGC4380', 'NGC4383', 'NGC4388', 'NGC4394', 'NGC4405', 'NGC4419', 'NGC4522',
            'NGC4532', 'NGC4533', 'NGC4568', 'NGC4606', 'NGC4607', 'NGC4651', 'NGC4713', 'NGC4808', 'NGC4396',
            'NGC4567', 'NGC4772', 'NGC4580', 'NGC4450', 'NGC4254', 'NGC4293', 'NGC4298', 'NGC4321', 'NGC4402',
            'NGC4424', 'NGC4457', 'NGC4535', 'NGC4536', 'NGC4548', 'NGC4569', 'NGC4579', 'NGC4654', 'NGC4689',
            'NGC4691', 'NGC4698', 'NGC4694']

#galaxies = ['IC3392', 'NGC4064', 'NGC4189', 'NGC4192', 'NGC4216', 'NGC4222', 'NGC4294', 'NGC4299', 'NGC4302',
#            'NGC4330', 'NGC4351', 'NGC4380', 'NGC4383', 'NGC4388', 'NGC4394', 'NGC4405', 'NGC4419', 'NGC4522',
#            'NGC4532', 'NGC4533', 'NGC4568', 'NGC4606', 'NGC4607', 'NGC4651', 'NGC4713', 'NGC4808', 'NGC4396',
#            'NGC4567', 'NGC4772', 'NGC4580', 'NGC4450', 'NGC4694', 'NGC4561']

#galaxies = ['NGC4254', 'NGC4293', 'NGC4298', 'NGC4321', 'NGC4402',
#            'NGC4424', 'NGC4457', 'NGC4535', 'NGC4536', 'NGC4548', 'NGC4569', 'NGC4579', 'NGC4654', 'NGC4689']

#galaxies = ['NGC4064', 'NGC4222', 'NGC4294', 'NGC4330', 'NGC4388', 'NGC4394', 'NGC4402', 'NGC4405', 'NGC4419',
#            'NGC4522', 'NGC4533', 'NGC4567', 'NGC4606', 'NGC4607', 'NGC4772']  # These are the 7m only detections

galaxies = ['NGC4380']

for i in range(len(galaxies)):

    galaxy = galaxies[i]

    print(galaxy)

    path = '/home/nikki/Documents/Data/VERTICO/'

    if resolution == 15:
        readpath = path + 'ReducedData/v1_1/15arcsec/' + galaxy + '/'
    elif resolution == 9:
        readpath = path + 'ReducedData/v1_1/9arcsec/' + galaxy + '/'
    else:
        readpath = path + 'ReducedData/v1_1/native/' + galaxy + '/'

    if sun:
        if resolution == 15:
            if not os.path.exists(path + 'products_v' + version + '/15arcsec/' + galaxy + '/'):
                os.mkdir(path + 'products_v' + version + '/15arcsec/' + galaxy + '/')
            savepath_temp = path + 'products_v' + version + '/15arcsec/' + galaxy + '/'
        elif resolution == 9:
            if not os.path.exists(path + 'products_v' + version + '/9arcsec/' + galaxy + '/'):
                os.mkdir(path + 'products_v' + version + '/9arcsec/' + galaxy + '/')
            savepath_temp = path + 'products_v' + version + '/9arcsec/' + galaxy + '/'
        else:
            if not os.path.exists(path + 'products_v' + version + '/native/' + galaxy + '/'):
                os.mkdir(path + 'products_v' + version + '/native/' + galaxy + '/')
            savepath_temp = path + 'products_v' + version + '/native/' + galaxy + '/'
    else:
        if resolution == 15:
            if not os.path.exists(path + 'products_v' + version + '/15_arcsec/dame11_method/' + galaxy + '/'):
                os.mkdir(path + 'products_v' + version + '/15_arcsec/dame11_method/' + galaxy + '/')
            savepath_temp = path + 'products_v' + version + '/15_arcsec/dame11_method/' + galaxy + '/'
        elif resolution == 9:
            if not os.path.exists(path + 'products_v' + version + '/9_arcsec/dame11_method/' + galaxy + '/'):
                os.mkdir(path + 'products_v' + version + '/9_arcsec/dame11_method/' + galaxy + '/')
            savepath_temp = path + 'products_v' + version + '/9_arcsec/dame11_method/' + galaxy + '/'
        else:
            if not os.path.exists(path + 'products_v' + version + '/native/dame11_method/' + galaxy + '/'):
                os.mkdir(path + 'products_v' + version + '/native/dame11_method/' + galaxy + '/')
            savepath_temp = path + 'products_v' + version + '/native/dame11_method/' + galaxy + '/'

    if galaxy == 'NGC4606' or galaxy == 'NGC4351':
        import matplotlib
        matplotlib.rcParams['text.usetex'] = False

    TP = True
    if resolution == 15:
        if sun:
            savepath = savepath_temp + galaxy + '_7m+tp_co21_pbcorr_round_k_15_arcsec_'
        else:
            savepath = savepath_temp + galaxy + '_7m+tp_co21_pbcorr_round_k_15_arcsec_dame11_'
    elif resolution == 9:
        if sun:
            savepath = savepath_temp + galaxy + '_7m+tp_co21_pbcorr_round_k_9_arcsec_'
        else:
            savepath = savepath_temp + galaxy + '_7m+tp_co21_pbcorr_round_k_9_arcsec_dame11_'
    else:
        if sun:
            savepath = savepath_temp + galaxy + '_7m+tp_co21_pbcorr_round_k_'
        else:
            savepath = savepath_temp + galaxy + '_7m+tp_co21_pbcorr_round_k_dame11_'
    try:
        file_pbcorr = readpath + galaxy + '_7m+tp_co21_pbcorr_round_k.fits'
        file_uncorr = readpath + galaxy + '_7m+tp_co21_flat_round_k.fits'
        cube_corr, cube_uncorr = ClipCube(galaxy, file_pbcorr, file_uncorr).readfits()
    except:
        if resolution == 15:
            if sun:
                savepath = savepath_temp + galaxy + '_7m_co21_pbcorr_round_k_15_arcsec_'
            else:
                savepath = savepath_temp + galaxy + '_7m_co21_pbcorr_round_k_15_arcsec_dame11_'
        elif resolution == 9:
            if sun:
                savepath = savepath_temp + galaxy + '_7m_co21_pbcorr_round_k_9_arcsec_'
            else:
                savepath = savepath_temp + galaxy + '_7m_co21_pbcorr_round_k_9_arcsec_dame11_'
        else:
            if sun:
                savepath = savepath_temp + galaxy + '_7m_co21_pbcorr_round_k_'
            else:
                savepath = savepath_temp + galaxy + '_7m_co21_pbcorr_round_k_dame11_'
        try:
            file_pbcorr = readpath + galaxy + '_7m_co21_pbcorr_round_k.fits'
            file_uncorr = readpath + galaxy + '_7m_co21_flat_round_k.fits'
            cube_corr, cube_uncorr = ClipCube(galaxy, file_pbcorr, file_uncorr).readfits()
        except:
            continue
        TP = False

    #savepath = '/home/nikki/Documents/Data/VERTICO/mom0_unclipped/' + galaxy + '_'

    #galaxy = 'NGC4222'
    #file_pbcorr = '/home/nikki/Documents/Data/VERTICO/ReducedData/test/NGC4222_tapered/NGC4222_7m_co21_pbcorr_round_k.fits'
    #file_uncorr = '/home/nikki/Documents/Data/VERTICO/ReducedData/test/NGC4222_tapered/NGC4222_7m_co21_flat_round_k.fits'
    #savepath = '/home/nikki/Documents/Data/VERTICO/ReducedData/test/NGC4222_tapered/Products/Dame_method/'

    '''
    try:
        data = fits.open(file_pbcorr)[0]
        pb = fits.open('/home/nikki/Documents/Data/VERTICO/ReducedData/' + galaxy + '/' +
                           galaxy + '_7m_co21_pb_rebin.fits')[0]

        if not pb.shape[1:3] == data.shape[1:3]:
            print('Wrong PB shape:')
            print(galaxy)
            continue
    except:
        print('No PB info: ')
        print(galaxy)
    '''

    '''
    if not pbcor:
        file_pbcorr = file_uncorr
        cube_corr = cube_uncorr.copy()
        if not os.path.exists(savepath_temp + 'PB_uncorrected/'):
            os.mkdir(savepath_temp + 'PB_uncorrected/')
        if TP:
            savepath = savepath_temp + 'PB_uncorrected/' + galaxy + '_7m+tp_co21_flat_round_k.fits'
        else:
            savepath = savepath_temp + 'PB_uncorrected/' + galaxy + '_7m_co21_flat_round_k.fits'
    '''

    # Have a first look at the cube to figure out some parameters
    #plt.figure()
    #plt.imshow(np.sum(cube_corr.data, axis=0), origin='lower')
    #plt.figure()
    #spec = np.sum(cube_corr.data, axis=(1, 2))
    #plt.plot(spec)
    #std = np.std(spec)
    #x = np.arange(0, len(spec), 1)
    #plt.plot(x, std * np.ones(len(x)))


    # Moment maps
    '''
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                  sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).moment_zero(units='K km/s')
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                  sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).moment_zero()
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                  sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).moment_zero(peak=True)
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).moment_1_2()
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).moment_1_2(moment=2)

    # Uncertainty maps

    try:
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                    sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).mom0_noise_maps()
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                    sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).mom1_2_noise_maps()
    except:
        pass

    #if galaxy == 'NGC4561': continue
    '''
    # PVDs
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).\
       PVD(axis='major', find_angle=False, check_slit=True)
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                  sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).\
        PVD(axis='minor', check_slit=False)
    '''
    # Spectra
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                  sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).spectrum(x_axis='vel_offset')
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).spectrum(x_axis='velocity')
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).spectrum(x_axis='frequency')
    '''
    # Radial profiles
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).radial_profile(x_units='arcsec', y_units='M_Sun pc^-2',
                                    alpha_co=5.4, table_path='/home/nikki/Documents/Data/VERTICO/VERTICO_master.fits',
                                                                            check_aperture=False)
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).radial_profile(x_units='kpc', y_units='M_Sun pc^-2',
                                   alpha_co=5.4, table_path='/home/nikki/Documents/Data/VERTICO/VERTICO_master.fits',
                                                                            check_aperture=False)
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).radial_profile(y_units='K km/s', x_units='kpc',
                                   alpha_co=5.4, table_path='/home/nikki/Documents/Data/VERTICO/VERTICO_master.fits',
                                                                            check_aperture=False)
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).radial_profile(y_units='K km/s', x_units='arcsec',
                                   alpha_co=5.4, table_path='/home/nikki/Documents/Data/VERTICO/VERTICO_master.fits',
                                                                            check_aperture=False)
    #'''
