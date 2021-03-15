from image_moments import CreateImages
from create_moments import *
import warnings; warnings.filterwarnings("ignore")
import os
from glob import glob

sample = 'heracles'
version = '1_1'
path = '/home/nikki/Documents/Data/VERTICO/heracles/'

refresh = True
overwrite = True
sun = True
tosave = True
resolution = 720

if resolution == 'native':
    data = [f for f in glob(path + 'cubes_10kms/native/' + '*hans.fits')]
elif resolution == 1200:
    data = [f for f in glob(path + 'cubes_10kms/exact_1200pc/*.fits')]
elif resolution == 720:
    data = [f for f in glob(path + 'cubes_10kms/exact_720pc/*.fits')]
elif resolution == 'nearest_720':
    data = [f for f in glob(path + 'cubes_10kms/nearest_aniano_720pc/*.fits')]
elif resolution == 'nearest_1200':
    data = [f for f in glob(path + 'cubes_10kms/nearest_aniano_1200pc/*.fits')]

for file in data:

    if resolution == 'native':

        long_str = file.split('/')[-1].split('hans.fits')[0]
        galaxy = long_str.split('_')[0]
        print(galaxy)

        #if not galaxy == 'NGC4579':
        #    continue

        file_pbcorr = file
        file_uncorr = file

        if not os.path.exists(path + 'products_v' + version + '/native/' + galaxy + '/'):
            os.mkdir(path + 'products_v' + version + '/native/' + galaxy + '/')
        savepath = path + 'products_v' + version + '/native/' + galaxy + '/' + galaxy + '_'

        #if os.path.exists(savepath + 'mom0_Kkms-1.pdf'):
        #    continue

    elif resolution == 1200 or resolution == 720:

        long_str = file.split('/')[-1].split('.fits')[0]
        galaxy = long_str.split('_')[0]
        print(galaxy)

        if (galaxy == 'ngc4536') | (galaxy == 'NGC4536'):
            continue

        if (galaxy == 'ngc2903') | (galaxy == 'NGC2903'):
            continue

        file_pbcorr = file
        file_uncorr = file

        if not os.path.exists(path + 'products_v' + version + '/' + str(resolution) + 'pc/' + galaxy + '/'):
            os.mkdir(path + 'products_v' + version + '/' + str(resolution) + 'pc/' + galaxy + '/')

        if resolution == 1200:
            savepath = path + 'products_v' + version + '/' + str(resolution) + 'pc/' + galaxy + '/' + galaxy + '_exact_1200pc_'
        else:
            savepath = path + 'products_v' + version + '/' + str(resolution) + 'pc/' + galaxy + '/' + galaxy + '_exact_720pc_'

        if os.path.exists(savepath + 'mom0_Kkms-1.fits'):
            continue

    elif resolution == 'nearest_720' or resolution == 'nearest_1200':

        long_str = file.split('/')[-1].split('.fits')[0]
        galaxy = long_str.split('_')[0]
        print(galaxy)

        if galaxy == 'ngc4536':
            continue

        if galaxy == 'ngc2903':
            continue

        file_pbcorr = file
        file_uncorr = file

        if resolution == 'nearest_1200':
            if not os.path.exists(path + 'products_v' + version + '/nearest_aniano_1200pc/' + galaxy + '/'):
                os.mkdir(path + 'products_v' + version + '/nearest_aniano_1200pc/' + galaxy + '/')
            savepath = path + 'products_v' + version + '/nearest_aniano_1200pc/' + galaxy + '/' + long_str + '_'
        elif resolution == 'nearest_720':
            if not os.path.exists(path + 'products_v' + version + '/nearest_aniano_720pc/' + galaxy + '/'):
                os.mkdir(path + 'products_v' + version + '/nearest_aniano_720pc/' + galaxy + '/')
            savepath = path + 'products_v' + version + '/nearest_aniano_720pc/' + galaxy + \
                       '/' + long_str + '_'

        if os.path.exists(savepath + 'mom0_Kkms-1.fits'):
            continue


    # Moment maps
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                  sun=sun, tosave=tosave, sample=sample).moment_zero(units='K km/s')
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                  sun=sun, tosave=tosave, sample=sample).moment_zero(units='M_Sun/pc^2', alpha_co=5.4)
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                  sun=sun, tosave=tosave, sample=sample).moment_zero(peak=True)
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, sample=sample).moment_1_2()
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, sample=sample).moment_1_2(moment=2)

    # Uncertainty maps
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, sample=sample).mom0_noise_maps()
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, sample=sample).mom1_2_noise_maps()

    # PVDs
    try:
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                     sun=sun, tosave=tosave, sample=sample).PVD(axis='major', find_angle=False, check_slit=False)
        if not galaxy == 'ngc2841':
            CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                          sun=sun, tosave=tosave, sample=sample).PVD(axis='minor', check_slit=False)
    except:
        pass

    # Spectra
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                  sun=sun, tosave=tosave, sample=sample).spectrum(x_axis='vel_offset')
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, sample=sample).spectrum(x_axis='velocity')
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, sample=sample).spectrum(x_axis='frequency')

    # Radial profiles
    try:
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                    sun=sun, tosave=tosave, sample=sample).radial_profile(x_units='arcsec', y_units='M_Sun pc^-2',
                                        alpha_co=5.4, table_path=None, check_aperture=False)
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                    sun=sun, tosave=tosave, sample=sample).radial_profile(x_units='kpc', y_units='M_Sun pc^-2',
                                       alpha_co=5.4, table_path=None, check_aperture=False)
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                    sun=sun, tosave=tosave, sample=sample).radial_profile(y_units='K km/s', x_units='kpc',
                                       alpha_co=5.4, table_path=None, check_aperture=False)
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                    sun=sun, tosave=tosave, sample=sample).radial_profile(y_units='K km/s', x_units='arcsec',
                                       alpha_co=5.4, table_path=None, check_aperture=False)
    except:
        pass

    #break