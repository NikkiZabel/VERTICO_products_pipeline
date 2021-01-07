from image_moments import CreateImages
from create_moments import *
import warnings; warnings.filterwarnings("ignore")
import os
from glob import glob

sample = 'heracles'
version = '1_0'
path = '/home/nikki/Documents/Data/VERTICO/heracles/'

refresh = True
overwrite = True
sun = True
tosave = True
resolution = 'native'

if resolution == 'native':
    data = [f for f in glob(path + 'native/' + '*hans.fits.gz')]
elif resolution == 1200:
    data = [f for f in glob(path + '*1200pc.fits')]
elif resolution == 720:
    data = [f for f in glob(path + '*720pc.fits')]
elif resolution == 'nearest_720':
    data = [f for f in glob(path + 'smoothed_cubes_nearest_720pc/' + '*.fits')]
elif resolution == 'nearest_1200':
    data = [f for f in glob(path + 'smoothed_cubes_nearest_1200pc/' + '*.fits')]

for file in data:

    if resolution == 'native':

        long_str = file.split('/')[-1].split('hans.fits.gz')[0]
        galaxy = long_str.split('_')[0]
        print(galaxy)

        file_pbcorr = file
        file_uncorr = file

        if not os.path.exists(path + 'products_v' + version + '/native/sun18_method/' + galaxy + '/'):
            os.mkdir(path + 'products_v' + version + '/native/sun18_method/' + galaxy + '/')
        savepath = path + 'products_v' + version + '/native/sun18_method/' + galaxy + '/' + long_str + '_'

    elif resolution == 1200 or resolution == 720:

        long_str = file.split('/')[-1].split('.fits')[0]
        galaxy = long_str.split('_')[0]
        print(galaxy)

        if galaxy == 'ngc4536':
            continue

        if galaxy == 'ngc2903':
            continue

        file_pbcorr = file
        file_uncorr = file

        if not os.path.exists(path + 'products_v' + version + '/' + str(resolution) + 'pc/sun18_method/' + galaxy + '/'):
            os.mkdir(path + 'products_v' + version + '/' + str(resolution) + 'pc/sun18_method/' + galaxy + '/')

        if resolution == 1200:
            savepath = path + 'products_v' + version + '/' + str(resolution) + 'pc/sun18_method/' + galaxy + '/' + galaxy + '_exact_1200pc_'
        else:
            savepath = path + 'products_v' + version + '/' + str(resolution) + 'pc/sun18_method/' + galaxy + '/' + galaxy + '_exact_720pc_'

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
            if not os.path.exists(path + 'products_v' + version + '/nearest_aniano_1200pc/sun18_method/' + galaxy + '/'):
                os.mkdir(path + 'products_v' + version + '/nearest_aniano_1200pc/sun18_method/' + galaxy + '/')
            savepath = path + 'products_v' + version + '/nearest_aniano_1200pc/sun18_method/' + galaxy + '/' + long_str + '_'
        elif resolution == 'nearest_720':
            if not os.path.exists(path + 'products_v' + version + '/nearest_aniano_720pc/sun18_method/' + galaxy + '/'):
                os.mkdir(path + 'products_v' + version + '/nearest_aniano_720pc/sun18_method/' + galaxy + '/')
            savepath = path + 'products_v' + version + '/nearest_aniano_720pc/sun18_method/' + galaxy + \
                       '/' + long_str + '_'

        if os.path.exists(savepath + 'mom0_Kkms-1.fits'):
            continue

    # Moment maps
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                  sun=sun, tosave=tosave, sample=sample).moment_zero(units='K km/s')
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                  sun=sun, tosave=tosave, sample=sample).moment_zero()
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                  sun=sun, tosave=tosave, sample=sample).moment_zero(peak=True)
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, sample=sample).moment_1_2()
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, sample=sample).moment_1_2(moment=2)

    # Uncertainty maps
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave).mom0_noise_maps()
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave).mom1_2_noise_maps()

    # PVDs
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave).PVD(axis='major', find_angle=False, check_slit=False)
    if not galaxy == 'ngc2841':
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                      sun=sun, tosave=tosave).PVD(axis='minor', check_slit=False)

    # Spectra
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                  sun=sun, tosave=tosave).spectrum(x_axis='vel_offset')
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave).spectrum(x_axis='velocity')
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave).spectrum(x_axis='frequency')

    # Radial profiles
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                sun=sun, tosave=tosave, sample=sample).radial_profile(x_units='arcsec', y_units='M_Sun pc^-2',
                                    alpha_co=6.25, table_path=None, check_aperture=False)
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                sun=sun, tosave=tosave, sample=sample).radial_profile(x_units='kpc', y_units='M_Sun pc^-2',
                                   alpha_co=6.25, table_path=None, check_aperture=False)
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                sun=sun, tosave=tosave, sample=sample).radial_profile(y_units='K km/s', x_units='kpc',
                                   alpha_co=6.25, table_path=None, check_aperture=False)
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                sun=sun, tosave=tosave, sample=sample).radial_profile(y_units='K km/s', x_units='arcsec',
                                   alpha_co=6.25, table_path=None, check_aperture=False)

    #break