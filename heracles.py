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
resolution = 'native'

if resolution == '1200_nyquist':
    data = [f for f in glob(path + 'cubes_10kms/1200pc/*_np_round_k.fits')]
elif resolution == '720_nyquist':
    data = [f for f in glob(path + 'cubes_10kms/720pc/*_np_round_k.fits')]
elif resolution == 'native':
    data = [f for f in glob(path + 'cubes_10kms/' + resolution + '/*hans.fits')]
else:
    data = [f for f in glob(path + 'cubes_10kms/' + resolution + '/*' + resolution + '_round_k.fits')]

for file in data:

    long_str = file.split('/')[-1].split('hans.fits')[0]
    galaxy = long_str.split('_')[0]
    print(galaxy)

    if (galaxy == 'NGC5474') or (galaxy == 'NGC4236'):
        continue

    if not galaxy == 'NGC2903':
        continue

    file_pbcorr = file
    file_uncorr = file

    if resolution == 'native':
        if not os.path.exists(path + 'products_v' + version + '/' + resolution + '/' + galaxy + '/'):
            os.mkdir(path + 'products_v' + version + '/' + resolution + '/' + galaxy + '/')
        savepath = path + 'products_v' + version + '/' + resolution + '/' + galaxy + '/' + galaxy + '_heracles_co21_round_'
    elif resolution == '1200_nyquist':
        if not os.path.exists(path + 'products_v' + version + '/1200pc/' + galaxy + '/'):
            os.mkdir(path + 'products_v' + version + '/1200pc/' + galaxy + '/')
        savepath = path + 'products_v' + version + '/1200pc/' + galaxy + '/' + galaxy + '_heracles_co21_1200pc_np_round_'
    elif resolution == '720_nyquist':
        if not os.path.exists(path + 'products_v' + version + '/720pc/' + galaxy + '/'):
            os.mkdir(path + 'products_v' + version + '/720pc/' + galaxy + '/')
        savepath = path + 'products_v' + version + '/720pc/' + galaxy + '/' + galaxy + '_heracles_co21_720pc_np_round_'
    else:
        if not os.path.exists(path + 'products_v' + version + '/' + resolution + '/' + galaxy + '/'):
            os.mkdir(path + 'products_v' + version + '/' + resolution + '/' + galaxy + '/')
        savepath = path + 'products_v' + version + '/' + resolution + '/' + galaxy + '/' + galaxy + '_heracles_co21_' + resolution + '_round_'

    #if os.path.exists(savepath + 'mom0_Kkms-1.pdf'):
    #    continue
    '''
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
    #try:
    if not galaxy == 'NGC2841':
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                     sun=sun, tosave=tosave, sample=sample).PVD(axis='major', find_angle=False, check_slit=False)
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                     sun=sun, tosave=tosave, sample=sample).PVD(axis='minor', check_slit=False)
    #except:
    #    pass
    
    # Spectra
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                  sun=sun, tosave=tosave, sample=sample).spectrum(x_axis='vel_offset')
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, sample=sample).spectrum(x_axis='velocity')
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, sample=sample).spectrum(x_axis='frequency')
    '''
    # Radial profiles
    #try:
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
