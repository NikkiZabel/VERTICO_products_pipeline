from image_moments import CreateImages
from create_moments import *
import warnings; warnings.filterwarnings("ignore")
import os
from glob import glob

sample = 'things'
version = '1_0'
path = '/home/nikki/Documents/Data/VERTICO/THINGS/'

refresh = True
overwrite = True
sun = True
tosave = True
redo_clip = False

data = [f for f in glob(path + '*.fits')]

for file in data:

    # Split off the galaxy's name
    galaxy = file.split('/')[7].split('_')[0]

    print(galaxy)

    # There is only one cube, so pb un/corrected files are the same
    file_pbcorr = file
    file_uncorr = file

    if not os.path.exists(path + 'products_v' + version + '/' + galaxy + '/'):
        os.mkdir(path + 'products_v' + version + '/' + galaxy + '/')
    else:
        continue
    savepath = path + 'products_v' + version + '/' + galaxy + '/' + galaxy + '_'

    # Moment maps

    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                  sun=sun, tosave=tosave, sample=sample, redo_clip=redo_clip).moment_zero(units='K km/s')
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                  sun=sun, tosave=tosave, sample=sample, redo_clip=redo_clip).moment_zero()
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                  sun=sun, tosave=tosave, sample=sample, redo_clip=redo_clip).moment_zero(peak=True)
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, sample=sample, redo_clip=redo_clip).moment_1_2()
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, sample=sample, redo_clip=redo_clip).moment_1_2(moment=2)

    # Uncertainty maps
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).mom0_noise_maps()
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).mom1_2_noise_maps()

    # PVDs
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, sample=sample).PVD(axis='major', find_angle=False, check_slit=False)
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                     sun=sun, tosave=tosave, sample=sample).PVD(axis='minor', check_slit=False)

    # Spectra
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).spectrum(x_axis='vel_offset')
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).spectrum(x_axis='velocity')
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).spectrum(x_axis='frequency')

    # Radial profiles
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                sun=sun, tosave=tosave, sample=sample, redo_clip=redo_clip).radial_profile(x_units='arcsec', y_units='M_Sun pc^-2',
                                    alpha_co=6.25, table_path=None, check_aperture=False)
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                sun=sun, tosave=tosave, sample=sample, redo_clip=redo_clip).radial_profile(x_units='kpc', y_units='M_Sun pc^-2',
                                   alpha_co=6.25, table_path=None, check_aperture=False)
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                sun=sun, tosave=tosave, sample=sample, redo_clip=redo_clip).radial_profile(y_units='K km/s', x_units='kpc',
                                   alpha_co=6.25, table_path=None, check_aperture=False)
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave, sample=sample, redo_clip=redo_clip).radial_profile(y_units='K km/s', x_units='arcsec',
                                                                       alpha_co=6.25, table_path=None,
                                                                       check_aperture=False)

