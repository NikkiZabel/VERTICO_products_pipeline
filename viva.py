from image_moments import CreateImages
from create_moments import *
import warnings; warnings.filterwarnings("ignore")
import os
from glob import glob

sample = 'viva'
version = '1_0'
path = '/home/nikki/Documents/Data/VERTICO/VIVA/Reprojected_15_arcsec_new/'
table_path = '/home/nikki/Documents/Data/VERTICO/VERTICO_master.fits'

refresh = True
overwrite = True
sun = True
tosave = True
redo_clip = False

data = [f for f in glob(path + '*.fits')]

for file in data:

    # Split off the galaxy's name
    galnames = file.split('/')[8].split('_cube')[0].split('_')

    # Sometimes a cube contains two galaxies, loop over both to create individual products
    for i in range(len(galnames)):

        # Add 'ngc' in front of the second galaxy name, as it is only listed as a number
        if i > 0:
            galaxy = 'ngc' + galnames[i]
        else:
            galaxy = galnames[i]

        print(galaxy)
        if galaxy == 'ngc4293':
            continue

        # There is only one cube, so pb un/corrected files are the same
        file_pbcorr = file
        file_uncorr = file

        if not os.path.exists(path + 'products_v' + version + '/sun18_method/' + galaxy + '/'):
            os.mkdir(path + 'products_v' + version + '/sun18_method/' + galaxy + '/')
        else:
            continue
        savepath = path + 'products_v' + version + '/sun18_method/' + galaxy + '/' + galaxy + '_reprojected_15_arcsec_'

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
                     sun=sun, tosave=tosave, redo_clip=redo_clip).mom0_noise_maps()
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                     sun=sun, tosave=tosave, redo_clip=redo_clip).mom1_2_noise_maps()

        # PVDs
        if not galaxy == 'ngc4579':
            CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                         sun=sun, tosave=tosave).PVD(axis='major', find_angle=False, check_slit=False)
        if not galaxy == 'ngc4450' or galaxy == 'ngc4536':
            CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                         sun=sun, tosave=tosave).PVD(axis='minor', check_slit=False)

        # Spectra
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                     sun=sun, tosave=tosave, redo_clip=redo_clip).spectrum(x_axis='vel_offset')
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                     sun=sun, tosave=tosave, redo_clip=redo_clip).spectrum(x_axis='velocity')
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                     sun=sun, tosave=tosave, redo_clip=redo_clip, sample=sample).spectrum(x_axis='frequency')

        # Radial profiles
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                    sun=sun, tosave=tosave, sample=sample, redo_clip=redo_clip).radial_profile(x_units='arcsec', y_units='M_Sun pc^-2',
                                        alpha_co=6.25, table_path=table_path, check_aperture=False)
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                    sun=sun, tosave=tosave, sample=sample, redo_clip=redo_clip).radial_profile(x_units='kpc', y_units='M_Sun pc^-2',
                                       alpha_co=6.25, table_path=table_path, check_aperture=False)
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                    sun=sun, tosave=tosave, sample=sample, redo_clip=redo_clip).radial_profile(y_units='K km/s', x_units='kpc',
                                       alpha_co=6.25, table_path=table_path, check_aperture=False)
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                     sun=sun, tosave=tosave, sample=sample, redo_clip=redo_clip).radial_profile(y_units='K km/s', x_units='arcsec',
                                                                           alpha_co=6.25, table_path=table_path,
                                                                           check_aperture=False)
