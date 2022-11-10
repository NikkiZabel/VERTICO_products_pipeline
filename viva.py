from image_moments import CreateImages
from create_moments import *
import warnings; warnings.filterwarnings("ignore")
import os
from glob import glob

sample = 'viva'
version = '1_1'
path = '/home/nikki/Documents/Data/VERTICO/VIVA/Reprojected_15_arcsec_latest/'
table_path = '/home/nikki/Documents/Data/VERTICO/VERTICO_master.fits'

refresh = True
overwrite = True
sun = True
tosave = True
redo_clip = True

data = [f for f in glob(path + '*.fits')]

for file in data:

    # Split off the galaxy's name
    galnames = file.split('/')[8].split('_cube')[0].split('_')

    # Sometimes a cube contains two galaxies, loop over both to create individual products
    for i in range(len(galnames)):

        # Add 'ngc' in front of the second galaxy name, as it is only listed as a number
        if i > 0:
            if galnames[i] == 'viva':
                break
            galaxy = 'NGC' + galnames[i]
        else:
            galaxy = galnames[i]

        if not galaxy == 'NGC4606':
            continue

        print(galaxy)

        if galaxy == 'NGC4293':
            continue
        if galaxy == 'NGC4419':
            continue
        #if galaxy == 'NGC4533':
        #    continue
        #if galaxy == 'NGC4606':
        #    continue

        # There is only one cube, so pb un/corrected files are the same
        file_pbcorr = file
        file_uncorr = file

        if not os.path.exists(path + 'products_v' + version + '/' + galaxy + '/'):
            os.mkdir(path + 'products_v' + version + '/' + galaxy + '/')
        #else:
        #    continue
        savepath = path + 'products_v' + version + '/' + galaxy + '/' + galaxy + '_viva_hi_15as_np_round_reproj_'

        #if os.path.exists(savepath + 'mom0_Jyb-1kms-1.pdf'):
        #    continue

        # Moment maps
        #redo_clip = True
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                      sun=sun, tosave=tosave, sample=sample, redo_clip=redo_clip).moment_zero(units='K km/s')
        redo_clip = False
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

        # Spectra
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                     sun=sun, tosave=tosave, sample=sample, redo_clip=redo_clip).spectrum(x_axis='vel_offset')
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                     sun=sun, tosave=tosave, sample=sample, redo_clip=redo_clip).spectrum(x_axis='velocity')
        CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                     sun=sun, tosave=tosave, sample=sample, redo_clip=redo_clip).spectrum(x_axis='frequency')

        # PVDs
        if not galaxy == 'NGC4579':
            CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                         sun=sun, tosave=tosave, sample=sample).PVD(axis='major', find_angle=False, check_slit=False)
        if not galaxy == 'NGC4450' or galaxy == 'NGC4536':
            CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                         sun=sun, tosave=tosave, sample=sample).PVD(axis='minor', check_slit=False)
        
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