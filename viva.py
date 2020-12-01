from image_moments import CreateImages
from create_moments import *
import warnings; warnings.filterwarnings("ignore")
import os
from glob import glob

sample = 'viva'
version = '1_0'
path = '/home/nikki/Documents/Data/VERTICO/VIVA/Reprojected_15_arcsec/'

refresh = True
overwrite = True
sun = True
tosave = True

data = [f for f in glob(path + 'cubes/' + '*.fits')]

for file in data:

    # Split off the galaxy's name
    galnames = file.split('/')[9].split('_cube')[0].split('_')

    # Sometimes a cube contains two galaxies, loop over both to create individual products
    for i in range(len(galnames)):

        # Add 'ngc' in front of the second galaxy name, as it is only listed as a number
        if i > 0:
            galaxy = 'ngc' + galnames[i]
        else:
            galaxy = galnames[i]

        print(galaxy)

        # There is only one cube, so pb un/corrected files are the same
        file_pbcorr = file
        file_uncorr = file

        if not os.path.exists(path + 'products_v' + version + '/sun18_method/' + galaxy + '/'):
            os.mkdir(path + 'products_v' + version + '/sun18_method/' + galaxy + '/')
        #else:
        #    continue
        savepath = path + 'products_v' + version + '/sun18_method/' + galaxy + '/'

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

        # Spectra
        try:
            CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                         sun=sun, tosave=tosave).spectrum(x_axis='vel_offset')
            CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                         sun=sun, tosave=tosave).spectrum(x_axis='velocity')
            CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                         sun=sun, tosave=tosave).spectrum(x_axis='frequency')
        except:
            pass

        #break









