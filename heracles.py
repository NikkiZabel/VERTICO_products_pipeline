from image_moments import CreateImages
from create_moments import *
import warnings; warnings.filterwarnings("ignore")
import os
from glob import glob
from astropy.io import fits

sample = 'heracles'
version = '1_0'
path = '/home/nikki/Documents/Data/VERTICO/heracles/'

refresh = True
overwrite = True
sun = True
tosave = True
resolution = 1200

def get_inc_pa(galaxy, ):
    table = fits.open(path + 'heracles_sdss_r_properties.fits')[1]
    gal_name_table = table.data['Galaxy']
    gal_num_table = np.array([n.split('C')[1] for n in gal_name_table])

    galaxy_num = galaxy.split('c')[1]
    incl = table.data['inclination'][gal_num_table == galaxy_num]
    pa = table.data['pa'][gal_num_table == galaxy_num]

    return incl[0], pa[0]

if resolution == 1200:
    data = [f for f in glob(path + '*1200pc.fits')]
elif resolution == 720:
    data = [f for f in glob(path + '*720pc.fits')]
elif resolution == 'nearest_720':
    data = [f for f in glob(path + 'smoothed_cubes_nearest_720pc/' + '*.fits')]
elif resolution == 'nearest_1200':
    data = [f for f in glob(path + 'smoothed_cubes_nearest_1200pc/' + '*.fits')]

for file in data:

    if resolution == 1200 or resolution == 720:

        galaxy = file.split('/')[7].split('_')[0]
        print(galaxy)

        file_pbcorr = file
        file_uncorr = file

        inc, pa = get_inc_pa(galaxy)

        if not os.path.exists(path + 'products_v' + version + '/' + str(resolution) + 'pc/sun18_method/' + galaxy + '/'):
            os.mkdir(path + 'products_v' + version + '/' + str(resolution) + 'pc/sun18_method/' + galaxy + '/')
        #else:
        #    continue
        savepath = path + 'products_v' + version + '/' + str(resolution) + 'pc/sun18_method/' + galaxy + '/'

    elif resolution == 'nearest_720' or resolution == 'nearest_1200':

        galaxy = file.split('/')[8].split('_')[0]
        print(galaxy)

        file_pbcorr = file
        file_uncorr = file

        if not os.path.exists(path + 'products_v' + version + '/smoothed_cubes_' + resolution + 'pc/sun18_method/' + galaxy + '/'):
            os.mkdir(path + 'products_v' + version + '/smoothed_cubes_' + resolution + 'pc/sun18_method/' + galaxy + '/')
        #else:
        #    continue
        savepath = path + 'products_v' + version + '/smoothed_cubes_' + resolution + 'pc/sun18_method/' + galaxy + '/'

    #CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
    #              sun=sun, tosave=tosave, sample=sample).moment_zero(units='K km/s')
    #CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
    #              sun=sun, tosave=tosave, sample=sample).moment_zero()
    #CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
    #              sun=sun, tosave=tosave, sample=sample).moment_zero(peak=True)
    #CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
    #             sun=sun, tosave=tosave, sample=sample).moment_1_2()
    #CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
    #             sun=sun, tosave=tosave, sample=sample).moment_1_2(moment=2)

    break