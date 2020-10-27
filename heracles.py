from image_moments import CreateImages
from create_moments import *
import warnings; warnings.filterwarnings("ignore")
import os
from glob import glob

version = '1_0'
path = '/home/nikki/Documents/Data/VERTICO/heracles/'

refresh = True
overwrite = True
sun = True
tosave = True
resolution = 1200

if resolution == 1200:
    data = [f for f in glob(path + '*1200pc.fits')]
elif resolution == 720:
    data = [f for f in glob(path + '*1200pc.fits')]

for file in data:
    galaxy = file.split('/')[7].split('_')[0]
    print(galaxy)

    file_pbcorr = file
    file_uncorr = file

    from matplotlib import pyplot as plt
    #x = fits.open(file_pbcorr)[0]
    #spectrum = np.nansum(x.data, axis=(1, 2))
    #print(np.where(spectrum != 0)[0])
    #break

    if not os.path.exists(path + 'products_v' + version + '/' + str(resolution) + 'pc/sun18_method/' + galaxy + '/'):
        os.mkdir(path + 'products_v' + version + '/' + str(resolution) + 'pc/sun18_method/' + galaxy + '/')
    savepath = path + 'products_v' + version + '/' + str(resolution) + 'pc/sun18_method/' + galaxy + '/'

    #CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
    #              sun=sun, tosave=tosave).moment_zero(units='K km/s')
    #CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
    #              sun=sun, tosave=tosave).moment_zero()
    #CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
    #              sun=sun, tosave=tosave).moment_zero(peak=True)
    #CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
    #             sun=sun, tosave=tosave).moment_1_2()
    CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
                 sun=sun, tosave=tosave).moment_1_2(moment=2)

    break