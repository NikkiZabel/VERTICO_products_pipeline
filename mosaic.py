from matplotlib import pyplot as plt
import os
from glob import glob
import numpy as np
from matplotlib.gridspec import GridSpec
from astropy.io import fits
from sauron_colormap import register_sauron_colormap; register_sauron_colormap()
from gal_params import parameters
import yt; cmap_name = 'RED TEMPERATURE_r'

# Read in the Smooth version of the moment 0 map for each galaxy
path = '/home/nikki/Documents/Data/VERTICO/Products/'
files_mom0 = [y for x in os.walk(path) for y in glob(os.path.join(x[0], '*/Sun_method/*moment0_K.fits'))]
files_mom1 = [y for x in os.walk(path) for y in glob(os.path.join(x[0], '*/Sun_method/*moment1.fits'))]
files_mom2 = [y for x in os.walk(path) for y in glob(os.path.join(x[0], '*/Sun_method/*moment2.fits'))]
files_vel = [y for x in os.walk(path) for y in glob(os.path.join(x[0], '*/Sun_method/*spectrum.csv'))]
files_peak = [y for x in os.walk(path) for y in glob(os.path.join(x[0], '*/Sun_method/*peak_temperature.fits'))]

def make_square(img):
    shape_diff = np.max(img.shape) - np.min(img.shape)
    square_img = np.zeros((np.max(img.shape), np.max(img.shape)))

    if img.shape[0] > img.shape[1]:
        first = True
    else:
        first = False
    if first:
        square_img[:, int(shape_diff/2):int(shape_diff/2+img.shape[1])] = img
    else:
        square_img[int(shape_diff/2):int(shape_diff/2+img.shape[0]), :] = img
    return square_img

count = 0

moment = 'peak'

f, axarr = plt.subplots(7, 7, figsize=(10, 10))
gs = GridSpec(7, 7, figure=f)

for i in range(axarr.shape[0]):
    for j in range(axarr.shape[1]):

        galaxy = files_mom0[count].split('/')[7]

        name, vrange, vrange2, cliplevel, stokes, start, stop, sysvel_offset, angle, \
        full_width, distance, nchan_low, cliplevel_low, nchan_high, cliplevel_high, prune_by_npix, \
        prune_by_fracbeam, expand_by_fracbeam, expand_by_npix, expand_by_nchan, inclination, eccentricity, \
        figsize = parameters(galaxy)

        axarr[i, j].axis('off')

        if (j == 3 and i == 3) or (j == 3 and i == 3): #or (j == 2 and i == 5) or (j == 2 and i == 6):
            continue

        if moment == 0:
            img = make_square(fits.open(files_mom0[count])[0].data)
            levels = np.linspace(np.std(img) / 1e5, np.max(img), 50)
            axarr[i, j].contourf(img, cmap='magma_r', levels=levels)
        elif moment == 1:
            img = make_square(fits.open(files_mom1[count])[0].data)
            img[img == 0] = np.nan
            vel_array = np.genfromtxt(files_vel[count], delimiter=',', skip_header=1)[:, 2]
            try:
                axarr[i, j].contourf(img, cmap='sauron', levels=np.linspace(-vrange, vrange, len(vel_array)),
                                     vmin=-vrange, vmax=vrange)
            except:
                pass
        elif moment == 2:
            img = make_square(fits.open(files_mom2[count])[0].data)
            img[img == 0] = np.nan
            vel_array = np.genfromtxt(files_vel[count], delimiter=',', skip_header=1)[:, 2]
            try:
                axarr[i, j].contourf(img, cmap='sauron', levels=np.linspace(0, vrange2, len(vel_array)), vmin=0, vmax=vrange2)
            except:
                pass
        elif moment == 'peak':
            img = make_square(fits.open(files_mom0[count])[0].data)
            levels = np.linspace(np.std(img) / 1e5, np.max(img), 50)
            axarr[i, j].contourf(img, cmap='RED TEMPERATURE_r', levels=levels)

        count += 1

ax1 = f.add_subplot(gs[3:4, 3:4])
ax1.axis('off')
logo = plt.imread('/home/nikki/Documents/Data/VERTICO/Z_materials/concept_logo.png')
ax1.imshow(logo)

plt.tight_layout()