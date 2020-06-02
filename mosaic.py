from matplotlib import pyplot as plt
import os
from glob import glob
import numpy as np
from matplotlib.gridspec import GridSpec
from astropy.io import fits

# Read in the Smooth version of the moment 0 map for each galaxy
path = '/home/nikki/Documents/Data/VERTICO/Products/'
files = [y for x in os.walk(path) for y in glob(os.path.join(x[0], '*/Smooth_method/*moment0_K.fits'))]

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

f, axarr = plt.subplots(7, 7, figsize=(10, 10))
gs = GridSpec(7, 7, figure=f)

for i in range(axarr.shape[0]):
    for j in range(axarr.shape[1]):

        axarr[i, j].axis('off')

        if (j == 3 and i == 3) or (j == 3 and i == 3): #or (j == 2 and i == 5) or (j == 2 and i == 6):
            continue

        img = make_square(fits.open(files[count])[0].data)

        levels = np.linspace(np.std(img)/1e5, np.max(img), 50)
        axarr[i, j].contourf(img, cmap='magma_r', levels=levels)

        count += 1

ax1 = f.add_subplot(gs[3:4, 3:4])
ax1.axis('off')
logo = plt.imread('/home/nikki/Documents/Data/VERTICO/Z_materials/concept_logo.png')
ax1.imshow(logo)

plt.tight_layout()