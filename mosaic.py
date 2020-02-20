from matplotlib import pyplot as plt
import os
import numpy as np
from matplotlib.gridspec import GridSpec

path = '/home/nikki/Documents/Data/VERTICO/Mosaic/'
files = os.listdir(path)
count = 0

f, axarr = plt.subplots(5, 6, figsize=(10, 10))
gs = GridSpec(5, 6, figure=f)

for i in range(axarr.shape[0]):
    for j in range(axarr.shape[1]):

        axarr[i, j].axis('off')

        if (i == 2 and j == 1) or (i == 2 and j == 2) or (i == 2 and j == 3) or (i == 2 and j == 4):
            continue

        img = plt.imread(path + files[count])

        img = img[:, ~np.all(img==1, axis=(0, 2)), :]
        img = img[~np.all(img==1, axis=(1, 2)), :, :]

        axarr[i, j].imshow(img)

        count += 1

ax1 = f.add_subplot(gs[2, 1:5])
ax1.axis('off')
logo = plt.imread('/home/nikki/Documents/Data/VERTICO/logo-full-blackblue.png')
ax1.imshow(logo)

plt.tight_layout()