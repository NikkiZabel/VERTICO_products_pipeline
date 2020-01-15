from matplotlib import pyplot as plt
import aplpy as apl
from astropy import wcs
from astropy.nddata.utils import Cutout2D
from astropy.io import fits
import numpy as np
import config_figs
from sauron_colormap import register_sauron_colormap; register_sauron_colormap()


class create_images:

    def __init__(self, path, centre_y, centre_x, size=500, galaxy=None, distance=None, make_cutout=False, savepath=None, tosave=False):
        self.path = path
        self.centre_y = centre_y
        self.centre_x = centre_x
        self.size = size
        self.galaxy = galaxy
        self.distance = distance
        self.tosave = tosave
        self.savepath = savepath
        self.make_cutout = make_cutout

    def moment_zero(self):

        image = fits.open(self.path + 'moment0.fits')[0]

        f = plt.figure(figsize=(10, 8))

        if self.make_cutout:
            # make a smaller cutout of the CO emission
            w = wcs.WCS(image.header)
            cutout = Cutout2D(image.data, (self.centre_y, self.centre_x), self.size, wcs=w, mode='partial')
            image.header['CRPIX1'] = cutout.wcs.wcs.crpix[0]
            image.header['CRPIX2'] = cutout.wcs.wcs.crpix[1]
            image = fits.PrimaryHDU(cutout.data, image.header)

        # show the image in colour
        fig = apl.FITSFigure(image, figure=f)
        fig.set_theme('publication')

        # add the galaxy name in the upper right corner
        if self.galaxy:
            fig.add_label(0.8, 0.9, self.galaxy, relative=True)

        fig.show_contour(image, cmap='magma_r', levels=np.linspace(0, np.amax(image.data), 20), filled=True,
                         overlap=True)

        # axes and ticks
        fig.ticks.set_color('black')
        fig.ticks.set_length(10)
        fig.ticks.set_linewidth(2)
        fig.tick_labels.set_xformat('hh:mm:ss')
        fig.tick_labels.set_yformat('dd:mm:ss')
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        fig.axis_labels.hide_y()
        fig.ticks.set_minor_frequency(5)

        # add a colourbar
        colors = plt.contourf([[0, 0], [0, 0]],
                              levels=np.linspace(0, np.amax(image.data) + np.amax(image.data) * 0.05, 20),
                              cmap='magma_r')

        if np.amax(image.data) < 0.1:
            ticks = np.arange(0, np.amax(image.data) + 0.015, 0.005)
        if np.amax(image.data) < 0.5:
            ticks = np.arange(0, np.amax(image.data) + 0, 0.05)
        elif np.amax(image.data) < 1:
            ticks = np.arange(0, np.amax(image.data) + 0.1, 0.1)
        elif np.amax(image.data) < 2:
            ticks = np.arange(0, np.amax(image.data) + 0.2, 0.2)
        elif np.amax(image.data) < 5:
            ticks = np.arange(0, np.amax(image.data) + 0.5, 0.5)
        elif np.amax(image.data) < 10:
            ticks = np.arange(0, np.amax(image.data) + 1, 1)
        elif np.amax(image.data) < 20:
            ticks = np.arange(0, np.amax(image.data) + 1, 2)
        elif np.amax(image.data) < 100:
            ticks = np.arange(0, np.amax(image.data) + 5, 10)
        else:
            ticks = np.arange(0, np.amax(image.data) + 10, 20)

        cbar = f.colorbar(colors, ticks=ticks)
        #cbar.set_label(r'Surface brightness [Jy beam$^{-1}$ km s$^{-1}$]')
        cbar.set_label(r'Surface brightness [K]')

        # show the beam of the observations
        fig.add_beam(frame=False)  # automatically imports BMAJ, BMIN, and BPA
        fig.beam.set_edgecolor('k')
        fig.beam.set_color('k')
        fig.beam.set_borderpad(1)

        # show a scalebar
        if self.distance:
            length = np.degrees(1.e-3 / self.distance)  # length of the scalebar in degrees, corresponding to 1 kpc
            fig.add_scalebar(length=length, label='1 kpc', frame=False)
            fig.scalebar.set_linewidth(5)

        # Make sure the axis labels don't fall off the figure
        plt.tight_layout()

        if self.tosave:
            plt.savefig(self.savepath + 'mom0.png', bbox_inches='tight')

        return

    def moment_1_2(self, moment=1, vrange=100, vrange_2=100):

        if moment == 1:
            image = fits.open(self.path + 'moment1.fits')[0]
        elif moment == 2:
            image = fits.open(self.path + 'moment2.fits')[0]

        vel_array = np.loadtxt(self.path + 'vel_array.txt')
        sysvel = np.loadtxt(self.path + 'sysvel.txt')
        sysvel = (sysvel+5)//10*10

        f = plt.figure(figsize=(10, 8))

        if self.make_cutout:
            #make a smaller cutout of the CO emission
            w = wcs.WCS(image.header)
            cutout = Cutout2D(image.data, (self.centre_y, self.centre_x), self.size, wcs=w, mode='partial')
            image.header['CRPIX1'] = cutout.wcs.wcs.crpix[0]
            image.header['CRPIX2'] = cutout.wcs.wcs.crpix[1]
            image = fits.PrimaryHDU(cutout.data, image.header)

        # show the image in colour
        fig = apl.FITSFigure(image, figure=f)

        # axes and ticks
        fig.ticks.set_color('black')
        fig.ticks.set_length(10)
        fig.ticks.set_linewidth(2)
        fig.tick_labels.set_xformat('hh:mm:ss')
        fig.tick_labels.set_yformat('dd:mm:ss')
        fig.ticks.show()
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        fig.ticks.set_minor_frequency(5)

        #add a colourbar
        if moment == 2:
            fig.show_contour(image, cmap='sauron', levels=np.linspace(0, vrange_2, len(vel_array)), vmin=0, vmax=vrange_2, \
            extend='both', filled=True, overlap=True)
            colors = plt.contourf([[0, 0], [0, 0]], levels=np.linspace(0, vrange_2, len(vel_array)), cmap='sauron')

            if vrange_2 < 11:
                ticks = np.arange(0, vrange_2 + 1, 1)
            elif vrange_2 < 100:
                ticks = np.arange(0, vrange_2+10, 10)
            else:
                ticks = np.arange(0, vrange_2+20, 20)
            cbar = f.colorbar(colors, ticks=ticks)
            cbar.set_label(r'Line width [km s$^{-1}$]')

        else:
            fig.show_contour(image, cmap='sauron', levels=np.linspace(-vrange,vrange,len(vel_array)), \
                             vmin=-vrange, vmax=vrange, extend='both', filled=True, overlap=True)
            colors = plt.contourf([[0, 0], [0, 0]], levels=np.linspace(-vrange,vrange,len(vel_array)), cmap='sauron')
            if vrange < 16:
                tickarr = np.arange(-vrange, 0, 3)
            elif vrange < 60:
                tickarr = np.arange(-vrange, 0, 10)
            elif vrange < 130:
                tickarr = np.arange(-vrange, 0, 20)
            else:
                tickarr = np.arange(-vrange, 0, 40)

            ticks = np.concatenate((tickarr, [0], abs(tickarr)))
            cbar = f.colorbar(colors, ticks=ticks)
            cbar.set_label(r'Velocity [km s$^{-1}$]')

        # show the beam of the observations
        fig.add_beam(frame=False)  # automatically imports BMAJ, BMIN, and BPA
        fig.beam.set_edgecolor('k')
        fig.beam.set_color('k')
        fig.beam.set_borderpad(1)

        # show a scalebar
        if self.distance:
            length = np.degrees(1.e-3 / self.distance)  # length of the scalebar in degrees, corresponding to 1 kpc
            fig.add_scalebar(length=length, label='1 kpc', frame=False)
            fig.scalebar.set_linewidth(5)

        #Make sure the axis labels don't fall off the figure
        plt.tight_layout()

        if self.tosave:
            if moment == 2:
                plt.savefig(self.savepath+'moment2.png', bbox_inches='tight',)
            else:
                plt.savefig(self.savepath+'moment1.png', bbox_inches='tight')
        return