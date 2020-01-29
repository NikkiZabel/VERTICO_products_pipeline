from matplotlib import pyplot as plt
import aplpy as apl
from astropy import wcs
from astropy.nddata.utils import Cutout2D
from astropy.io import fits
import numpy as np
from scipy import ndimage
import config_figs
from sauron_colormap import register_sauron_colormap; register_sauron_colormap()
from targets import galaxies
from create_moments import moment_maps

class create_images:

    def __init__(self, galname, path_pbcorr, path_uncorr, savepath=None, refresh=False, overwrite=True,
                 make_cutout=False, sun=True, tosave=False):
        self.galaxy = galaxies(galname)
        self.path_pbcorr = path_pbcorr
        self.path_uncorr = path_uncorr
        self.refresh = refresh
        self.overwrite = overwrite
        self.savepath = savepath
        self.make_cutout = make_cutout
        self.sun = sun
        self.tosave = tosave

    def moment_zero(self):

        if self.refresh:
            if self.overwrite:
                _, image, _, _, _ = moment_maps(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                                tosave=True).calc_moms()
            else:
                _, image, _, _, _ = moment_maps(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                                tosave=False).calc_moms()
        else:
            image = fits.open(self.path + 'moment0.fits')[0]

        f = plt.figure(figsize=(12, 7))

        if self.make_cutout:
            # make a smaller cutout of the CO emission
            w = wcs.WCS(image.header)
            cutout = Cutout2D(image.data, (self.galaxy.centre_y, self.galaxy.centre_x), self.galaxy.size, wcs=w,
                              mode='partial')
            image.header['CRPIX1'] = cutout.wcs.wcs.crpix[0]
            image.header['CRPIX2'] = cutout.wcs.wcs.crpix[1]
            image = fits.PrimaryHDU(cutout.data, image.header)

        # show the image in colour
        fig = apl.FITSFigure(image, figure=f)
        fig.set_theme('publication')

        # add the galaxy name in the upper right corner
        fig.add_label(0.88, 0.9, self.galaxy.name, relative=True)

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
        length = np.degrees(1.e-3 / self.galaxy.distance)  # length of the scalebar in degrees, corresponding to 1 kpc
        fig.add_scalebar(length=length, label='1 kpc', frame=False)
        fig.scalebar.set_linewidth(5)

        plt.tight_layout()

        if self.tosave:
            plt.savefig(self.savepath + 'moment0.pdf', bbox_inches='tight')

        return

    def moment_1_2(self, moment=1):

        if moment == 1:
            if self.refresh:
                if self.overwrite:
                    _, _, image, _, sysvel = moment_maps(self.galaxy.name, self.path_pbcorr, self.path_uncorr,
                                                         sun=self.sun, tosave=True).calc_moms()
                else:
                    _, _, image, _, sysvel = moment_maps(self.galaxy.name, self.path_pbcorr, self.path_uncorr,
                                                         sun=self.sun, tosave=False).calc_moms()
            else:
                image = fits.open(self.path + 'moment1.fits')[0]

        elif moment == 2:
            if self.refresh:
                if self.overwrite:
                    _, _, _, image, sysvel = moment_maps(self.galaxy.name, self.path_pbcorr, self.path_uncorr,
                                                         sun=self.sun, tosave=True).calc_moms()
                else:
                    _, _, _, image, sysvel = moment_maps(self.galaxy.name, self.path_pbcorr, self.path_uncorr,
                                                         sun=self.sun, tosave=False).calc_moms()
            else:
                image = fits.open(self.path + 'moment2.fits')[0]

        cube_pbcorr, cube_uncorr = moment_maps(self.galaxy.name, self.path_pbcorr, self.path_uncorr,
                                               sun=self.sun, tosave=False).readfits()
        vel_array, _, _ = moment_maps(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                      tosave=False).create_vel_array(cube_pbcorr)

        sysvel = (sysvel + 5) // 10 * 10

        f = plt.figure(figsize=(12, 7))

        if self.make_cutout:
            #make a smaller cutout of the CO emission
            w = wcs.WCS(image.header)
            cutout = Cutout2D(image.data, (self.galaxy.centre_y, self.galaxy.centre_x), self.galaxy.size, wcs=w,
                              mode='partial')
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
            fig.show_contour(image, cmap='sauron', levels=np.linspace(0, self.galaxy.vrange2, len(vel_array)), vmin=0,
                             vmax=self.galaxy.vrange2, extend='both', filled=True, overlap=True)
            colors = plt.contourf([[0, 0], [0, 0]], levels=np.linspace(0, self.galaxy.vrange2, len(vel_array)),
                                  cmap='sauron')

            if self.galaxy.vrange2 < 11:
                ticks = np.arange(0, self.galaxy.vrange2 + 1, 1)
            elif self.galaxy.vrange2 < 100:
                ticks = np.arange(0, self.galaxy.vrange2 + 10, 10)
            else:
                ticks = np.arange(0, self.galaxy.vrange2 + 20, 20)
            cbar = f.colorbar(colors, ticks=ticks)
            cbar.set_label(r'Line width [km s$^{-1}$]')

        else:
            fig.show_contour(image, cmap='sauron', levels=np.linspace(-self.galaxy.vrange, self.galaxy.vrange,
                len(vel_array)), vmin=-self.galaxy.vrange, vmax=self.galaxy.vrange, extend='both', filled=True,
                             overlap=True)
            colors = plt.contourf([[0, 0], [0, 0]], levels=np.linspace(-self.galaxy.vrange, self.galaxy.vrange,
                                                                       len(vel_array)), cmap='sauron')
            if self.galaxy.vrange < 16:
                tickarr = np.arange(-self.galaxy.vrange, 0, 3)
            elif self.galaxy.vrange < 60:
                tickarr = np.arange(-self.galaxy.vrange, 0, 10)
            elif self.galaxy.vrange < 130:
                tickarr = np.arange(-self.galaxy.vrange, 0, 20)
            else:
                tickarr = np.arange(-self.galaxy.vrange, 0, 40)

            ticks = np.concatenate((tickarr, [0], abs(tickarr)))
            cbar = f.colorbar(colors, ticks=ticks)
            cbar.set_label(r'Velocity [km s$^{-1}$]')

        # show the beam of the observations
        fig.add_beam(frame=False)  # automatically imports BMAJ, BMIN, and BPA
        fig.beam.set_edgecolor('k')
        fig.beam.set_color('k')
        fig.beam.set_borderpad(1)

        # show a scalebar
        length = np.degrees(1.e-3 / self.galaxy.distance)  # length of the scalebar in degrees, corresponding to 1 kpc
        fig.add_scalebar(length=length, label='1 kpc', frame=False)
        fig.scalebar.set_linewidth(5)

        #Make sure the axis labels don't fall off the figure
        plt.tight_layout()

        if self.tosave:
            if moment == 2:
                plt.savefig(self.savepath+'moment2.pdf', bbox_inches='tight',)
            else:
                plt.savefig(self.savepath+'moment1.pdf', bbox_inches='tight')

        return

    def PVD(self, findcentre=False, find_velcentre=False, full_width=False):
        """
        Plot a position-velocity diagram of the galaxy
        :param mom0: moment 0 image of the galaxy
        :param mom1: moment 1 image of the galaxy
        :param angle: angle by which the images should be rotated to obtain a horizontal projection of the galaxy
        :param centre: the vertical offset around which the emission is centred in the rotated image
        """

        #matplotlib.rcParams['axes.linewidth'] = 0.5

        if self.refresh:
            if self.overwrite:
                clipped_cube, _, _, _, _ = moment_maps(self.galaxy.name, self.path_pbcorr, self.path_uncorr,
                                                       sun=self.sun, tosave=True).calc_moms()
            else:
                clipped_cube, _, _, _, _ = moment_maps(self.galaxy.name, self.path_pbcorr, self.path_uncorr,
                                                       sun=self.sun, tosave=False).calc_moms()
        else:
            clipped_cube = fits.open(self.path + 'clipped_cube.fits')[0]

        vres = clipped_cube.header['CDELT3'] / 1000.  # velocity resolution
        res = clipped_cube.header['CDELT2']  # degrees per pixel
        beampix = clipped_cube.header['BMAJ'] / res  # beam size in pixels
        slitsize = np.ceil(beampix / 2)

        # Rotate the cube along the spatial axes, so that the galaxy lies horizontal
        cube_rot = ndimage.interpolation.rotate(clipped_cube.data, self.galaxy.angle, axes=(1, 2), reshape=True)

        # If you still have to determine where the slit should be, show a projection of the rotated cube, and return
        if findcentre:
            plt.imshow(np.sum(cube_rot, axis=0))
            return

        # Define a slit around the centre of the galaxy with a width of the beam size (or use the full width of the galaxy)
        if full_width:
            slit = cube_rot[:, self.galaxy.centre_pvd - self.galaxy.size/2:self.galaxy.centre_pvd + self.galaxy.size/2, :]
        else:
            slit = cube_rot[:, int(self.galaxy.centre_pvd - slitsize):int(self.galaxy.centre_pvd + slitsize), :]

        # Collapse along the slit to create the PV diagram
        PV = np.sum(slit, axis=1)

        # There is a lot of garbage because of the interpolation used by the rotation function, define a lower limit to get rid of that
        PV[PV < 0.001] = 0

        # Show the PVD to determine where the velocity centre is by eye
        if find_velcentre:
            plt.imshow(PV)
            return

        # OR DO THIS??
        # find the middle of the galaxy, where the relative velocity is closest to zero
        #mid_v = np.where((vel_array - sysvel) > 0)[0][0]  # index of velocity axis where the velocity is closest to the system velocity
        #line = PV[mid_v, :]  # line along the spatial axis where this is true
        #emission = line[line > 0.01]  # chunk of that line where there's actually emission
        #mid_emis = emission[int(len(emission) / 2.)]  # the middle of the chunck with emission
        #middle_PV = np.where(PV[mid_v, :] == mid_emis)[0]  # index of the middle of the chunck with emission

        # Create a physical position axis using the middle found and the spatial resolution, in arc seconds
        position = np.arange(0, PV.shape[1], 1)
        offset = (position - self.galaxy.centre_pvd) * res * 3600

        # Plot the PVD
        fig, ax = plt.subplots(figsize=(10, 7))
        ax_kpc = ax.twiny()  # add second x-axis at the top of the figure (position in kpc)
        ax_rel = ax.twinx()  # add second y-axis to the right of the figure (relative velocity)

        # Create the extra y axis in absolute velocities
        def absolute_yax(ax):
            """
            Convert velocities to relative velocities for the second x axis
            """

            def absolute_vel(v):
                """
                Convert the relative velocities to absolute velocities
                """
                return v + self.galaxy.sysvel

            y1, y2 = ax.get_ylim()
            ax_rel.set_ylim(absolute_vel(y1), absolute_vel(y2))
            ax_rel.figure.canvas.draw()

        # Create the extra x axis in kpc
        def x_kpc(ax):
            """
            Convert velocities to relative velocities for the second x axis.
            arcseconds to degrees to radians, and Mpc to kpc
            """

            def deg_to_kpc(arcsecs):
                return self.galaxy.distance * arcsecs / 3600. * (np.pi / 180.) * 1000.

            x1, x2 = ax.get_xlim()
            ax_kpc.set_xlim(deg_to_kpc(x1), deg_to_kpc(x2))
            ax_kpc.figure.canvas.draw()

        # Automatically update ylim of the relative velocity axis when ylim of ax changes.
        ax.callbacks.connect("ylim_changed", absolute_yax)
        ax.callbacks.connect("xlim_changed", x_kpc)

        velocity, _, _ = moment_maps(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                     tosave=False).create_vel_array(clipped_cube)

        levels = list(np.linspace(np.amin(PV), np.amax(PV), 20))
        levels = 20

        # Contours in black
        ax.contour(PV,
                   extent=[np.amin(offset), np.amax(offset), velocity[0] - self.galaxy.sysvel, velocity[-1] - self.galaxy.sysvel],
                   colors='k', levels=levels,
                   linewidths=1)

        # Filling with coloured contours
        ax.contourf(PV,
                    extent=[np.amin(offset), np.amax(offset), velocity[0] - self.galaxy.sysvel, velocity[-1] - self.galaxy.sysvel],
                    cmap='afmhot_r', vmin=0.1 * np.amax(PV), vmax=0.8 * np.amax(PV),
                    levels=levels)

        # Make the plot pretty
        #ax.set_xlim(-1.5 * self.galaxy.size * res * 3600, 1.5 * self.galaxy.size * res * 3600)
        ax.set_ylim(1.1 * -self.galaxy.vrange, 1.1 * self.galaxy.vrange)
        ax.set_xlabel('Offset [arcsec]')
        ax_kpc.set_xlabel('Offset [kpc]')
        ax.set_ylabel(r'Relative velocity [km s$^{-1}$]')
        ax_rel.set_ylabel(r'Velocity [km s$^{-1}$]')
        x1, x2 = ax.get_xlim()
        y1, y2 = ax.get_ylim()
        ax.errorbar(0.8 * x2, 0.7 * y2, xerr=clipped_cube.header['BMAJ'] * 3600 / 2., yerr=vres / 2., ecolor='k', capsize=2.5)
        ax.annotate('PA = ' + str(self.galaxy.angle) + '$^o$', xy=(-0.4 * x2, -0.7 * y2), fontsize=20)
        plt.tight_layout()

        if self.tosave:
            plt.savefig(self.savepath + 'PVD.pdf', bbox_inches='tight')

    def spectrum(self):

        spectrum = moment_maps(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                               tosave=False).spectrum()
        cube_pbcorr, cube_uncorr = moment_maps(self.galaxy.name, self.path_pbcorr, self.path_uncorr,
                                               sun=self.sun, tosave=False).readfits()
        _, _, velocity = moment_maps(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                     tosave=False).create_vel_array(cube_pbcorr)

        spectrum = spectrum[self.galaxy.start - 5:self.galaxy.stop + 5]
        velocity = velocity[self.galaxy.start - 5:self.galaxy.stop + 5]

        fig, ax = plt.subplots(figsize=(7, 7))
        ax.plot(velocity, spectrum, color='k', drawstyle='steps')

        # Line through zero
        x = np.arange(np.amin(velocity) - 100, np.amax(velocity) + 100, 1)
        zeroline = np.zeros(len(x))
        plt.plot(x, zeroline, linestyle=':', c='r', linewidth=1)

        # Set various parameters
        ax.set_xlim(velocity[len(velocity) - 1] + 5, velocity[0] - 5)
        ax.set_xlabel(r'Velocity [km s$^{-1}$]')
        #ax.set_ylabel('Flux density [mJy]')
        ax.set_ylabel('Flux density [K]')

        plt.tight_layout()

        if self.tosave:
            plt.savefig(self.savepath + 'spectrum.pdf', bbox_inches='tight')
