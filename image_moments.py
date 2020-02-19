from matplotlib import pyplot as plt
import aplpy as apl
from astropy import wcs
from astropy.nddata.utils import Cutout2D
from astropy.io import fits
import numpy as np
import config_figs
from sauron_colormap import register_sauron_colormap; register_sauron_colormap()
from targets import galaxies
from create_moments import MomentMaps
from clip_cube import ClipCube


class CreateImages:

    def __init__(self, galname, path_pbcorr, path_uncorr, savepath=None, refresh=False, overwrite=True,
                 make_cutout=False, sun=True, tosave=False):
        self.galaxy = galaxies(galname)
        self.path_pbcorr = path_pbcorr
        self.path_uncorr = path_uncorr
        self.refresh = refresh
        self.overwrite = overwrite
        self.savepath = savepath + self.galaxy.name + '_' or './' + self.galaxy.name + '_'
        self.make_cutout = make_cutout
        self.sun = sun
        self.tosave = tosave

    def moment_zero(self, units='M_Sun/pc^2', path=''):

        if self.refresh:
            if units == 'M_Sun/pc^2':
                if self.overwrite:
                    _, image, _, _, _ = MomentMaps(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                                    savepath=self.savepath, tosave=True).calc_moms(units='M_Sun/pc^2')
                else:
                    _, image, _, _, _ = MomentMaps(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                                    tosave=False).calc_moms(units='M_Sun/pc^2')
            elif units == 'K km/s':
                if self.overwrite:
                    _, image, _, _, _ = MomentMaps(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                                    savepath=self.savepath, tosave=True).calc_moms(units='K km/s')
                else:
                    _, image, _, _, _ = MomentMaps(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                                    tosave=False).calc_moms(units='K km/s')
            else:
                raise AttributeError('Please choose between "K km/s" and "M_Sun/pc^2"')
        else:
            image = fits.open(path + 'moment0.fits')[0]

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

        fig.show_contour(image, cmap='magma_r', levels=np.linspace(np.amax(image.data)*1e-9, np.amax(image.data), 20),
                         filled=True, overlap=True)

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
        elif np.amax(image.data) < 200:
            ticks = np.arange(0, np.amax(image.data) + 10, 20)
        elif np.amax(image.data) < 1000:
            ticks = np.arange(0, np.amax(image.data) + 20, 40)
        else:
            ticks = np.arange(0, np.amax(image.data) + 100, 200)

        cbar = f.colorbar(colors, ticks=ticks)
        #cbar.set_label(r'Surface brightness [Jy beam$^{-1}$ km s$^{-1}$]')
        if units == 'K km/s':
            cbar.set_label(r'Integrated intensity [K km s$^{-1}$]')
        elif units == 'M_Sun/pc^2':
            cbar.set_label(r'Surface density [M$_\odot$ pc$^{-2}$]')
        else:
            raise AttributeError('Please choose between "K km/s" and "M_Sun/pc^2"')

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
            if units == 'K km/s':
                plt.savefig(self.savepath + 'moment0_K.pdf', bbox_inches='tight')
            if units == 'M_Sun/pc^2':
                plt.savefig(self.savepath + 'moment0_M_Sun.pdf', bbox_inches='tight')
        return

    def moment_1_2(self, moment=1):

        if moment == 1:
            if self.refresh:
                if self.overwrite:
                    _, _, image, _, sysvel = MomentMaps(self.galaxy.name, self.path_pbcorr, self.path_uncorr,
                                                         savepath=self.savepath, sun=self.sun, tosave=True).calc_moms()
                else:
                    _, _, image, _, sysvel = MomentMaps(self.galaxy.name, self.path_pbcorr, self.path_uncorr,
                                                         sun=self.sun, tosave=False).calc_moms()
            else:
                image = fits.open(self.path + 'moment1.fits')[0]

        elif moment == 2:
            if self.refresh:
                if self.overwrite:
                    _, _, _, image, sysvel = MomentMaps(self.galaxy.name, self.path_pbcorr, self.path_uncorr,
                                                         savepath=self.savepath, sun=self.sun, tosave=True).calc_moms()
                else:
                    _, _, _, image, sysvel = MomentMaps(self.galaxy.name, self.path_pbcorr, self.path_uncorr,
                                                         sun=self.sun, tosave=False).calc_moms()
            else:
                image = fits.open(self.path + 'moment2.fits')[0]

        cube_pbcorr, cube_uncorr = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr,
                                               sun=self.sun, tosave=False).readfits()
        vel_array, _, _ = MomentMaps(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
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
            cbar.set_label(r'Linewidth [km s$^{-1}$]')

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

    def PVD(self, axis='major', findcentre=False, find_velcentre=False, full_width=False):
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
                PV, shift_x = MomentMaps(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                 savepath=self.savepath, tosave=True).\
                    PVD(axis=axis, full_width=full_width)
            else:
                PV, shift_x = MomentMaps(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun, tosave=False).\
                    PVD(axis=axis, full_width=full_width)
        else:
            _, shift_x = MomentMaps(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun, tosave=False). \
                PVD(axis=axis, full_width=full_width)
            PV = fits.open(self.path + 'PVD.fits')[0]

        clipped_cube, _, _, _, sysvel = MomentMaps(self.galaxy.name, self.path_pbcorr, self.path_uncorr,
                                               sun=self.sun, tosave=False).calc_moms()

        # OR DO THIS??
        # find the middle of the galaxy, where the relative velocity is closest to zero
        #mid_v = np.where((vel_array - sysvel) > 0)[0][0]  # index of velocity axis where the velocity is closest to the system velocity
        #line = PV[mid_v, :]  # line along the spatial axis where this is true
        #emission = line[line > 0.01]  # chunk of that line where there's actually emission
        #mid_emis = emission[int(len(emission) / 2.)]  # the middle of the chunck with emission
        #middle_PV = np.where(PV[mid_v, :] == mid_emis)[0]  # index of the middle of the chunck with emission

        # Create a physical position axis using the middle found and the spatial resolution, in arc seconds
        res = clipped_cube.header['CDELT2']
        vres = clipped_cube.header['CDELT3'] / 1000.  # velocity resolution
        position = np.arange(0, PV.shape[1], 1)
        offset = (position - len(position) / 2 + shift_x) * res * 3600
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
                return v + sysvel

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

        velocity, _, _ = MomentMaps(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                     tosave=False).create_vel_array(clipped_cube)

        levels = list(np.linspace(np.amin(PV.data), np.amax(PV.data), 20))
        levels = 20

        # Contours in black
        ax.contour(PV.data,
                   extent=[np.amin(offset), np.amax(offset), velocity[0] - sysvel, velocity[-1] - sysvel],
                   colors='k', levels=levels,
                   linewidths=1)

        # Filling with coloured contours
        ax.contourf(PV.data,
                    extent=[np.amin(offset), np.amax(offset), velocity[0] - sysvel, velocity[-1] - sysvel],
                    cmap='afmhot_r', vmin=0.1 * np.amax(PV.data), vmax=0.8 * np.amax(PV.data),
                    levels=levels)

        # Make the plot pretty
        #ax.set_xlim(-1.5 * self.galaxy.size * res * 3600, 1.5 * self.galaxy.size * res * 3600)
        ax.set_ylim(velocity[-1] - sysvel + 20, velocity[0] - sysvel - 20)
        ax.set_xlabel('Offset [arcsec]')
        ax_kpc.set_xlabel('Offset [kpc]')
        ax.set_ylabel(r'Relative velocity [km s$^{-1}$]')
        ax_rel.set_ylabel(r'Velocity [km s$^{-1}$]')
        x1, x2 = ax.get_xlim()
        y1, y2 = ax.get_ylim()
        ax.errorbar(0.8 * x2, 0.7 * y2, xerr=clipped_cube.header['BMAJ'] * 3600 / 2., yerr=vres / 2., ecolor='k', capsize=2.5)
        ax.annotate('PA = ' + str((self.galaxy.angle)) + '$^o$', xy=(-0.8 * x2, -0.7 * y2), fontsize=20)

        plt.tight_layout()

        if self.tosave:
            if axis == 'major':
                plt.savefig(self.savepath + 'PVD_major.pdf', bbox_inches='tight')
            if axis == 'minor':
                plt.savefig(self.savepath + 'PVD_minor.pdf', bbox_inches='tight')

    def spectrum(self, x_axis='velocity'):

        if self.refresh:
            if self.overwrite:
                spectrum, velocity, frequency = MomentMaps(self.galaxy.name, self.path_pbcorr, self.path_uncorr,
                                                    sun=self.sun, savepath=self.savepath, tosave=True).spectrum()
            else:
                spectrum, velocity, frequency = MomentMaps(self.galaxy.name, self.path_pbcorr, self.path_uncorr,
                                                    sun=self.sun, savepath=self.savepath, tosave=False).spectrum()

        else:
            temp = np.loadtxt(self.savepath + 'spectrum.csv', delimiter=',')
            spectrum = temp[:, 0]
            velocity = temp[:, 1]
            frequency = temp[:, 2]

        fig, ax = plt.subplots(figsize=(7, 7))

        if x_axis == 'velocity':
            ax.plot(velocity, spectrum, color='k', drawstyle='steps')
            x = np.arange(np.amin(velocity) - 100, np.amax(velocity) + 100, 1)
            ax.set_xlim(velocity[len(velocity) - 1] + 5, velocity[0] - 5)
            ax.set_xlabel(r'Velocity [km s$^{-1}$]')
        elif x_axis == 'frequency':
            ax.plot(frequency, spectrum, color='k', drawstyle='steps')
            x = np.arange(np.amin(frequency) - 5, np.amax(frequency) + 5, 1)
            ax.set_xlim(frequency[len(frequency) - 1], frequency[0])
            ax.set_xlabel(r'Frequency [GHz]')
        else:
            raise AttributeError('Please choose between "velocity" and "frequency" for "x-axis"')

        # Line through zero
        zeroline = np.zeros(len(x))
        plt.plot(x, zeroline, linestyle=':', c='r', linewidth=1)

        ax.set_ylabel('Brightness temperature [K]')

        plt.tight_layout()

        if self.tosave:
            if x_axis == 'frequency':
                plt.savefig(self.savepath + 'spectrum_frequency.pdf', bbox_inches='tight')
            if x_axis == 'velocity':
                plt.savefig(self.savepath + 'spectrum_velocity.pdf', bbox_inches='tight')

    def radial_profile(self, units='kpc', alpha_co=6.25, table_path=None, check_aperture=False):

        if self.refresh:
            if self.overwrite:
                rad_prof, radii_arcsec, radii_kpc, error = MomentMaps(self.galaxy.name, self.path_pbcorr,
                        self.path_uncorr, sun=self.sun, savepath=self.savepath, tosave=True).\
                    radial_profile(alpha_co=alpha_co, table_path=table_path, check_aperture=check_aperture)
            else:
                rad_prof, radii_arcsec, radii_kpc, error = MomentMaps(self.galaxy.name, self.path_pbcorr,
                        self.path_uncorr, sun=self.sun, savepath=self.savepath, tosave=False).\
                    radial_profile(alpha_co=alpha_co, table_path=table_path, check_aperture=check_aperture)

        else:
            temp = np.loadtxt(self.savepath + 'radial_profile.csv', delimiter=',')
            rad_prof = temp[:, 0]
            radii_arcsec = temp[:, 1]
            radii_kpc = temp[:, 2]

        plt.figure(figsize=(10, 7))

        if units == 'kpc':
            plt.errorbar(radii_kpc, np.log10(rad_prof), yerr=error/rad_prof * 0.434, c='k', linestyle='None', marker='o')
            plt.xlabel('Radius [kpc]')
        elif units == 'arcsec':
            plt.errorbar(radii_arcsec, np.log10(rad_prof), yerr=error/rad_prof * 0.434, c='k', linestyle='None', marker='o')
            plt.xlabel(r'Radius [$^{\prime\prime}$]')
        else:
            raise AttributeError('Please choose between "kpc" and "arcsec" for the keyword "units".')

        plt.ylabel(r'log(Surface density [$M_\odot$ pc$^{-2}$])')
        plt.xlim(-0.01)

        plt.tight_layout()

        if self.tosave:
            if units == 'kpc':
                plt.savefig(self.savepath + 'radial_profile_kpc.pdf', bbox_inches='tight')
            if units == 'arcsec':
                plt.savefig(self.savepath + 'radial_profile_arcsec.pdf', bbox_inches='tight')