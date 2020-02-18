from astropy.io import fits
import numpy as np
import scipy.ndimage as ndimage
from targets import galaxies
from clip_cube import ClipCube
from photutils import EllipticalAnnulus
from photutils import aperture_photometry


class MomentMaps:

    def __init__(self, galname, path_pbcorr, path_uncorr, savepath=None, sun=True, tosave=False):
        self.galaxy = galaxies(galname)
        self.path_pbcorr = path_pbcorr
        self.path_uncorr = path_uncorr
        self.savepath = savepath or './'
        self.sun = sun
        self.tosave = tosave

    def makebeam(self, xpixels, ypixels, header, rot=0, cent=0):
        """
        Creates the psf from the beam size given a custom location in the image and rotation, assuming a 2D Gaussian.
        :param xpixels (int or float): number of pixels in the x-axis
        :param ypixels (int or float): number of pixels in the y-axis
        :param header (fits-header): fits header of the corresponding spectral cube
        :param rot (int or float, optional): degrees over which the psf should be rotated
        :param cent (list or ndarray of length 2, optional): custom centre of the psf
        :return (2D array): psf of the beam
        """

        # Extract relevant information from header
        res = header['CDELT2']  # degrees per pixel
        bmaj = header['BMAJ']  # degrees
        bmin = header['BMIN']  # degrees
        beam = np.array([bmaj, bmin]) / res  # Beam FWHM in pixels

        # Convert the FWHM of the beam to the std dev in the Gaussian distribution
        sigma = beam / np.sqrt(8. * np.log(2.))

        # If centre is not defined, make it the centre of the image
        if not cent: cent = [xpixels / 2., ypixels / 2]

        # Implement a rotation if specified
        if np.tan(np.radians(rot)) == 0:
            dirfac = 1
        else:
            dirfac = np.sign(np.tan(np.radians(rot)))

        # Calculate the psf
        x, y = np.indices((int(xpixels), int(ypixels)), dtype='float')
        x -= cent[0]
        y -= cent[1]

        a = (np.cos(np.radians(rot)) ** 2) / (2.0 * (sigma[0] ** 2)) + (np.sin(np.radians(rot)) ** 2) / (
        2.0 * (sigma[1] ** 2))
        b = ((dirfac) * (np.sin(2.0 * np.radians(rot)) ** 2) / (4.0 * (sigma[0] ** 2))) + (
        (-1 * dirfac) * (np.sin(2.0 * np.radians(rot)) ** 2) / (4.0 * (sigma[1] ** 2)))
        c = (np.sin(np.radians(rot)) ** 2) / (2.0 * (sigma[0] ** 2)) + (np.cos(np.radians(rot)) ** 2) / (
        2.0 * (sigma[1] ** 2))

        psf = np.exp(-1 * (a * (x ** 2) - 2.0 * b * (x * y) + c * (y ** 2)))

        return psf

    def new_header(self, header):
        """
        Remove the velocity axis from a HDU header, so it corresponds to the 2D version of the corresponding data cube.
        :param header (HDU header): header of the original data cube
        :return: input header, but with velocity axis related keywords removed.
        """

        header = header.copy()

        try:
            header.pop('PC3_1')
            header.pop('PC3_2')
            header.pop('PC1_3')
            header.pop('PC2_3')
            header.pop('PC3_3')
        except:
            header.pop('PC03_01')
            header.pop('PC03_03')
            header.pop('PC03_02')
            header.pop('PC01_03')
            header.pop('PC02_03')

        header.pop('CTYPE3')
        header.pop('CRVAL3')
        header.pop('CDELT3')
        header.pop('CRPIX3')
        header.pop('CUNIT3')
        header.pop('NAXIS3')
        header.pop('OBSGEO-Z')

        header['NAXIS'] = 2

        return header

    def create_vel_array(self, cube):
        """
        From the relevant header keywords, create an array with the velocities in km/s corresponding to the spectral
        axis of the spectral cube
        :param cube (HDU file): HDU file of the spectral cube for which we want to make the velocity array
        :return: three arrays: the velocity array in one dimension, the length equals the numbers of channels that
        contain emission, the same velocity array but in the shape of the spectral cube, and the one-dimensional
        velocity array corresponding to the entire velocity axis of the cube (including line-free channels)
        """
        v_val = cube.header['CRVAL3'] / 1000.  # Velocity in the reference channel, m/s to km/s
        v_step = cube.header['CDELT3'] / 1000.  # Velocity step in each channel, m/s to km/s
        v_ref = cube.header['CRPIX3']  # Location of the reference channel

        # Construct the velocity arrays (keep in mind that fits-files are 1 indexed)
        vel_array = (np.arange(0, len(cube.data[:, 0, 0])) - v_ref + 1 + self.galaxy.start) * v_step + v_val
        vel_narray = np.tile(vel_array, (len(cube.data[0, 0, :]), len(cube.data[0, :, 0]), 1)).transpose()
        vel_array_full = (np.arange(0, len(cube.data[:, 0, 0])) - v_ref + 1) * v_step + v_val

        return vel_array, vel_narray, vel_array_full

    def add_clipping_keywords(self, header):
        if self.sun:
            header['CLIPL_L'] = self.galaxy.cliplevel_low
            header.comments['CLIPL_L'] = 'Lower clip level (Sun method)'
            header['CLIPL_H'] = self.galaxy.cliplevel_high
            header.comments['CLIPL_H'] = 'Higher clip level (Sun method)'
            header['NCHAN_L'] = self.galaxy.cliplevel_low
            header.comments['NCHAN_L'] = 'Lower number of consec. chans (Sun method)'
            header['NCHAN_H'] = self.galaxy.cliplevel_high
            header.comments['NCHAN_H'] = 'Higher number of consec. chans (Sun method)'
        else:
            header['CLIPL'] = self.galaxy.cliplevel
            header.comments['CLIPL'] = 'SNR for smooth clip (Dame11)' \

        return header

    def calc_moms(self, units='M_Sun/pc^2', alpha_co=6.25):
        """
        Clip the spectral cube according to the desired method, and create moment 0, 1, and 2 maps. Save them as fits
        files if so desired. Also calculate the systemic velocity from the moment 1 map.
        :param units (string): desired units for the moment 0 map. Default is M_Sun/pc^2, the alternative is K km/s.
        :param alpha_co (float): in case units == 'M_Sun/pc^2', multiply by alpha_co to obtain these units. Default
        value for CO(2-1) from https://arxiv.org/pdf/1805.00937.pdf.
        :return: clipped spectral cube, HDUs of the moment 0, 1, and 2 maps, and the systemic velocity in km/s
        """

        cube = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun, savepath=self.savepath,
                        tosave=self.tosave).do_clip()

        vel_array, vel_narray, vel_fullarray = self.create_vel_array(cube)

        mom0 = np.sum((cube.data * abs(cube.header['CDELT3']) / 1000), axis=0)
        if units == 'M_Sun/pc^2':
            mom0 = mom0 / cube.header['JTOK'] * 91.9 * alpha_co * (cube.header['BMAJ'] * 3600 * cube.header[
                'BMIN'] * 3600) ** (-1)
        elif units == 'K km/s':
            pass
        else:
            raise AttributeError('Please choose between "K km/s" and "M_Sun/pc^2"')

        mom1 = np.sum(cube.data * vel_narray, axis=0) / np.sum(cube.data, axis=0)
        mom2 = np.sqrt(np.sum(abs(cube.data) * (vel_narray - mom1) ** 2., axis=0) / np.sum(abs(cube.data), axis=0))

        # Calculate the systemic velocity from the spatial inner part of the cube (to avoid PB effects)
        inner_cube = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, savepath=self.savepath,
                              tosave=self.tosave).innersquare(mom1)

        sysvel = np.nanmean(inner_cube) + self.galaxy.sysvel_offset
        mom1 -= sysvel

        mom0_hdu = fits.PrimaryHDU(mom0, self.new_header(cube.header))
        mom1_hdu = fits.PrimaryHDU(mom1, self.new_header(cube.header))
        mom2_hdu = fits.PrimaryHDU(mom2, self.new_header(cube.header))

        # Change or add any (additional) keywords to the headers
        if units == 'M_Sun/pc^2': mom0_hdu.header['BTYPE'] = 'Column density'
        else: mom0_hdu.header['BTYPE'] = 'Integrated intensity'
        mom1_hdu.header['BTYPE'] = 'Velocity'
        mom2_hdu.header['BTYPE'] = 'Linewidth'
        mom0_hdu.header['BUNIT'] = units; mom0_hdu.header.comments['BUNIT'] = ''
        mom1_hdu.header['BUNIT'] = 'km/s'; mom1_hdu.header.comments['BUNIT'] = ''
        mom2_hdu.header['BUNIT'] = 'km/s'; mom2_hdu.header.comments['BUNIT'] = ''
        mom0_hdu.header['ALPHA_CO'] = alpha_co; mom0_hdu.header.comments['ALPHA_CO'] = 'Assuming a line ratio of 0.7'
        mom1_hdu.header['SYSVEL'] = sysvel; mom1_hdu.header.comments['SYSVEL'] = 'km/s'
        self.add_clipping_keywords(mom0_hdu.header)
        self.add_clipping_keywords(mom1_hdu.header)
        self.add_clipping_keywords(mom2_hdu.header)

        if self.tosave:
            cube.writeto(self.savepath + 'clipped_cube.fits', overwrite=True)
            if units == 'M_Sun/pc^2':
                mom0_hdu.writeto(self.savepath + 'moment0_M_Sun.fits', overwrite=True)
            if units == 'K km/s':
                mom0_hdu.writeto(self.savepath + 'moment0_K.fits', overwrite=True)
            mom1_hdu.writeto(self.savepath + 'moment1.fits', overwrite=True)
            mom2_hdu.writeto(self.savepath + 'moment2.fits', overwrite=True)

        return cube, mom0_hdu, mom1_hdu, mom2_hdu, sysvel

    def PVD(self, axis='major', full_width=False):

        clipped_cube, _, _, _, sysvel = self.calc_moms()

        res = clipped_cube.header['CDELT2']  # degrees per pixel
        beampix = clipped_cube.header['BMAJ'] / res  # beam size in pixels
        slitsize = np.ceil(beampix / 2)
        shift_x = self.galaxy.centre_x - clipped_cube.shape[1] / 2
        shift_y = self.galaxy.centre_y - clipped_cube.shape[2] / 2

        # Rotate the cube along the spatial axes, so that the galaxy lies horizontal
        if axis == 'major':
            rot_angle = self.galaxy.angle
        elif axis == 'minor':
            rot_angle = self.galaxy.angle + 90
        else:
            raise AttributeError('Please choose between "major" and "minor" for the "axis" keyword')

        if shift_x > 0:
            temp = np.zeros((clipped_cube.shape[0], clipped_cube.shape[1] + int(abs(shift_x)), clipped_cube.shape[2]))
            temp[:, :-int(abs(shift_x)), :] = clipped_cube.data
        elif shift_x < 0:
            temp = np.zeros((clipped_cube.shape[0], clipped_cube.shape[1] + int(abs(shift_x)), clipped_cube.shape[2]))
            temp[:, int(abs(shift_x)):, :] = clipped_cube.data
        else:
            temp = clipped_cube.data

        if shift_y > 0:
            pvdcube = np.zeros((temp.shape[0], temp.shape[1], temp.shape[2] + int(abs(shift_y))))
            pvdcube[:, :, :-int(abs(shift_y))] = temp
        elif shift_y < 0:
            pvdcube = np.zeros((temp.shape[0], temp.shape[1], temp.shape[2] + int(abs(shift_y))))
            pvdcube[:, :, int(abs(shift_y)):] = temp
        else:
            pvdcube = temp

        cube_rot = ndimage.interpolation.rotate(pvdcube, rot_angle, axes=(1, 2), reshape=True)

        # Define a slit around the centre of the galaxy with a width of the beam size (or use the full width of the galaxy)
        if full_width:
            slit = cube_rot[:, cube_rot.shape[1] / 2 - self.galaxy.size / 2:cube_rot.shape[1] / 2 +
                                                                                self.galaxy.size / 2, :]
        else:
            slit = cube_rot[:, int(cube_rot.shape[1] / 2 - slitsize):int(cube_rot.shape[1] / 2 + slitsize), :]

        # Collapse along the slit to create the PV diagram
        PV = np.sum(slit, axis=1)

        # There is a lot of garbage because of the interpolation used by the rotation function, define a lower limit to get rid of that
        PV[PV < 0.001] = 0

        # Create an appropriate header
        pvd_header = fits.Header()
        pvd_header['SIMPLE'] = clipped_cube.header['SIMPLE']
        pvd_header.comments['SIMPLE'] = clipped_cube.header.comments['SIMPLE']
        pvd_header['BITPIX'] = clipped_cube.header['BITPIX']
        pvd_header.comments['BITPIX'] = clipped_cube.header.comments['BITPIX']
        pvd_header['NAXIS'] = 2
        pvd_header.comments['NAXIS'] = clipped_cube.header.comments['NAXIS']
        pvd_header['NAXIS1'] = PV.shape[0]
        pvd_header['NAXIS2'] = PV.shape[1]
        pvd_header['PVD_AXIS'] = axis
        pvd_header['PA'] = -(self.galaxy.angle - 360 - 90)
        pvd_header['BMAJ'] = clipped_cube.header['BMAJ']
        pvd_header['BMIN'] = clipped_cube.header['BMIN']
        pvd_header['BPA'] = clipped_cube.header['BPA']
        pvd_header['OBJECT'] = clipped_cube.header['OBJECT']
        pvd_header['EQUINOX'] = 2000
        pvd_header['RADESYS'] = 'FK5'
        pvd_header['LONPOLE'] = clipped_cube.header['LONPOLE']
        pvd_header['LATPOLE'] = clipped_cube.header['LATPOLE']
        pvd_header['CTYPE1'] = 'OFFSET'
        pvd_header['CRVAL1'] = 0
        pvd_header['CDELT1'] = clipped_cube.header['CDELT1'] / (cube_rot.shape[2] / pvdcube.shape[2])
        pvd_header['CRPIX1'] = np.ceil(PV.shape[1] / 2)
        pvd_header['CUNIT1'] = clipped_cube.header['CUNIT1']
        pvd_header['CTYPE2'] = 'VRAD'
        pvd_header['CRVAL2'] = clipped_cube.header['CRVAL3']
        pvd_header['CDELT2'] = clipped_cube.header['CDELT3']
        pvd_header['CRPIX2'] = clipped_cube.header['CRPIX3']
        pvd_header['CUNIT2'] = 'km/s'
        pvd_header['PC1_1'] = pvd_header['CDELT1']
        pvd_header['PC2_1'] = clipped_cube.header['PC2_1']
        pvd_header['PC1_2'] = clipped_cube.header['PC1_2']
        pvd_header['PC2_2'] = pvd_header['CDELT2']
        pvd_header['SYSVEL'] = sysvel + self.galaxy.sysvel_offset
        pvd_header['RESTFRQ'] = clipped_cube.header['RESTFRQ']
        pvd_header.comments['RESTFRQ'] = clipped_cube.header.comments['RESTFRQ']
        pvd_header['SPECSYS'] = clipped_cube.header['SPECSYS']
        pvd_header.comments['SPECSYS'] = clipped_cube.header.comments['SPECSYS']
        pvd_header['ALTRVAL'] = clipped_cube.header['ALTRVAL']
        pvd_header['ALTRPIX'] = clipped_cube.header['ALTRPIX']
        pvd_header['VELREF'] = clipped_cube.header['VELREF']
        pvd_header['USEWEIGH'] = clipped_cube.header['USEWEIGH']
        pvd_header['JTOK'] = clipped_cube.header['JTOK']
        pvd_header['OBSRA'] = clipped_cube.header['OBSRA']
        pvd_header['OBSDEC'] = clipped_cube.header['OBSDEC']
        pvd_header['OBSGEO-X'] = clipped_cube.header['OBSGEO-X']
        pvd_header['OBSGEO-Y'] = clipped_cube.header['OBSGEO-Y']
        pvd_header['OBSGEO-Z'] = clipped_cube.header['OBSGEO-Z']
        pvd_header['DISTANCE'] = self.galaxy.distance
        pvd_header.comments['DISTANCE'] = 'Mpc'
        pvd_header['ORIGIN'] = clipped_cube.header['ORIGIN']

        pvd_hdu = fits.PrimaryHDU(PV, pvd_header)

        if self.tosave:
            pvd_hdu.writeto(self.savepath + 'PVD.fits', overwrite=True)

        return pvd_hdu, shift_x

    def spectrum(self):
        """
        Calculate the spectrum from the spectral cube.
        :return: array containing the spectrum in whatever units the cube is in, without the beam^-1 (so probably K or
        (m)Jy)
        """

        cube_pbcorr, cube_uncorr = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun, savepath=self.savepath,
                        tosave=self.tosave).readfits()

        #clipped_cube, _, _, _, sysvel = self.calc_moms()
        #cube_pbcorr = clipped_cube

        # Calculate the beam size, so we can divide by this to get rid of the beam^-1 in the units.
        #psf = self.makebeam(cube_pbcorr.shape[1], cube_pbcorr.shape[2], cube_pbcorr.header)
        #beamsize = np.sum(psf)

        # Make a cutout around the emission in the spatial direction, to reduce noise
        cutout = cube_pbcorr.data[:, self.galaxy.centre_y - self.galaxy.size:self.galaxy.centre_y + self.galaxy.size,
                 self.galaxy.centre_x - self.galaxy.size:self.galaxy.centre_x + self.galaxy.size]

        ### THIS OVERWRITES THE CUTOUT TO USE THE ENTIRE CUBE ###
        cutout = cube_pbcorr.data                       # <----------------------------------------------------------
        ### REMOVE IF WE DO WANT TO USE THE CUTOUT ###

        # Make this work if necessary
        #if custom_region:
            #region = pyregion.open(path + 'ds9.reg')
            #mask = region.get_mask(hdu=mom0_hdu)
            #mask_3d = np.tile(mask, (len(cube[:, 0, 0]), 1, 1))
            #cutout = np.where(mask_3d, cube, 0)

        spectrum = np.nansum(cutout, axis=(1, 2))
        _, _, vel_array_full = self.create_vel_array(cube_pbcorr)

        spectrum_velocities = vel_array_full[self.galaxy.start - 5:self.galaxy.stop + 5]
        spectrum = spectrum[self.galaxy.start - 5:self.galaxy.stop + 5]
        spectrum_frequencies = cube_pbcorr.header['RESTFRQ'] * (1 - spectrum_velocities / 299792.458) / 1e9

        if self.tosave:
            np.savetxt(self.savepath + 'spectrum.csv',
                       np.column_stack((spectrum, spectrum_velocities, spectrum_frequencies)), delimiter=',',
                       header='Spectrum (K), Velocity (km/s), Frequency (GHz)')

        # Estimate the rms in the spectrum
        #emis, noise = self.splitCube(cutout, self.galaxy.start, self.galaxy.stop)
        #rms = np.std(np.sum(noise, axis=(1, 2)))
        #np.savetxt(path + 'specnoise.txt', [rms / beamsize])

        psf = self.makebeam(cube_pbcorr.shape[1], cube_pbcorr.shape[2], cube_pbcorr.header)
        beamsize = np.sum(psf)
        spectrum /= beamsize
        print(np.log10(3.93e-17 * 16.5 ** 2. * 2e20 / 0.7 * np.trapz(np.flip(spectrum), np.flip(spectrum_velocities)) / cube_pbcorr.header['JTOK']))

        return spectrum, spectrum_velocities, spectrum_frequencies # / beamsize

    def radial_profile(self, alpha_co=6.25, table_path=None, check_aperture=False):

        _, mom0_hdu, _, _, _ = self.calc_moms(units='M_Sun/pc^2', alpha_co=alpha_co)
        beam_pix = mom0_hdu.header['BMAJ'] / mom0_hdu.header['CDELT2']

        # Estimate the rms from the spatial inner part of the cube
        cube_pbcorr, cube_uncorr = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                            savepath=self.savepath,
                                            tosave=self.tosave).readfits()
        emiscube, noisecube = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                            savepath=self.savepath,
                                            tosave=self.tosave).split_cube(cube_pbcorr)
        inner = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                            savepath=self.savepath,
                                            tosave=self.tosave).innersquare(noisecube.data)
        rms = np.nanstd(inner)
        rms = rms / cube_pbcorr.header['JTOK']
        rms = rms * abs(cube_pbcorr.header['CDELT3']) / 1000 * 91.9 * alpha_co * (cube_pbcorr.header['BMAJ'] * 3600 *
                cube_pbcorr.header['BMIN'] * 3600) ** (-1)

        if self.galaxy.eccentricity:
            e = self.galaxy.eccentricity
        elif self.galaxy.inclination:
            e = np.sin(np.deg2rad(self.galaxy.inclination))
        elif not table_path:
            raise AttributeError('Please provide the inclination of the galaxy or its projected eccentricity,'
                                 'or provide the path to the VERTICO master table to read it from there.')
        else:
            print('Reading the inclination from the master table.')
            table = fits.open(table_path)
            inc = table[1].data['inclination'][table[1].data['Galaxy'] == self.galaxy.name]
            e = np.sin(np.deg2rad(inc))[0]

        centre = (self.galaxy.centre_y, self.galaxy.centre_x)
        b_in = -beam_pix + 0.000000000001
        b_out = 0
        theta = self.galaxy.angle - 180

        rad_prof = []
        area = []
        emission = 2112
        area_temp = 1

        if check_aperture:
            from matplotlib import pyplot as plt
            plt.figure()
            plt.imshow(mom0_hdu.data)

        while emission / area_temp > 1:
            b_in += beam_pix
            b_out += beam_pix
            if emission == 2112:
                a_in = 0.0000000000001
            else:
                a_in = a_out

            if e == 1:
                a_out = b_out
            else:
                a_out = b_out / np.sqrt(1 - e ** 2)

            aperture = EllipticalAnnulus(centre, a_in, a_out, b_out, theta)

            if check_aperture:
                aperture.plot(color='red')

            emission = aperture_photometry(mom0_hdu.data, aperture)['aperture_sum'][0] #* \
                      # (np.deg2rad(mom0_hdu.header['CDELT2']) * self.galaxy.distance * 1e6) ** 2

            #psf = self.makebeam(mom0_hdu.shape[0], mom0_hdu.shape[1], mom0_hdu.header)
            #beamsize = np.sum(psf)
            #emission = aperture_photometry(mom0_hdu.data, aperture)['aperture_sum'][0] * 3.93e-17 * 16.5 ** 2. * 2e20 \
            #/ 0.7 / mom0_hdu.header['JTOK'] / beamsize

            area_temp = aperture.area
            area.append(area_temp)
            rad_prof.append(emission / area_temp)

        rad_prof = rad_prof[:-1]
        area = area[:-1]
        radii_deg = (np.arange(len(rad_prof)) + 1) * mom0_hdu.header['CDELT2']
        radii_kpc = np.deg2rad(radii_deg) * self.galaxy.distance * 1000

        beam_area_pc = np.pi * np.deg2rad(mom0_hdu.header['BMAJ']) * np.deg2rad(mom0_hdu.header['BMIN']) * \
                       (self.galaxy.distance * 1e6) ** 2

        if self.tosave:
            np.savetxt(self.savepath + 'radial_profile.csv',
                       np.column_stack((rad_prof, np.ones(len(rad_prof)) * rms, radii_deg * 3600, radii_kpc)), delimiter=',',
                       header='Surface density (M_Sun pc^-2), RMS error (M_Sun pc^-2), Radii (arcsec), Radii (kpc)')

        #print(np.log10(np.amax(rad_prof_cum)))

        return rad_prof, radii_deg * 3600, radii_kpc, rms
