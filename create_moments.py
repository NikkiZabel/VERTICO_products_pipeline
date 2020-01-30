from astropy.io import fits
import numpy as np
import scipy.ndimage as ndimage
from targets import galaxies
from clip_cube import ClipCube


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
        #header.pop('OBSGEO-Z')

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
            header['CLIPLEVEL_LOW'] = self.galaxy.cliplevel_low
            header.comments['CLIPLEVEL_LOW'] = 'Lower clip level using the Sun clipping method'
            header['CLIPLEVEL_HIGH'] = self.galaxy.cliplevel_high
            header.comments['CLIPLEVEL_HIGH'] = 'Higher clip level using the Sun clipping method'
            header['NCHAN_LOW'] = self.galaxy.cliplevel_low
            header.comments['NCHAN_LOW'] = 'Lower number of consecutive channels using the Sun clipping method'
            header['NCHAN_HIGH'] = self.galaxy.cliplevel_high
            header.comments['NCHAN_HIGH'] = 'Higher number of consecutive channels using the Sun clipping method'
        else:
            header['CLIPLEVEL'] = self.galaxy.cliplevel
            header.comments['CLIPLEVEL'] = 'SNR to which the spectral cube was smooth clipped using the method of' \
                                           'Dame 2011'

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
            mom0 *= alpha_co
        elif units == 'K km/s':
            pass
        else:
            raise AttributeError('Please choose between "K km/s" and "M_Sun/pc^2"')

        mom1 = np.sum(cube.data * vel_narray, axis=0) / np.sum(cube.data, axis=0)
        mom2 = np.sqrt(np.sum(abs(cube.data) * (vel_narray - mom1) ** 2., axis=0) / np.sum(abs(cube.data), axis=0))

        # Calculate the systemic velocity from the spatial inner part of the cube (to avoid PB effects)
        inner_cube = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, savepath=self.savepath,
                              tosave=self.tosave).innersquare(mom1)

        if not self.galaxy.sysvel:
            sysvel = np.nanmean(inner_cube)
        mom1 -= sysvel

        mom0_hdu = fits.PrimaryHDU(mom0, self.new_header(cube.header))
        mom1_hdu = fits.PrimaryHDU(mom1, self.new_header(cube.header))
        mom2_hdu = fits.PrimaryHDU(mom2, self.new_header(cube.header))

        # Change or add any (additional) keywords to the headers
        if units == 'M_Sun/pc^2': mom0_hdu.header['BTYPE'] = 'Column density'
        else: mom0_hdu.header['BTYPE'] = 'Surface brightness'
        mom1_hdu.header['BTYPE'] = 'Velocity'
        mom2_hdu.header['BTYPE'] = 'Linewidth'
        mom0_hdu.header['BUNIT'] = units; mom0_hdu.header.comments['BUNIT'] = ''
        mom1_hdu.header['BUNIT'] = 'km/s'; mom1_hdu.header.comments['BUNIT'] = ''
        mom2_hdu.header['BUNIT'] = 'km/s'; mom2_hdu.header.comments['BUNIT'] = ''
        mom0_hdu.header['ALPHA_CO'] = alpha_co
        mom1_hdu.header['SYSVEL'] = sysvel; mom1_hdu.header.comments['SYSVEL'] = 'km/s'
        self.add_clipping_keywords(mom0_hdu.header)
        self.add_clipping_keywords(mom1_hdu.header)
        self.add_clipping_keywords(mom2_hdu.header)

        if self.tosave:
            cube.writeto(self.savepath + 'clipped_cube.fits', overwrite=True)
            mom0_hdu.writeto(self.savepath + 'moment0.fits', overwrite=True)
            mom1_hdu.writeto(self.savepath + 'moment1.fits', overwrite=True)
            mom2_hdu.writeto(self.savepath + 'moment2.fits', overwrite=True)

        return cube, mom0_hdu, mom1_hdu, mom2_hdu, sysvel

    def PVD(self, axis='major', findcentre=False, find_velcentre=False, full_width=False):

        clipped_cube, _, _, _, _ = self.calc_moms()

        res = clipped_cube.header['CDELT2']  # degrees per pixel
        beampix = clipped_cube.header['BMAJ'] / res  # beam size in pixels
        slitsize = np.ceil(beampix / 2)

        # Rotate the cube along the spatial axes, so that the galaxy lies horizontal
        if axis == 'major':
            rot_angle = self.galaxy.angle
        elif axis == 'minor':
            rot_angle = self.galaxy.angle + 90
        else:
            raise AttributeError('Please choose between "major" and "minor" for the "axis" keyword')

        cube_rot = ndimage.interpolation.rotate(clipped_cube.data, rot_angle, axes=(1, 2), reshape=True)

        # If you still have to determine where the slit should be, show a projection of the rotated cube, and return
        if findcentre:
            from matplotlib import pyplot as plt
            plt.imshow(np.sum(cube_rot, axis=0))
            return

        # Define a slit around the centre of the galaxy with a width of the beam size (or use the full width of the galaxy)
        if full_width:
            slit = cube_rot[:, self.galaxy.centre_pvd - self.galaxy.size / 2:self.galaxy.centre_pvd + self.galaxy.size / 2,
                   :]
        else:
            slit = cube_rot[:, int(self.galaxy.centre_pvd - slitsize):int(self.galaxy.centre_pvd + slitsize), :]

        # Collapse along the slit to create the PV diagram
        PV = np.sum(slit, axis=1)

        # There is a lot of garbage because of the interpolation used by the rotation function, define a lower limit to get rid of that
        PV[PV < 0.001] = 0

        # Show the PVD to determine where the velocity centre is by eye
        if find_velcentre:
            from matplotlib import pyplot as plt
            plt.imshow(PV)
            return

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

        # I'M NOT SURE WHAT HAPPENS TO THESE ???
        pvd_header['LONPOLE'] = clipped_cube.header['LONPOLE']
        pvd_header['LATPOLE'] = clipped_cube.header['LATPOLE']
        pvd_header['PC1_1'] = clipped_cube.header['PC1_1']
        pvd_header['PC2_1'] = clipped_cube.header['PC2_1']
        pvd_header['PC1_2'] = clipped_cube.header['PC1_2']
        pvd_header['PC2_2'] = clipped_cube.header['PC2_2']
        pvd_header['CTYPE1'] = '?????'
        pvd_header['CRVAL1'] = '?????'
        pvd_header['CDELT1'] = '?????'
        pvd_header['CRPIX'] = '?????'
        pvd_header['CUNIT1'] = '?????'
        pvd_header['CTYPE2'] = 'VRAD'
        pvd_header['CRVAL2'] = '????'
        pvd_header['CDELT2'] = '????'
        pvd_header['CRPIX2'] = '????'
        pvd_header['CUNIT2'] = 'km/s'

        pvd_header['RESTFRQ'] = clipped_cube.header['RESTFRQ']
        pvd_header.comments['RESTFRQ'] = clipped_cube.header.comments['RESTFRQ']
        pvd_header['SPECSYS'] = clipped_cube.header['SPECSYS']
        pvd_header.comments['SPECSYS'] = clipped_cube.header.comments['SPECSYS']


        ### WHAT ARE ALTRVAL AND ALTRPIX AND VELREF AND USEWEIGH AND JTOK???


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

        return pvd_hdu

    def spectrum(self):
        """
        Calculate the spectrum from the spectral cube.
        :return: array containing the spectrum in whatever units the cube is in, without the beam^-1 (so probably K or
        (m)Jy)
        """

        cube_pbcorr, cube_uncorr = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun, savepath=self.savepath,
                        tosave=self.tosave).readfits()

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
            np.savetxt(self.savepath + 'spectrum.txt', spectrum)
            np.savetxt(self.savepath + 'spectrum_velocities.txt', spectrum_velocities)
            np.savetxt(self.savepath + 'spectrum_frequencies.txt', spectrum_frequencies)

        # Estimate the rms in the spectrum
        #emis, noise = self.splitCube(cutout, self.galaxy.start, self.galaxy.stop)
        #rms = np.std(np.sum(noise, axis=(1, 2)))
        #np.savetxt(path + 'specnoise.txt', [rms / beamsize])

        return spectrum, spectrum_velocities, spectrum_frequencies # / beamsize