from astropy.io import fits
import numpy as np
import scipy.ndimage as ndimage
from targets import galaxies
from scipy.ndimage import binary_dilation, label


class moment_maps:

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

    def remove_stokes(self, cube):
        """
        If the fits file containing the spectral cube has a Stokes axis, remove it and remove the corresponding
        header keywords.
        :param cube (HDU file): HDU file containing the spectral cube and its header
        :return (HDU file): the input HDU file with the Stokes axis and corresponding header keywords removed.
        """

        cube.data = np.squeeze(cube.data)

        try:
            cube.header.pop('PC01_04')
            cube.header.pop('PC02_04')
            cube.header.pop('PC03_04')
            cube.header.pop('PC04_04')
            cube.header.pop('PC04_01')
            cube.header.pop('PC04_02')
            cube.header.pop('PC04_03')
        except:
            cube.header.pop('PC1_4')
            cube.header.pop('PC2_4')
            cube.header.pop('PC3_4')
            cube.header.pop('PC4_4')
            cube.header.pop('PC4_1')
            cube.header.pop('PC4_2')
            cube.header.pop('PC4_3')

        cube.header.pop('CTYPE4')
        cube.header.pop('CRVAL4')
        cube.header.pop('CRPIX4')
        cube.header.pop('CUNIT4')
        cube.header.pop('CDELT4')

        return cube

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

    def readfits(self):
        """
        Read in the fits files containing the primary beam corrected and uncorrected specral cubes.
        :return: Two HDU files containing the spectral cube with and without primary beam correction applied,
        respectively, and their headers.
        """

        cube_pbcorr = fits.open(self.path_pbcorr)[0]
        cube_uncorr = fits.open(self.path_uncorr)[0]

        if self.galaxy.stokes:
            cube_pbcorr = self.remove_stokes(cube_pbcorr)
            cube_uncorr = self.remove_stokes(cube_uncorr)

        # Get rid of nans
        cube_pbcorr.data[~np.isfinite(cube_pbcorr.data)] = 0
        cube_uncorr.data[~np.isfinite(cube_uncorr.data)] = 0

        return cube_pbcorr, cube_uncorr

    def splitCube(self, cube):
        """
        Split a cube into a cube containing the channels with emission and a cube containing the channels without.
        :param cube (HDU file): input HDU file containing the spectral cube and its header
        :return: two HDU files containing the cube channels with emission and the line-free channels, respectively
        """
        emiscube = cube.data[self.galaxy.start:self.galaxy.stop, :, :]
        noisecube = np.concatenate((cube.data[:self.galaxy.start, :, :], cube.data[self.galaxy.stop:, :, :]), axis=0)

        emiscube_hdu = fits.PrimaryHDU(emiscube, cube.header)
        emiscube_hdu.header['NAXIS3'] = emiscube.shape[0]

        noisecube_hdu = fits.PrimaryHDU(noisecube, cube.header)
        noisecube_hdu.header['NAXIS3'] = noisecube.shape[0]

        return emiscube_hdu, noisecube_hdu

    def innersquare(self, cube):
        """
        Get the central square (in spatial directions) of the spectral cube (useful for calculating the rms in a PB
        corrected spectral cube). Can be used for 2 and 3 dimensions, in the latter case the velocity axis is left
        unchanged.
        :param cube (2D or 3D array): input cube or 2D image
        :return: 2D or 3D array of the inner 1/8 of the cube in the spatial directions
        """

        start = int(len(cube[1]) * (7 / 16))
        stop = int(len(cube[1]) * (9 / 16))

        if len(cube.shape) == 3:
            return cube[:, start:stop, start:stop]
        elif len(cube.shape) == 2:
            return cube[start:stop, start:stop]
        else:
            raise AttributeError('Please provide a 2D or 3D array.')


    def clip(self, cube):
        """
        Creates a mask of the input cube where spaxels above the desired SNR are 1 and spaxels below the desired SNR
        are 0.
        :param cube (HDU file): spectral cube with which to create the mask
        :return: boolean mask with the spaxels above the provided level set to 1 and the spaxels below to 0
        """
        emiscube, noisecube = self.splitCube(cube)
        inner_noisecube = self.innersquare(noisecube.data)
        rms = np.nanstd(inner_noisecube)
        emiscube.data[emiscube.data < self.galaxy.cliplevel * rms] = 0
        emiscube.data[emiscube.data > self.galaxy.cliplevel * rms] = 1

        return emiscube.data.astype(bool)

    def prune_small_detections(self, cube, mask):
        """
        Mask structures in the spectral cube that are smaller than the desired size specified by "prune_by_npix" or
        "prune_by_fracbeam" in the galaxy parameters.
        :param cube (HDU file): the cube we are working on, to extract relevant information about the beam
        :param mask (3D array): the mask we have created thus far using the SUn clipping method
        :return: updated mask with the small detections masked out
        """

        if self.galaxy.prune_by_npix:
            prune_by_npix = self.galaxy.prune_by_npix
        else:
            res = cube.header['CDELT2']  # deg. / pix.
            bmaj_pix = cube.header['BMAJ'] / res  # deg. / (deg. / pix.)
            bmin_pix = cube.header['BMIN'] / res  # deg. / (deg. / pix.)
            beam_area_pix = np.pi * bmaj_pix * bmin_pix
            prune_by_npix = beam_area_pix * self.galaxy.prune_by_fracbeam

        labels, count = label(mask)
        for idx in np.arange(count) + 1:
            if (labels == idx).any(axis=0).sum() / idx < prune_by_npix:
                mask[labels == idx] = False

        return mask

    def expand_along_spatial(self, cube, mask):
        """
        Expand the mask along spatial dimensions by an amount provided by either "expand_by_npix" or
        "expand_by_fracbeam" in the galaxy parameters.
        :param cube (HDU file): cube that we are working on, to extract the relevant information from its header
        :param mask (3D array): mask that we have created so far with the Sun clipping method
        :return: updated, expanded mask
        """

        if self.galaxy.expand_by_npix:
            expand_by_npix = int(self.galaxy.expand_by_npix)
        else:
            res = cube.header['CDELT2']  # deg. / pix.
            bmaj = cube.header['BMAJ']  # deg.
            bmin = cube.header['BMIN']  # deg.
            beam_hwhm_pix = np.average([bmaj, bmin]) / res / 2  # deg. / (deg. / pix.)
            expand_by_npix = int(beam_hwhm_pix * self.galaxy.expand_by_fracbeam)

        structure = np.zeros([3, expand_by_npix * 2 + 1, expand_by_npix * 2 + 1])
        Y, X = np.ogrid[:expand_by_npix * 2 + 1, :expand_by_npix * 2 + 1]
        R = np.sqrt((X - expand_by_npix) ** 2 + (Y - expand_by_npix) ** 2)
        structure[1, :] = R <= expand_by_npix
        mask = binary_dilation(mask, iterations=1, structure=structure)

        return mask

    def expand_along_spectral(self, mask):
        """
        Expand the mask along the velocity direction as provided by "expand_by_nchan" in the galaxy parameters.
        :param mask: mask that we have created so far with the Sun clipping method
        :return: updated, expanded mask
        """
        for i in range(self.galaxy.expand_by_nchan):
            tempmask = np.roll(mask, shift=1, axis=0)
            tempmask[0, :] = False
            mask |= tempmask
            tempmask = np.roll(mask, shift=-1, axis=0)
            tempmask[-1, :] = False
            mask |= tempmask

        return mask

    def sun_method(self, emiscube, noisecube):
        """
        Apply the clipping method from Sun, possibly prune detections with small areas on the sky and/or expand the
        mask in the spatial/velocity directions.
        :param emiscube (HDU file): HDU containing the cube with only the channels with emission in them, from which
        the mask will be created
        :param noisecube (HDU file): HDU containing the cube with line-free channels, from which the rms will be
        estimated
        :return: mask with the same shape as "emiscube" where spaxels with a too low SNR are set to 0 according to the
        Sun method, and spaxels with a high enough SNR are set to 1.
        """

        # Check if the necessary parameters are provided
        if not (
                self.galaxy.nchan_low and self.galaxy.cliplevel_low and self.galaxy.nchan_high and
                self.galaxy.cliplevel_high):
            raise AttributeError('If you want to use Sun\'s method, please provide "nchan_low", "cliplevel_low", '
                                 '"nchan_high", and "cliplevel_high".')

        # Estimate the rms from the spatial inner part of the cube
        inner = self.innersquare(noisecube.data)
        rms = np.nanstd(inner)

        snr = emiscube.data / rms

        # Generate core mask
        mask_core = (snr > self.galaxy.cliplevel_high).astype(bool)
        for i in range(self.galaxy.nchan_high - 1):
            mask_core &= np.roll(mask_core, shift=1, axis=0)
        mask_core[:self.galaxy.nchan_high - 1] = False
        for i in range(self.galaxy.nchan_high - 1):
            mask_core |= np.roll(mask_core, shift=-1, axis=0)

        # Generate wing mask
        mask_wing = (snr > self.galaxy.cliplevel_low).astype(bool)
        for i in range(self.galaxy.nchan_low - 1):
            mask_wing &= np.roll(mask_wing, shift=1, axis=0)
        mask_wing[:self.galaxy.nchan_low - 1] = False
        for i in range(self.galaxy.nchan_low - 1):
            mask_wing |= np.roll(mask_wing, shift=-1, axis=0)

        # Dilate core mask inside wing mask
        mask = binary_dilation(mask_core, iterations=0, mask=mask_wing)

        # Prune detections with small projected areas on the sky
        if self.galaxy.prune_by_fracbeam or self.galaxy.prune_by_npix:
            mask = self.prune_small_detections(emiscube, mask)

        # Expand along spatial dimensions by a fraction of the beam FWHM
        if self.galaxy.expand_by_fracbeam or self.galaxy.expand_by_npix:
            mask = self.expand_along_spatial(emiscube, mask)

        # Expand along spectral dimension by a number of channels
        if self.galaxy.expand_by_nchan:
            mask = self.expand_along_spectral(mask)

        return mask

    def smooth_clip(self, cube):
        """
        Apply a Gaussian blur, using sigma = 4 in the velocity direction (seems to work best), to the uncorrected cube.
        The mode 'nearest' seems to give the best results.
        :return: (ndarray) mask to apply to the un-clipped cube
        """

        res = cube.header['CDELT2']  # deg/pixel
        bmaj = cube.header['BMAJ']  # degrees
        beam = bmaj / res  # beam size in pixels, use the major axis
        sigma = 1.5 * beam / np.sqrt(8. * np.log(2.))
        smooth_cube = ndimage.filters.gaussian_filter(cube.data, (4., sigma, sigma), order=0, mode='nearest')
        smooth_hdu = fits.PrimaryHDU(smooth_cube, cube.header)
        mask = self.clip(smooth_hdu)

        return mask

    def do_clip(self, cube_pbcorr, cube_uncorr):
        """
        Clip the array, either according to the Sun method (if self.sun == True, which is default) or the smooth
        clipping method from Dame.
        :param cube_pbcorr (HDU file): primary beam corrected spectral cube, which we want to clip
        :param cube_uncorr (HDU file): primary beam UNcorrected spectral cube, from which we want to make the mask
        :return: HDU file with the clipped, primary beam corrected spectral cube
        """
        # copy the non-PB corrected datacube
        cube_uncorr_copy = cube_uncorr.copy()

        # get part of the cube that has emission from the PB corrected datacube
        emiscube_pbcorr, noisecube_pbcorr = self.splitCube(cube_pbcorr)
        emiscube_uncorr, noisecube_uncorr = self.splitCube(cube_uncorr_copy)

        if self.sun:
            mask = self.sun_method(emiscube_uncorr, noisecube_pbcorr)
        else:
            mask = self.smooth_clip(cube_uncorr_copy)

        emiscube_pbcorr.data *= mask
        clipped_hdu = fits.PrimaryHDU(emiscube_pbcorr.data, cube_pbcorr.header)

        return clipped_hdu

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

        cube_pbcorr, cube_uncorr = self.readfits()
        cube = self.do_clip(cube_pbcorr, cube_uncorr)
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
        inner_cube = self.innersquare(mom1)

        if not self.galaxy.sysvel:
            sysvel = np.nanmean(inner_cube)
        mom1 -= sysvel

        mom0_hdu = fits.PrimaryHDU(mom0, self.new_header(cube_pbcorr.header))
        mom1_hdu = fits.PrimaryHDU(mom1, self.new_header(cube_pbcorr.header))
        mom2_hdu = fits.PrimaryHDU(mom2, self.new_header(cube_pbcorr.header))

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

        cube_pbcorr, cube_uncorr = self.readfits()

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