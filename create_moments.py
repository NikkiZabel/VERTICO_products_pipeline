from astropy.io import fits
import numpy as np
import scipy.ndimage as ndimage
from targets import galaxies
from scipy.ndimage import binary_dilation, label

class moment_maps:

    def __init__(self, galname, path_pbcorr, path_uncorr, savepath='./', sun=True, tosave=False):
        self.galaxy = galaxies(galname)
        self.path_pbcorr = path_pbcorr
        self.path_uncorr = path_uncorr
        self.savepath = savepath
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
        beam = np.average([bmaj, bmin]) / res  # Average beam FWHM in pixels

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
        :param cube:
        :param mask:
        :return:
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
        for i in range(self.galaxy.expand_by_nchan):
            tempmask = np.roll(mask, shift=1, axis=0)
            tempmask[0, :] = False
            mask |= tempmask
            tempmask = np.roll(mask, shift=-1, axis=0)
            tempmask[-1, :] = False
            mask |= tempmask

        return mask

    def sun_method(self, emiscube, noisecube):

        if not (
                self.galaxy.nchan_low and self.galaxy.cliplevel_low and self.galaxy.nchan_high and
                self.galaxy.cliplevel_high):
            raise AttributeError('If you want to use Sun\'s method, please provide "nchan_low", "cliplevel_low", '
                                 '"nchan_high", and "cliplevel_high".')

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
        smooth_cube = ndimage.filters.gaussian_filter(cube, (4., sigma, sigma), order=0, mode='nearest')
        smooth_hdu = fits.PrimaryHDU(smooth_cube, cube.header)
        mask = self.clip(smooth_hdu)

        return mask

    def do_clip(self, cube_pbcorr, cube_uncorr):
        '''
        :param cube: data cube
        :return: smooth-clipped version of the part of the data cube containing the emission line
        '''
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
        v_val = cube.header['CRVAL3'] / 1000.  # velocity in the reference channel, m/s to km/s
        v_step = cube.header['CDELT3'] / 1000.  # velocity step in each channel, m/s to km/s
        v_ref = cube.header['CRPIX3']  # location of the reference channel

        # Construct the velocity array
        vel_array = ((np.arange(0, len(cube.data[:, 0, 0])) - v_ref - 1 + self.galaxy.start) * v_step + v_val)  # fits index starts from 1, start at the beginning of the emission cube
        vel_narray = np.tile(vel_array, (len(cube.data[0, 0, :]), len(cube.data[0, :, 0]), 1)).transpose()
        vel_array_full = ((np.arange(0, len(cube.data[:, 0, 0])) - v_ref - 1) * v_step + v_val)

        return vel_array, vel_narray, vel_array_full

    def calc_moms(self):
        """
        Calculate the moments 0,1, and 2 by calculating sum(I), sum(v*I)/sum(I), and sum(v*I^2)/sum(I)
        :param cube: raw data cube
        :param sigma: rms in raw data cube
        :return:
        """
        cube_pbcorr, cube_uncorr = self.readfits()
        cube = self.do_clip(cube_pbcorr, cube_uncorr)
        vel_array, vel_narray, vel_fullarray = self.create_vel_array(cube)

        zeroth = np.sum((cube.data * abs(cube.header['CDELT3']) / 1000), axis=0)
        first = np.sum(cube.data * vel_narray, axis=0) / np.sum(cube.data, axis=0)
        second = np.sqrt(np.sum(abs(cube.data) * (vel_narray - first) ** 2., axis=0) / np.sum(abs(cube.data), axis=0))

        # calculate the systemic velocity from the spatial inner part of the cube (to avoid PB effects)
        inner_cube = self.innersquare(first)
        sysvel = np.mean(inner_cube[np.isfinite(inner_cube)])
        first -= sysvel

        mom0_hdu = fits.PrimaryHDU(zeroth, self.new_header(cube_pbcorr.header))
        mom1_hdu = fits.PrimaryHDU(first, self.new_header(cube_pbcorr.header))
        mom2_hdu = fits.PrimaryHDU(second, self.new_header(cube_pbcorr.header))

        if self.tosave:
            cube.writeto(self.savepath + 'clipped_cube.fits', overwrite=True)
            mom0_hdu.writeto(self.savepath + 'moment0.fits', overwrite=True)
            mom1_hdu.writeto(self.savepath + 'moment1.fits', overwrite=True)
            mom2_hdu.writeto(self.savepath + 'moment2.fits', overwrite=True)

        return cube, mom0_hdu, mom1_hdu, mom2_hdu, sysvel

    def spectrum(self):

        cube_pbcorr, cube_uncorr = self.readfits()

        # convert from Jy/beam to Jy using by calculating the integral of the psf
        psf = self.makebeam(cube_pbcorr.shape[1], cube_pbcorr.shape[2], cube_pbcorr.header)
        beamsize = np.sum(psf)

        cutout = cube_pbcorr.data[:, self.galaxy.centre_y - self.galaxy.size:self.galaxy.centre_y + self.galaxy.size,
                 self.galaxy.centre_x - self.galaxy.size:self.galaxy.centre_x + self.galaxy.size]

        cutout = cube_pbcorr.data

        # Make this work if necessary
        #if custom_region:
            #region = pyregion.open(path + 'ds9.reg')
            #mask = region.get_mask(hdu=mom0_hdu)
            #mask_3d = np.tile(mask, (len(cube[:, 0, 0]), 1, 1))
            #cutout = np.where(mask_3d, cube, 0)

        spectrum = np.nansum(cutout, axis=(1, 2))

        # Estimate the rms in the spectrum
        #emis, noise = self.splitCube(cutout, self.galaxy.start, self.galaxy.stop)
        #rms = np.std(np.sum(noise, axis=(1, 2)))
        #np.savetxt(path + 'specnoise.txt', [rms / beamsize])

        return spectrum / beamsize