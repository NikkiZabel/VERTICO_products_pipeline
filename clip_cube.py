from astropy.io import fits
import numpy as np
from scipy import ndimage
from scipy.ndimage import binary_dilation, label
from targets import galaxies


class ClipCube:

    def __init__(self, galname, path_pbcorr, path_uncorr, savepath=None, sun=True, tosave=True):
        self.galaxy = galaxies(galname)
        self.path_pbcorr = path_pbcorr
        self.path_uncorr = path_uncorr
        self.savepath = savepath or './'
        self.sun = sun
        self.tosave = tosave

    @staticmethod
    def remove_stokes(cube):
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

    @staticmethod
    def innersquare(cube):
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

    def split_cube(self, cube):
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

    def clip(self, cube):
        """
        Creates a mask of the input cube where spaxels above the desired SNR are 1 and spaxels below the desired SNR
        are 0.
        :param cube (HDU file): spectral cube with which to create the mask
        :return: boolean mask with the spaxels above the provided level set to 1 and the spaxels below to 0
        """
        emiscube, noisecube = self.split_cube(cube)
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

    def do_clip(self):
        """
        Clip the array, either according to the Sun method (if self.sun == True, which is default) or the smooth
        clipping method from Dame.
        :param cube_pbcorr (HDU file): primary beam corrected spectral cube, which we want to clip
        :param cube_uncorr (HDU file): primary beam UNcorrected spectral cube, from which we want to make the mask
        :return: HDU file with the clipped, primary beam corrected spectral cube
        """
        # copy the non-PB corrected datacube

        cube_pbcorr, cube_uncorr = self.readfits()

        cube_uncorr_copy = cube_uncorr.copy()

        # get part of the cube that has emission from the PB corrected datacube
        emiscube_pbcorr, noisecube_pbcorr = self.split_cube(cube_pbcorr)
        emiscube_uncorr, noisecube_uncorr = self.split_cube(cube_uncorr_copy)

        if self.sun:
            mask = self.sun_method(emiscube_uncorr, noisecube_pbcorr)
        else:
            mask = self.smooth_clip(cube_uncorr_copy)

        emiscube_pbcorr.data *= mask
        clipped_hdu = fits.PrimaryHDU(emiscube_pbcorr.data, cube_pbcorr.header)

        if self.tosave:
            clipped_hdu.writeto(self.savepath + 'clipped_cube.fits', overwrite=True)

        return clipped_hdu