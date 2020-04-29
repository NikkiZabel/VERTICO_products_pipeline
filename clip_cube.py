from astropy.io import fits
import numpy as np
from scipy import ndimage
from scipy.ndimage import binary_dilation, label
from targets import galaxies
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.ned import Ned


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

    def innersquare(self, cube):
        """
        Get the central square (in spatial directions) of the spectral cube (useful for calculating the rms in a PB
        corrected spectral cube). Can be used for 2 and 3 dimensions, in the latter case the velocity axis is left
        unchanged.
        :param cube (2D or 3D array): input cube or 2D image
        :return: 2D or 3D array of the inner 1/8 of the cube in the spatial directions
        """

        if len(cube.shape) == 3:
            start_x = int(cube.shape[1] / 2 - 20)
            stop_x = int(cube.shape[1] / 2 + 20)
            start_y = int(cube.shape[2] / 2 - 20)
            stop_y = int(cube.shape[2] / 2 + 20)
            inner_square = cube[:, start_x:stop_x, start_y:stop_y]
            if (inner_square == inner_square).any():
                return inner_square
            else:
                return cube

        elif len(cube.shape) == 2:
            start_x = int(cube.shape[0] / 2 - 20)
            stop_x = int(cube.shape[0] / 2 + 20)
            start_y = int(cube.shape[1] / 2 - 20)
            stop_y = int(cube.shape[1] / 2 + 20)
            inner_square = cube[start_x:stop_x, start_y:stop_y]
            if (inner_square == inner_square).any():
                return inner_square
            else:
                return cube
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

    def cut_empty_rows(self, cube, noisecube=None):

        w = wcs.WCS(cube.header, naxis=2)
        centre_sky = SkyCoord(cube.header['OBSRA'], cube.header['OBSDEC'], unit=(u.deg, u.deg))
        centre_pix = wcs.utils.skycoord_to_pixel(centre_sky, w)

        beam = cube.header['BMAJ']  # deg
        res = cube.header['CDELT2']  # deg/pixel
        beam_pix = beam / res

        # Find empty rows
        empty_x = np.all(cube.data == 0, axis=(0, 2))

        # If the number of empty rows is larger than the beamsize, cut away rows and leave as many as the beam size
        if np.sum(empty_x) > beam_pix:
            beam_pix = int(np.round(beam_pix))
            idx_false = [i for i, x in enumerate(empty_x) if not x]
            first_false = idx_false[0]
            last_false = idx_false[-1]
            empty_x[first_false:last_false] = False  # Make sure empty rows in the middle of the data are not affected
            last_false += 1
            if first_false - beam_pix > 0:
                empty_x[first_false - beam_pix:first_false] = False
            else:
                empty_x[0:first_false] = False
            empty_x[last_false:last_false + beam_pix] = False

            cube.data = cube.data[:, ~empty_x, :]
            if noisecube:
                noisecube.data = noisecube.data[:, ~empty_x, :]

            # Adjust the central pixel in the image header correspondingly
            pix_shift = [i for i, x in enumerate(empty_x) if not x][1] - 1
            if empty_x[0]:
                cube.header['CRPIX2'] -= pix_shift
            cube.header['NAXIS2'] = cube.shape[1]
            if noisecube:
                noisecube.header['CRPIX2'] -= pix_shift
                noisecube.header['NAXIS2'] = noisecube.shape[1]

        if not noisecube:
            noisecube = np.array([0])

        return cube, noisecube

    def cut_empty_columns(self, cube, noisecube=None):

        if len(noisecube.shape) == 1:
            noisecube = None

        beam = cube.header['BMAJ']  # deg
        res = cube.header['CDELT2']  # deg/pixel
        beam_pix = beam / res

        # Find empty columns
        empty_y = np.all(cube.data == 0, axis=(0, 1))

        # If the number of empty columns is larger than the beamsize, cut away rows and leave as many as the beam size
        if np.sum(empty_y) > beam_pix:
            beam_pix = int(np.round(beam_pix))
            idx_false = [i for i, x in enumerate(empty_y) if not x]
            first_false = idx_false[0]
            last_false = idx_false[-1]
            empty_y[first_false:last_false] = False  # Make sure empty rows in the middle of the data are not affected
            last_false += 1
            if first_false - beam_pix > 0:
                empty_y[first_false - beam_pix:first_false] = False
            else:
                empty_y[0:first_false] = False
            empty_y[last_false:last_false + beam_pix] = False

            cube.data = cube.data[:, :, ~empty_y]
            if noisecube:
                noisecube.data = noisecube.data[:, :, ~empty_y]

            # Adjust the header
            pix_shift = [i for i, x in enumerate(empty_y) if not x][1] - 1
            if empty_y[0]:
                cube.header['CRPIX1'] -= pix_shift
            cube.header['NAXIS1'] = cube.shape[2]
            if noisecube:
                noisecube.header['CRPIX1'] -= pix_shift
                cube.header['NAXIS1'] = noisecube.shape[2]

        if not noisecube:
            noisecube = np.array([0])

        return cube, noisecube

    def centre_data(self, cube, noisecube=None):
        """
        Pad the image with zeros so that the centre of the galaxy overlaps with the centre of the
        cube
        :return:
        """

        # Read in central coordinates from NED
        ra = Ned.query_object(cube.header['OBJECT'])['RA'][0]
        dec = Ned.query_object(cube.header['OBJECT'])['DEC'][0]

        # Find the central pixel based on these coordinates and WCS info
        w = wcs.WCS(cube.header, naxis=2)
        centre_sky = SkyCoord(ra, dec, unit=(u.deg, u.deg))
        centre_pix = wcs.utils.skycoord_to_pixel(centre_sky, w)

        # The amount the centre needs to shift to overlap with the central coordinates
        shift_x = int(np.round(centre_pix[1] - cube.shape[1] / 2))
        shift_y = int(np.round(centre_pix[0] - cube.shape[2] / 2))

        # Pad the image with twice the amount it has to shift, so that the new centre overlaps with the coordinates
        if shift_x > 0:
            temp = np.zeros((cube.shape[0], cube.shape[1] + shift_x * 2, cube.shape[2]))
            temp[:, 0:cube.shape[1], :] = cube.data
            if noisecube:
                temp_noise = np.zeros((noisecube.shape[0], noisecube.shape[1] + shift_x * 2, noisecube.shape[2]))
                temp_noise[:, 0:noisecube.shape[1], :] = noisecube.data
        elif shift_x < 0:
            temp = np.zeros((cube.shape[0], cube.shape[1] + abs(shift_x) * 2, cube.shape[2]))
            temp[:, temp.shape[1] - cube.shape[1]:temp.shape[1], :] = cube.data
            if noisecube:
                temp_noise = np.zeros((noisecube.shape[0], noisecube.shape[1] + abs(shift_x) * 2, noisecube.shape[2]))
                temp_noise[:, temp_noise.shape[1] - noisecube.shape[1]:temp_noise.shape[1], :] = noisecube.data
        else:
            temp = cube.data
            if noisecube:
                temp_noise = noisecube.data

        # Same in the y-direction
        if shift_y > 0:
            cube_new = np.zeros((temp.shape[0], temp.shape[1], temp.shape[2] + shift_y * 2))
            cube_new[:, :, 0:temp.shape[2]] = temp
            if noisecube:
                noisecube_new = np.zeros((temp_noise.shape[0], temp_noise.shape[1], temp_noise.shape[2] + shift_y * 2))
                noisecube_new[:, :, 0:temp_noise.shape[2]] = temp_noise
        elif shift_y < 0:
            cube_new = np.zeros((temp.shape[0], temp.shape[1], temp.shape[2] + abs(shift_y) * 2))
            cube_new[:, :, cube_new.shape[2] - temp.shape[2]:cube_new.shape[2]] = temp
            if noisecube:
                noisecube_new = np.zeros((temp_noise.shape[0], temp_noise.shape[1], temp_noise.shape[2] + abs(shift_y) * 2))
                noisecube_new[:, :, noisecube_new.shape[2] - temp_noise.shape[2]:noisecube_new.shape[2]] = temp_noise
        else:
            cube_new = temp
            if noisecube:
                noisecube_new = temp_noise

        new_header = cube.header.copy()
        # Only change the CRPIX if the padding happens BEFORE the current CRPIX
        if shift_y < 0:
            new_header['CRPIX1'] = cube.header['CRPIX1'] - 2 * shift_y
        if shift_x < 0:
            new_header['CRPIX2'] = cube.header['CRPIX2'] - 2 * shift_x
        new_header['NAXIS1'] = cube_new.shape[2]
        new_header['NAXIS2'] = cube_new.shape[1]

        cube_hdu = fits.PrimaryHDU(cube_new, new_header)

        if not noisecube:
            noisecube_new = np.array([0])

        noisecube_hdu = fits.PrimaryHDU(noisecube_new, new_header)

        if noisecube:
            noisecube_hdu.header['NAXIS3'] = noisecube.header['NAXIS3']

        return cube_hdu, noisecube_hdu

    def preprocess(self, cube, noisecube=None):

        cube, noisecube = self.cut_empty_rows(cube, noisecube)
        cube, noisecube = self.cut_empty_columns(cube, noisecube)
        cube, noisecube = self.centre_data(cube, noisecube)

        return cube, noisecube

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

    def clip(self, smooth_cube):
        """
        Creates a mask of the input cube where spaxels above the desired SNR are 1 and spaxels below the desired SNR
        are 0.
        :param cube (HDU file): spectral cube with which to create the mask
        :return: boolean mask with the spaxels above the provided level set to 1 and the spaxels below to 0
        """

        cube_pbcorr, cube_uncorr = self.readfits()

        emiscube, noisecube = self.split_cube(cube_pbcorr)
        inner_noisecube = self.innersquare(noisecube.data)
        rms = np.nanstd(inner_noisecube)

        emiscube_smooth, noisecube_smooth = self.split_cube(smooth_cube)

        emiscube_smooth.data[emiscube_smooth.data < self.galaxy.cliplevel * rms] = 0
        emiscube_smooth.data[emiscube_smooth.data > self.galaxy.cliplevel * rms] = 1

        return emiscube_smooth.data.astype(bool)

    def prune_small_detections(self, cube, mask):
        """
        Mask structures in the spectral cube that are smaller than the desired size specified by "prune_by_npix" or
        "prune_by_fracbeam" in the galaxy parameters.
        :param cube (HDU file): the cube we are working on, to extract relevant information about the beam
        :param mask (3D array): the mask we have created thus far using the Sun clipping method
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
            if (labels == idx).any(axis=0).sum() < prune_by_npix:
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

    def smooth_mask(self, cube):
        """
        Apply a Gaussian blur, using sigma = 4 in the velocity direction (seems to work best), to the uncorrected cube.
        The mode 'nearest' seems to give the best results.
        :return: (ndarray) mask to apply to the un-clipped cube
        """

        beam = np.array([cube.header['BMAJ'], cube.header['BMIN']]) / cube.header['CDELT2']
        sigma = 1.5 * cube.header['BMAJ'] / cube.header['CDELT2']  # / np.sqrt(8. * np.log(2.))

        '''
        xpixels = cube.shape[1]
        ypixels = cube.shape[2]
        cent = [xpixels / 2, ypixels / 2]
        rot = 0
        dirfac = 1

        # Convolve with 1.5 the beam in the spatial direction
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
        '''
        #cube_spatial_smooth = np.empty(cube.shape)
        #for i in range(cube.shape[0]):
        #    cube_spatial_smooth[i, :, :] = convolve_fft(cube.data[i, :, :], psf)

        #smooth_cube = ndimage.filters.gaussian_filter1d(cube_spatial_smooth, 4, axis=0, order=0, mode='nearest')
        smooth_cube = ndimage.uniform_filter(cube.data, size=[4, sigma, sigma], mode='constant')  # mode='nearest'

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
            mask_hdu = fits.PrimaryHDU(mask.astype(int), cube_pbcorr.header)
            mask_hdu.writeto(self.savepath + 'mask_sun.fits', overwrite=True)
        else:
            mask = self.smooth_mask(cube_uncorr_copy)
            mask_hdu = fits.PrimaryHDU(mask.astype(int), cube_pbcorr.header)
            mask_hdu.writeto(self.savepath + 'mask_smooth.fits', overwrite=True)

        emiscube_pbcorr.data *= mask
        clipped_hdu = fits.PrimaryHDU(emiscube_pbcorr.data, cube_pbcorr.header)

        # Do some pre-processing to make the creation of the moments easier
        #pb = fits.open('/home/nikki/Documents/Data/VERTICO/ReducedData/' + str(self.galaxy.name) + '/' +
        #               str(self.galaxy.name) + '_7m_co21_pb_rebin.fits')[0]

        clipped_hdu, noisecube_hdu = self.preprocess(clipped_hdu, noisecube=None)

        if self.tosave:
            clipped_hdu.writeto(self.savepath + 'clipped_cube.fits', overwrite=True)
            if noisecube_hdu.shape[0] > 1:
                noisecube_hdu.writeto(self.savepath + 'trimmed_noisecube.fits', overwrite=True)

        return clipped_hdu, noisecube_hdu
