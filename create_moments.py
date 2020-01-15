from astropy.io import fits
import numpy as np
import scipy.ndimage as ndimage
from gal_params import *

class moment_maps:

    def __init__(self, galaxy, path, pbcor, dv, cliplevel, stokes=False, tosave=False):
        self.galaxy = galaxy
        self.path = path
        self.pbcor = pbcor
        self.stokes = stokes
        self.tosave = tosave
        self.start = 33
        self.stop = 73
        self.cliplevel = cliplevel
        self.dv = dv

    def makebeam(self, xpixels, ypixels, header, rot=0, cent=0):

        # extract relevant information from header (pixel size, beam size)
        res = header['CDELT2']  # degrees per pixel
        bmaj = header['BMAJ']  # degrees
        bmin = header['BMIN']  # degrees
        beam = np.array([bmaj, bmin]) / res  # beam size in pixels

        # convert the FWHM of the beam to the std dev in the Gaussian distribution
        sigma = beam / np.sqrt(8. * np.log(2.))

        if not cent: cent = [xpixels / 2., ypixels / 2]

        if np.tan(np.radians(rot)) == 0:
            dirfac = 1
        else:
            dirfac = np.sign(np.tan(np.radians(rot)))

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

    def get_beamsize(self, cube):
        # convert from Jy/beam to Jy using by calculating the integral of the psf
        x = cube.shape[1]  # sizes of the spatial axes of the cube
        y = cube.shape[2]
        psf = self.makebeam(x, y, cube.header)  # create the psf which is a 2D elliptical Gaussian
        return np.sum(psf)

    def readfits(self):
        '''
        Read in cube HDU and get rid of nans
        :param path: path to fits file
        :param pbcor: use pb corrected cube or not?
        :return: cube HDU with nans set to zero
        '''
        if self.pbcor:
            fitsfile = self.path + self.galaxy + '.co.image.pbcor.fits'
        else:
            fitsfile = self.path + self.galaxy + '/' + self.galaxy + '_7m+tp_co21_flat_round_k.fits'

        cube = fits.open(fitsfile)[0]

        if self.stokes:
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

        cube.data[np.isnan(cube.data)] = 0

        return cube

    def new_header(self, header):
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

    def splitCube(self, cube):
        '''
        :param cube: data cube
        :return: two cubes, one containing the channels with the line, and one containing the line-free channels
        '''
        emiscube = cube.data[self.start:self.stop, :, :]
        noisecube = np.concatenate((cube.data[:self.start, :, :], cube.data[self.stop:, :, :]), axis=0)

        emiscube_hdu = fits.PrimaryHDU(emiscube, cube.header)
        emiscube_hdu.header['NAXIS3'] = emiscube.shape[0]

        noisecube_hdu = fits.PrimaryHDU(noisecube, cube.header)
        noisecube_hdu.header['NAXIS3'] = noisecube.shape[0]

        return emiscube_hdu, noisecube_hdu

    def innersquare(self, cube, N):
        '''
        :param cube: data cube
        :param N: # dimensions of the data cube
        :return: the spatial inner part of the data cube, that contains no noise from the PB
        '''
        start = int(len(cube[1]) * (7. / 16.))
        stop = int(len(cube[1]) * (9. / 16.))
        if N == 3:
            return cube[:, start:stop, start:stop]
        elif N == 2:
            return cube[start:stop, start:stop]
        else:
            print("Enter a correct value for N")

    def get_rms(self, cube):
        '''
        :param cube: data cube
        :return: rms in the cube, based on the inner part to avoid PB contamination
        '''
        emiscube, noisecube = self.splitCube(cube)
        inner_noisecube = self.innersquare(noisecube.data, 3)

        return np.nanstd(inner_noisecube)

    def clip(self, cube):
        '''
        :param cube: data cube
        :return: hard-clipped version of the part of the cube containing the emission line
        '''
        emiscube, noisecube = self.splitCube(cube)
        sigma = self.get_rms(cube)
        emiscube.data[emiscube.data < self.cliplevel * sigma] = 0

        return emiscube.data

    def smoothclip(self, cube):
        '''
        :param cube: data cube
        :return: smooth-clipped version of the part of the data cube containing the emission line
        '''
        # copy the datacube
        cube_copy = cube.data.copy()

        # get part of the cube that has emission
        emiscube, noisecube = self.splitCube(cube)

        # extract relevant information from header (pixel size, beam size)
        res = cube.header['CDELT2']  # deg/pixel
        bmaj = cube.header['BMAJ']  # degrees
        beam = bmaj / res  # beam size in pixels, use the major axis

        # convert the FWHM of the beam to the std dev in the Gaussian distribution and use 1.5 times the beam size to convolve with
        sigma = 1.5 * beam / np.sqrt(8. * np.log(2.))

        # apply a Gaussian blur, using sigma = 4 in the velocity direction (seems to work best). The mode 'nearest' seems to give the best results.
        cube_smoothed = ndimage.filters.gaussian_filter(cube_copy, (4., sigma, sigma), order=0, mode='nearest')
        cube_smoothed_hdu = fits.PrimaryHDU(cube_smoothed, cube.header)

        # hard-clip the smoothed copy of the data cube
        clipped = self.clip(cube_smoothed_hdu)

        # mask the original cube with the smoothed cube
        emiscube.data[clipped == 0] = 0

        # give it the header of the original cube with updated info for the number of channels
        clipped_hdu = fits.PrimaryHDU(emiscube.data, cube.header)
        clipped_hdu.header['NAXIS3'] = emiscube.shape[0]

        return clipped_hdu

    def create_vel_array(self, cube):

        v_val = cube.header['CRVAL3'] / 1000.  # velocity in the reference channel, m/s to km/s
        v_step = cube.header['CDELT3'] / 1000.  # velocity step in each channel, m/s to km/s
        v_ref = cube.header['CRPIX3']  # location of the reference channel

        # Construct the velocity array
        vel_array = ((np.arange(0, len(cube.data[:, 0, 0])) - v_ref - 1 + self.start) * v_step + v_val)  # fits index starts from 1, start at the beginning of the emission cube
        vel_narray = np.tile(vel_array, (len(cube.data[0, 0, :]), len(cube.data[0, :, 0]), 1)).transpose()

        return vel_array, vel_narray

    def calc_moms(self):
        '''
        Calculate the moments 0,1, and 2 by calculating sum(I), sum(v*I)/sum(I), and sum(v*I^2)/sum(I)
        :param cube: raw data cube
        :param sigma: rms in raw data cube
        :param dV: channel width in km/s
        :return:
        '''
        rawcube = self.readfits()
        cube = self.smoothclip(rawcube)
        rms = self.get_rms(cube)

        sigma = rms * self.dv

        vel_array, vel_narray = self.create_vel_array(cube)
        np.savetxt(self.path + 'vel_array.txt', vel_array)

        beamsize = self.get_beamsize(cube)

        # calculate moment 0
        zeroth = np.sum((cube.data * self.dv), axis=0)

        # calculate moment 1
        first = np.sum(cube.data * vel_narray, axis=0) / np.sum(cube.data, axis=0)

        # calculate moment 2
        second = np.sqrt(np.sum(abs(cube.data) * (vel_narray - first) ** 2., axis=0) / np.sum(abs(cube.data), axis=0))

        # calculate the systemic velocity from the spatial inner part of the cube (to avoid PB effects)
        inner_cube = self.innersquare(first, 2)
        sysvel = np.mean(inner_cube[np.isfinite(inner_cube)])
        np.savetxt(self.path + 'sysvel.txt', [sysvel])
        first -= sysvel

        mom0_hdu = fits.PrimaryHDU(zeroth, self.new_header(rawcube.header))
        mom1_hdu = fits.PrimaryHDU(first, self.new_header(rawcube.header))
        mom2_hdu = fits.PrimaryHDU(second, self.new_header(rawcube.header))

        if self.tosave:
            mom0_hdu.writeto(self.path + 'moment0.fits', overwrite=True)
            mom1_hdu.writeto(self.path + 'moment1.fits', overwrite=True)
            mom2_hdu.writeto(self.path + 'moment2.fits', overwrite=True)

        return mom0_hdu, mom1_hdu, mom2_hdu