from astropy.io import fits
import numpy as np
import scipy.ndimage as ndimage
from targets import galaxies
from clip_cube import ClipCube
from photutils import EllipticalAnnulus
from photutils import aperture_photometry
from astropy import wcs
from astropy.stats import mad_std
import os


class MomentMaps:

    def __init__(self, galname, path_pbcorr, path_uncorr, savepath=None, sun=True, tosave=False, sample=None, redo_clip=False):
        self.galaxy = galaxies(galname, sample)
        self.path_pbcorr = path_pbcorr
        self.path_uncorr = path_uncorr
        self.savepath = savepath or './'
        self.sun = sun
        self.tosave = tosave
        self.sample = sample
        self.redo_clip = redo_clip
        self.galaxy.start, self.galaxy.stop = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                    savepath=self.savepath, tosave=self.tosave, sample=self.sample).do_clip(get_chans=True)


    def pixel_size_check(self, header, key="CDELT1", expected_pix_size=2, raise_exception=True):
        """
        Check if the pixel size is what is expected (modulo some floating point 
                                                     error).

        Parameters
        ----------
        header : FITS header
            Of the spectral cube of the galaxy in question.
        key : header keyword, optional
            Keyword containing the pixel size information. The default is 
            "CDELT1".
        expected_pix_size : float or int, optional
            The size the pixel should have. The default is 2.
        raise_exception : bool, optional
            Whether to raise an exception if the pixel size is incorrect. The 
            default is True.

        Raises
        ------
        Exception
            Raised if the pixel size of the data is not as expected.

        Returns
        -------
        None.

        """

        pixel_size_arcsec = abs(header[key] * 3600)
        pixel_size_error_message = self.galaxy.name + ": " + key + " = " + str(pixel_size_arcsec) + "arcsec" + \
                                   " not " + str(expected_pix_size) + "arcsec"
        if not np.max(np.abs(pixel_size_arcsec - expected_pix_size)) < 1e-6:
            if raise_exception == True:
                raise Exception(pixel_size_error_message)
            else:
                print(pixel_size_error_message)
                pass


    def calc_noise_in_cube(self,
            masking_scheme='simple', mask=None,
            spatial_average_npix=None, spatial_average_nbeam=5.0,
            spectral_average_nchan=5, verbose=False):
        """
        From Jiayi Sun's script: 
            https://github.com/astrojysun/Sun_Astro_Tools/blob/master/sun_astro_tools/spectralcube.py

        Estimate rms noise in a (continuum-subtracted) spectral cube.
        Parameters
        ----------
        masking_scheme : {'simple', 'user'}, optional
            Scheme for flagging signal in the cube. 'simple' means to flag
            all values above 3*rms or below -3*rms (default scheme);
            'user' means to use the user-specified mask (i.e., `mask`).
        mask : `np.ndarray` object, optional
            User-specified signal mask (this parameter is ignored if
            `masking_scheme` is not 'user')
        spatial_average_npix : int, optional
            Size of the spatial averaging box, in terms of pixel number
            If not None, `spatial_average_nbeam` will be ingored.
            (Default: None)
        spatial_average_nbeam : float, optional
            Size of the spatial averaging box, in the unit of beam FWHM
            (Default: 5.0)
        spectral_average_nchan : int, optional
            Size of the spectral averaging box, in terms of channel number
            (Default: 5)
        verbose : bool, optional
            Whether to print the detailed processing information in terminal
            Default is to not print.

        Returns
        -------
        rmscube : SpectralCube object
            Spectral cube containing the rms noise at each ppv location
        """

        if masking_scheme not in ['simple', 'user']:
            raise ValueError("'masking_scheme' should be specified as"
                             "either 'simple' or 'user'")
        elif masking_scheme == 'user' and mask is None:
            raise ValueError("'masking_scheme' set to 'user', yet "
                             "no user-specified mask found")

        # Calculate the noise in the line-free part of the PB UN/CORRECTED cube
        cube_pbcorr, cube_uncorr = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                    savepath=self.savepath, tosave=self.tosave, sample=self.sample).readfits()

        # Centre and clip empty rows and columns to get it in the same shape as the other products
        if not self.redo_clip:
            if os.path.exists(self.savepath + 'noise_subcube_slab.fits'):
                noisecube = fits.read(self.savepath + 'noise_subcube_slab.fits')[0]
            else:
                _, noisecube = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                        savepath=self.savepath, tosave=self.tosave, sample=self.sample).do_clip(clip_also_nat='noise')
        else:
            _, noisecube = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                    savepath=self.savepath, tosave=self.tosave, sample=self.sample).do_clip(clip_also_nat='noise')

        # extract negative values (only needed if masking_scheme='simple')
        if masking_scheme == 'simple':
            if verbose:
                print("Extracting negative values...")
            negdata = np.where(noisecube.data < 0, noisecube.data, np.nan)
            negdata = np.stack([negdata, -1 * negdata], axis=-1)
        else:
            negdata = None

        # find rms noise as a function of channel
        if verbose:
            print("Estimating rms noise as a function of channel...")
        if masking_scheme == 'user':
            mask_v = mask
        elif masking_scheme == 'simple':
            rms_v = mad_std(negdata, axis=(1, 2, 3), ignore_nan=True)
            uplim_v = (3 * rms_v).reshape(-1, 1, 1)
            lolim_v = (-3 * rms_v).reshape(-1, 1, 1)
            mask_v = (((noisecube.data - uplim_v) < 0) &
                      ((noisecube.data - lolim_v) > 0))

        rms_v = mad_std(np.where(mask_v == 1, noisecube.data, np.nan), axis=(1, 2), ignore_nan=True)
        rms_v = ndimage.generic_filter(rms_v, np.nanmedian,
                               mode='constant', cval=np.nan,
                               size=spectral_average_nchan)

        # find rms noise as a function of sightline
        if verbose:
            print("Estimating rms noise as a function of sightline...")
        if masking_scheme == 'user':
            mask_s = mask
        elif masking_scheme == 'simple':
            rms_s = mad_std(negdata, axis=(0, 3), ignore_nan=True)
            uplim_s = 3 * rms_s
            lolim_s = -3 * rms_s
            mask_s = (((noisecube.data - uplim_s) < 0) &
                      ((noisecube.data - lolim_s) > 0))
        rms_s = mad_std(np.where(mask_s == 1, noisecube.data, np.nan), axis=0, ignore_nan=True)
        if spatial_average_npix is None:
            beamFWHM_pix = cube_pbcorr.header['BMAJ'] / cube_pbcorr.header['CDELT2']
            beamFWHM_pix = np.max([beamFWHM_pix, 3])
            spatial_average_npix = int(spatial_average_nbeam * beamFWHM_pix)
        rms_s = ndimage.generic_filter(rms_s, np.nanmedian, mode='constant', cval=np.nan, size=spatial_average_npix)

        # create rms noise cube from the tensor product of rms_v and rms_s
        if verbose:
            print("Creating rms noise cube (direct tensor product)...")
        rmscube = np.einsum('i,jk', rms_v, rms_s)

        # correct the normalization of the rms cube
        if masking_scheme == 'user':
            mask_n = mask
        elif masking_scheme == 'simple':
            rms_n = mad_std(negdata, ignore_nan=True)
            uplim_n = 3 * rms_n
            lolim_n = -3 * rms_n
            mask_n = (((noisecube.data - uplim_n) < 0) &
                      ((noisecube.data - lolim_n) > 0))
        rms_n = mad_std(noisecube.data[mask_n])
        rmscube /= rms_n

        # Write as FITS file
        rmscube_hdu = fits.PrimaryHDU(rmscube, noisecube.header)
        rmscube_hdu.writeto(self.savepath + 'rms_cube.fits', overwrite=True)

        return rmscube_hdu


    def makebeam(self, xpixels, ypixels, header, rot=0, cent=0):
        """
        Creates the psf from the beam size given a custom location in the image 
        and rotation, assuming a 2D Gaussian.

        Parameters
        ----------
        xpixels : float or int
            number of pixels in the x-axis.
        ypixels : float or int
            number of pixels in the y-axis.
        header : FITS header
            FITS header of the corresponding spectral cube.
        rot : float or int, optional
            Number of degrees over which the psf should be rotated. 
            The default is 0.
        cent : list or ndarray of length 2, optional
            Alternative centre of the psf. The default is 0.

        Returns
        -------
        psf : 2D numpy array
            psf of the image

        """

        # Extract relevant information from header
        res = header['CDELT2']  # degrees per pixel
        bmaj = header['BMAJ']  # degrees
        bmin = header['BMIN']  # degrees
        beam = np.array([bmaj, bmin]) / res  # Beam FWHM in pixels

        # Convert the FWHM of the beam to the std dev in the Gaussian distribution
        sigma = beam / np.sqrt(8. * np.log(2.))

        # If centre is not defined, make it the centre of the image
        if not cent: cent = [xpixels / 2, ypixels / 2]

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
        Change the 3D header to the corresponding 2D one.

        Parameters
        ----------
        header : FITS header
            Corresponding to the data cube.

        Returns
        -------
        header : FITS header
            Corresponding to the 2D image created from the 3D data cube.

        """

        header = header.copy()

        try:
            header.pop('PC3_1')
            header.pop('PC3_2')
            header.pop('PC1_3')
            header.pop('PC2_3')
            header.pop('PC3_3')
        except:
            try:
                header.pop('PC03_01')
                header.pop('PC03_03')
                header.pop('PC03_02')
                header.pop('PC01_03')
                header.pop('PC02_03')
            except:
                pass

        header.pop('CTYPE3')
        header.pop('CRVAL3')
        header.pop('CDELT3')
        header.pop('CRPIX3')
        try:
            header.pop('CUNIT3')
        except:
         pass
        header.pop('NAXIS3')
        try:
            header.pop('OBSGEO-Z')
        except:
            pass

        header['NAXIS'] = 2
        try:
            if header['WCSAXES'] == 3:
                header['WCSAXES'] = 2
        except:
            pass
        try:
            header.pop('CROTA3')
        except:
            pass

        return header


    def create_vel_array(self, cube):
        """
        Creates the velocity array corresponding to the spectral axis
        of the cube in km/s.

        Parameters
        ----------
        cube : FITS file
            The spectral cube from which we will make the velocity array.

        Raises
        ------
        KeyError
            Raised if the units are something other than velocity or frequency.

        Returns
        -------
        vel_array : 1D numpy array
            Contains the velocities corresponding to the cube spectral axis.
            Only contains values from the velocity corresponding to the start of the 
        spectral line data onwards
        vel_narray : 3D numpy array
            Same as vel_array but in the shape of the data cube (i.e. tiled in
                                                    the spatial dimensions).
        vel_array_full : 1D numpy array
            Same as vel_array but for the entire spectral axis.

        """

        cube_orig, _ = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun, savepath=self.savepath,
                        tosave=self.tosave, sample=self.sample).readfits()

        v_ref = cube_orig.header['CRPIX3']  # Location of the reference channel
        if cube_orig.header['CTYPE3'] == 'VRAD' or cube_orig.header['CTYPE3'] == 'VELOCITY':
            v_val = cube_orig.header['CRVAL3'] / 1000  # Velocity in the reference channel, m/s to km/s
            v_step = cube_orig.header['CDELT3'] / 1000  # Velocity step in each channel, m/s to km/s
        elif cube_orig.header['CTYPE3'] == 'FREQ':
            if self.sample == 'viva' or self.sample == 'things':
                v_val = 299792.458 * (1 - (cube_orig.header['CRVAL3'] / 1e9) / 1.420405752)
                v_shift = 299792.458 * (1 - ((cube_orig.header['CRVAL3'] + cube_orig.header['CDELT3']) / 1e9) / 1.420405752)
                v_step = - (v_val - v_shift)
            else:
                v_val = 299792.458 * (1 - (cube_orig.header['CRVAL3'] / 1e9) / 230.538000)
                v_shift = 299792.458 * (1 - ((cube_orig.header['CRVAL3'] + cube_orig.header['CDELT3']) / 1e9) / 230.538000)
                v_step = - (v_val - v_shift)
        else:
            raise KeyError('Pipeline cannot deal with these units yet.')

        # Construct the velocity arrays (keep in mind that fits-files are 1 indexed)
        vel_array = (np.arange(0, len(cube.data[:, 0, 0])) - v_ref + 1 + self.galaxy.start) * v_step + v_val
        vel_narray = np.tile(vel_array, (len(cube.data[0, 0, :]), len(cube.data[0, :, 0]), 1)).transpose()
        vel_array_full = (np.arange(0, len(cube.data[:, 0, 0])) - v_ref + 1) * v_step + v_val

        return vel_array, vel_narray, vel_array_full


    def add_clipping_keywords(self, header):
        """
        Add information to the header specifying details about the clipping.

        Parameters
        ----------
        header : FITS header
            Header of the cube that was clipped.

        Returns
        -------
        header : FITS header
            Header of the cube with clipping-related keywords added.

        """
        if self.sun:
            try:
                header.add_comment('Cube was clipped using the Sun+18 masking method', before='BUNIT')
            except:
                header.add_comment('Cube was clipped using the Sun+18 masking method', after='NAXIS2')
            header['CLIPL_L'] = self.galaxy.cliplevel_low
            header.comments['CLIPL_L'] = 'S/N threshold specified for the "wing mask"'
            header['CLIPL_H'] = self.galaxy.cliplevel_high
            header.comments['CLIPL_H'] = 'S/N threshold specified for the "core mask"'
            header['NCHAN_L'] = self.galaxy.nchan_low
            header.comments['NCHAN_L'] = '# of consecutive channels specified for the "core mask"'
            header['NCHAN_H'] = self.galaxy.nchan_high
            header.comments['NCHAN_H'] = '# of consecutive channels specified for the "wing mask"'
        else:
            try:
                header.add_comment('Cube was clipped using the Dame11 masking method', before='BUNIT')
            except:
                header.add_comment('Cube was clipped using the Dame11 masking method', after='NAXIS2')
            header['CLIPL'] = self.galaxy.cliplevel
            header.comments['CLIPL'] = 'SNR used for clip (Dame11)'

        header['CLIP_RMS'] = self.uncertainty_maps(calc_rms=True)
        header.comments['CLIP_RMS'] = 'rms noise level used in masking (K km/s)'

        return header


    def calc_moms(self, units='M_Sun/pc^2', alpha_co=5.4):
        """
        Clip the spectral cube according to the desired method (either the 
        method Pythonified by Jiayi Sun or the more basic smooth + clip 
        strategy from Dame+11), and create moment 0, 1, and 2 maps. Saves maps 
        as fits files if so desired. Also calculate the systemic velocity from 
        the moment 1 map.

        Parameters
        ----------
        units : str, optional
            Preferred units (either 'K km/s' or 'M_Sun/pc^2'). The default is 
                             'M_Sun/pc^2'.
        alpha_co : float, optional
            In case units == 'M_Sun/pc^2', multiply the moment 0 map by this 
            factor to obtain these units. The default is 5.4, which is the 
            value for CO(2-1) quoted in https://arxiv.org/pdf/1805.00937.pdf.

        Raises
        ------
        AttributeError
            Raised if units is set to anything other than 'K km/s' or 
            'M_Sun/pc^2'.

        Returns
        -------
        cube : FITS file
            3D spectral line cube for which the moment maps will be calculated.
        mom0_hdu : FITS file
            Contains the moment 0 map + corresponding header.
        mom1_hdu : FITS file
            Contains the moment 1 map + corresponding header.
        mom2_hdu : FITS file
            Contains the moment 2 map + corresponding header.
        sysvel : float
            Systemic velocity of the gas in the system in km/s.

        """

        if self.redo_clip:
            cube, _ = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun, savepath=self.savepath,
                            tosave=self.tosave, sample=self.sample).do_clip()
        elif os.path.exists(self.savepath + 'subcube_slab.fits'):
            cube = fits.open(self.savepath + 'subcube_slab.fits')[0]
        else:
            cube, _ = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                               savepath=self.savepath, tosave=self.tosave, sample=self.sample).do_clip()

        vel_array, vel_narray, vel_fullarray = self.create_vel_array(cube)

        if cube.header['CTYPE3'] == 'VRAD' or cube.header['CTYPE3'] == 'VELOCITY':
            mom0 = np.sum((cube.data * abs(cube.header['CDELT3']) / 1000), axis=0)
        elif cube.header['CTYPE3'] == 'FREQ':
            if self.sample == 'viva' or self.sample == 'things':
                v_val = 299792.458 * (1 - (cube.header['CRVAL3'] / 1e9) / 1.420405752)
                v_shift = 299792.458 * (1 - ((cube.header['CRVAL3'] + cube.header['CDELT3']) / 1e9) / 1.420405752)
                v_step = - (v_val - v_shift)
            else:
                v_val = 299792.458 * (1 - (cube.header['CRVAL3'] / 1e9) / 230.538000)
                v_shift = 299792.458 * (1 - ((cube.header['CRVAL3'] + cube.header['CDELT3']) / 1e9) / 230.538000)
                v_step = - (v_val - v_shift)
            mom0 = np.sum(cube.data * abs(v_step), axis=0)
        else:
            raise AttributeError("Can't deal with these units yet.")

        if units == 'M_Sun/pc^2':
            #mom0 = mom0 / cube.header['JTOK'] * 91.7 * alpha_co * (cube.header['BMAJ'] * 3600 * cube.header[
            #    'BMIN'] * 3600) ** (-1) / 4

            if self.sample == 'viva' or self.sample == 'things':
                coldens_atom_cm = mom0 * 1.10e24 / (cube.header['BMAJ'] * 3600 * cube.header['BMIN'] * 3600)
                Msol_to_matom = 1.187883838e57
                pc_to_cm = 9.521e36
                coldens_Msol_pc = coldens_atom_cm / Msol_to_matom * pc_to_cm
                mom0 = coldens_Msol_pc
            elif self.sample == None:
                xco = 1.23361968e+20
                mom0 *= 91.9 * xco / (cube.header['BMAJ'] * 3600 * cube.header['BMIN'] * 3600) * 1.6014457E-20
            else:
                mom0 *= alpha_co
        elif units == 'K km/s':
            #if self.sample == 'viva':
            #    mom0 = mom0 * 1.222e3 / (cube.header['BMAJ'] * 3600 * cube.header['BMIN'] * 3600) / (1.420405751) ** 2
            pass
        else:
            raise AttributeError('Please choose between "K km/s" and "M_Sun/pc^2"')
        mom1 = np.sum(cube.data * vel_narray, axis=0) / np.sum(cube.data, axis=0)
        mom2 = np.sqrt(np.sum(abs(cube.data) * (vel_narray - mom1) ** 2, axis=0) / np.sum(abs(cube.data), axis=0))

        # Calculate the systemic velocity from the spatial inner part of the cube (to avoid PB effects)
        inner_cube = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, savepath=self.savepath,
                              tosave=self.tosave, sample=self.sample).innersquare(mom1)

        from matplotlib import pyplot as plt
        plt.imshow(inner_cube)

        sysvel = np.nanmean(inner_cube)# + self.galaxy.sysvel_offset
        mom1 -= sysvel

        mom0_hdu = fits.PrimaryHDU(mom0, self.new_header(cube.header))
        mom1_hdu = fits.PrimaryHDU(mom1, self.new_header(cube.header))
        mom2_hdu = fits.PrimaryHDU(mom2, self.new_header(cube.header))

        #self.pixel_size_check(header=mom0_hdu.header)

        # Change or add any (additional) keywords to the headers
        if units == 'M_Sun/pc^2':
            mom0_hdu.header['BTYPE'] = 'Column density'
            if not self.sample == 'viva' or self.sample == 'things':
                mom0_hdu.header.comments['BTYPE'] = 'Total molecular gas (H_2 + He)'
        else:
            mom0_hdu.header['BTYPE'] = 'Integrated intensity'

        mom1_hdu.header['BTYPE'] = 'Velocity'
        mom2_hdu.header['BTYPE'] = 'Linewidth'
        mom0_hdu.header['BUNIT'] = units; mom0_hdu.header.comments['BUNIT'] = ''
        mom1_hdu.header['BUNIT'] = 'km/s'; mom1_hdu.header.comments['BUNIT'] = ''
        mom2_hdu.header['BUNIT'] = 'km/s'; mom2_hdu.header.comments['BUNIT'] = ''
        if not self.sample == 'viva' or self.sample == 'things':
            mom0_hdu.header['ALPHA_CO'] = alpha_co; #mom0_hdu.header.comments['ALPHA_CO'] = 'Assuming a line ratio of 0.8'
        mom1_hdu.header['SYSVEL'] = sysvel; mom1_hdu.header.comments['SYSVEL'] = 'km/s'

        self.add_clipping_keywords(mom0_hdu.header)
        self.add_clipping_keywords(mom1_hdu.header)
        self.add_clipping_keywords(mom2_hdu.header)

        if self.tosave:
            if units == 'M_Sun/pc^2':
                mom0_hdu.writeto(self.savepath + 'mom0_Msolpc-2.fits', overwrite=True)
            if units == 'K km/s':
                if self.sample == 'viva' or self.sample == 'things':
                    mom0_hdu.writeto(self.savepath + 'mom0_Jyb-1kms-1.fits', overwrite=True)
                else:
                    mom0_hdu.writeto(self.savepath + 'mom0_Kkms-1.fits', overwrite=True)
            mom1_hdu.writeto(self.savepath + 'mom1.fits', overwrite=True)
            mom2_hdu.writeto(self.savepath + 'mom2.fits', overwrite=True)

        return cube, mom0_hdu, mom1_hdu, mom2_hdu, sysvel


    def PVD(self, axis='major', find_angle=False, check_slit=False):
        """
        Creates a position-velocity diagram from the clipped data cube.

        Parameters
        ----------
        axis : str, optional
            Specifies whether the PVD will be taken along the minor or major 
            axis. The default is 'major'.
        find_angle : bool, optional
            If True, display 2D image of the data rotated by "rot_angle" 
            degrees, for quick visual determination of the position angle. The 
            default is False.
        check_slit : bool, optional
            If True, shows a 2D image of the data with the area chosen for PVD
            calculation overplotted. The default is False.

        Raises
        ------
        AttributeError
            Raised if anything other than 'major' or 'minor' is provided for 
            "axis".

        Returns
        -------
        pvd_hdu : FITS file
            contains the PVD and corresponding header.

        """

        clipped_cube, _, _, _, sysvel = self.calc_moms()

        res = clipped_cube.header['CDELT2']  # degrees per pixel
        beampix = clipped_cube.header['BMAJ'] / res  # beam size in pixels
        slitsize = np.ceil(beampix / 2)

        # Rotate the cube along the spatial axes, so that the galaxy lies horizontal
        if axis == 'major':
            rot_angle = self.galaxy.angle + 90
        elif axis == 'minor':
            rot_angle = self.galaxy.angle
        else:
            raise AttributeError('Please choose between "major" and "minor" for the "axis" keyword')

        cube_rot = ndimage.interpolation.rotate(clipped_cube.data, rot_angle, axes=(1, 2), reshape=True)

        if find_angle:
            from matplotlib import pyplot as plt
            plt.imshow(np.sum(cube_rot, axis=0))
            return

        # Define a slit around the centre of the galaxy with a width of the beam size (or use the full width of the galaxy)
        if self.galaxy.full_width:
            slit = cube_rot
        else:
            slit = cube_rot[:, int(cube_rot.shape[1] / 2 - slitsize):int(cube_rot.shape[1] / 2 + slitsize), :]

        # Collapse along the slit to create the PV diagram
        PV = np.sum(slit, axis=1)

        if check_slit:
            from matplotlib import pyplot as plt
            test = cube_rot.copy()
            test[:, int(cube_rot.shape[1] / 2 - slitsize): int(cube_rot.shape[1] / 2 + slitsize), :] = np.nan
            plt.imshow(np.sum(test, axis=0))

        # There is a lot of garbage because of the interpolation used by the rotation function, define a lower limit to get rid of that
        PV[PV < 1e-3] = 0

        # Remove empty rows and columns
        PV = PV[:, ~np.all(PV == 0, axis=0)]
        PV = PV[~np.all(PV == 0, axis=1), :]

        # Find the sky coordinates of the central pixel
        w = wcs.WCS(clipped_cube.header, naxis=2)
        centre_sky = wcs.utils.pixel_to_skycoord(clipped_cube.shape[0] / 2, clipped_cube.shape[1] / 2, wcs=w)

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
        pvd_header['FULL_WIDTH'] = self.galaxy.full_width
        pvd_header.comments['FULL_WIDTH'] = 'True if all emission used.'
        try:
            pvd_header['OBJECT'] = clipped_cube.header['OBJECT']
        except:
            pass
        pvd_header['EQUINOX'] = 2000
        pvd_header['RADESYS'] = 'FK5'
        try:
            pvd_header['LONPOLE'] = clipped_cube.header['LONPOLE']
            pvd_header['LATPOLE'] = clipped_cube.header['LATPOLE']
        except:
            pass
        pvd_header['CTYPE1'] = 'OFFSET'
        pvd_header['CRVAL1'] = 0
        pvd_header['CDELT1'] = clipped_cube.header['CDELT1'] / (cube_rot.shape[2] / cube_rot.shape[2])
        pvd_header['CRPIX1'] = np.ceil(PV.shape[1] / 2)
        try:
            pvd_header['CUNIT1'] = clipped_cube.header['CUNIT1']
        except:
            pass
        pvd_header['CTYPE2'] = 'VRAD'
        pvd_header['CRVAL2'] = clipped_cube.header['CRVAL3']
        pvd_header['CDELT2'] = clipped_cube.header['CDELT3']
        pvd_header['CRPIX2'] = clipped_cube.header['CRPIX3']
        try:
            pvd_header['CUNIT2'] = 'km/s'
        except:
            pass
        pvd_header['PC1_1'] = pvd_header['CDELT1']
        try:
            pvd_header['PC2_1'] = clipped_cube.header['PC2_1']
            pvd_header['PC1_2'] = clipped_cube.header['PC1_2']
        except:
            pass
        pvd_header['PC2_2'] = pvd_header['CDELT2']
        pvd_header['CENTR_PIX'] = '(' + str(clipped_cube.shape[0] / 2) + ', ' + str(clipped_cube.shape[1] / 2) + ')'
        pvd_header.comments['CENTR_PIX'] = 'Central pix used for rot. + loc. slit'
        pvd_header['CENTR_PIX_SKY'] = '(' + str(np.round(centre_sky.ra.deg, 2)) + ', ' + str(np.round(centre_sky.dec.deg, 2)) + ')'
        pvd_header.comments['CENTR_PIX_SKY'] = 'Central pix in sky coords (deg)'
        pvd_header['SYSVEL'] = sysvel + self.galaxy.sysvel_offset
        try:
            pvd_header['RESTFRQ'] = clipped_cube.header['RESTFRQ']
            pvd_header.comments['RESTFRQ'] = clipped_cube.header.comments['RESTFRQ']
        except:
            pass
        try:
            pvd_header['SPECSYS'] = clipped_cube.header['SPECSYS']
            pvd_header.comments['SPECSYS'] = clipped_cube.header.comments['SPECSYS']
        except:
            pass
        try:
            pvd_header['ALTRVAL'] = clipped_cube.header['ALTRVAL']
            pvd_header['ALTRPIX'] = clipped_cube.header['ALTRPIX']
        except:
            pass
        try:
            pvd_header['VELREF'] = clipped_cube.header['VELREF']
        except:
            pass
        try:
            pvd_header['USEWEIGH'] = clipped_cube.header['USEWEIGH']
        except:
            pass
        try:
            pvd_header['JTOK'] = clipped_cube.header['JTOK']
        except:
            pass
        try:
            pvd_header['OBSRA'] = clipped_cube.header['OBSRA']
            pvd_header['OBSDEC'] = clipped_cube.header['OBSDEC']
        except:
            pass
        try:
            pvd_header['OBSGEO-X'] = clipped_cube.header['OBSGEO-X']
            pvd_header['OBSGEO-Y'] = clipped_cube.header['OBSGEO-Y']
            pvd_header['OBSGEO-Z'] = clipped_cube.header['OBSGEO-Z']
        except:
            pass
        pvd_header['DISTANCE'] = self.galaxy.distance
        pvd_header.comments['DISTANCE'] = 'Mpc'
        try:
            pvd_header['ORIGIN'] = clipped_cube.header['ORIGIN']
        except:
            pass
        if self.sample == 'viva' or self.sample == 'things':
            pvd_header['BUNIT'] = 'Jy/b km/s'
        else:
            pvd_header['BUNIT'] = 'K km/s'
        self.add_clipping_keywords(pvd_header)

        pvd_hdu = fits.PrimaryHDU(PV, pvd_header)

        if self.tosave:
            if axis == 'major':
                pvd_hdu.writeto(self.savepath + 'pvd_major.fits', overwrite=True)
            if axis == 'minor':
                pvd_hdu.writeto(self.savepath + 'pvd_minor.fits', overwrite=True)

        return pvd_hdu


    def spectrum(self, useclipped=False):
        """
        Calculate the spectrum from the spectral cube.
        :return: array containing the spectrum in whatever units the cube is in, without the beam^-1 (so probably K or
        (m)Jy)
        """

        cube_pbcorr, cube_uncorr = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun, savepath=self.savepath,
                        tosave=self.tosave, sample=self.sample).readfits()

        if useclipped:
            cutout = fits.open(self.savepath + 'unclipped_subcube.fits')[0]

        clipped_cube, _, _, _, sysvel = self.calc_moms()

        #cube_pbcorr = clipped_cube

        # Calculate the beam size, so we can divide by this to get rid of the beam^-1 in the units.
        #psf = self.makebeam(cube_pbcorr.shape[1], cube_pbcorr.shape[2], cube_pbcorr.header)
        #beamsize = np.sum(psf)

        # Make this work if necessary
        #if custom_region:
            #region = pyregion.open(path + 'ds9.reg')
            #mask = region.get_mask(hdu=mom0_hdu)
            #mask_3d = np.tile(mask, (len(cube[:, 0, 0]), 1, 1))
            #cutout = np.where(mask_3d, cube, 0)

        spectrum = np.nansum(cube_pbcorr.data, axis=(1, 2))
        _, _, vel_array_full = self.create_vel_array(cube_pbcorr)

        extra_chans = 5

        if self.galaxy.start - extra_chans < 0:
            spectrum_velocities = vel_array_full[0:self.galaxy.stop + extra_chans]
            spectrum = spectrum[0:self.galaxy.stop + extra_chans]
            spectrum_vel_offset = spectrum_velocities - sysvel + self.galaxy.sysvel_offset
        else:
            spectrum_velocities = vel_array_full[self.galaxy.start - extra_chans:self.galaxy.stop + extra_chans]
            spectrum = spectrum[self.galaxy.start - extra_chans:self.galaxy.stop + extra_chans]
            spectrum_vel_offset = spectrum_velocities - sysvel + self.galaxy.sysvel_offset

        try:
            rest_frequency = cube_pbcorr.header['RESTFRQ']
        except:
            if self.sample == 'viva' or self.sample == 'things':
                print('Warning: setting rest frequency to 1.420 GHz')
                rest_frequency = 1420405751
            else:
                print('Warning: assuming the CO(2-1) line was observed and setting rest frequency to 230.53800000 GHz.')
                rest_frequency = 230538000000

        spectrum_frequencies = rest_frequency * (1 - spectrum_velocities / 299792.458) / 1e9

        if self.tosave:
            clip_rms = self.uncertainty_maps(calc_rms=True)
            if self.sun:
                if self.sample == 'viva' or self.sample == 'things':
                    csv_header = 'Clipping method = Sun+18; rms = ' + str(clip_rms) + ' Jy/b km/s \n \n Spectrum (Jy/b), Velocity (km/s), Velocity offset (km/s), Frequency (GHz)'
                else:
                    csv_header = 'Clipping method = Sun+18; rms = ' + str(clip_rms) + ' K km/s \n \n Spectrum (K), Velocity (km/s), Velocity offset (km/s), Frequency (GHz)'

                np.savetxt(self.savepath + 'spectrum.csv',
                           np.column_stack((spectrum, spectrum_velocities, spectrum_vel_offset, spectrum_frequencies)),
                           delimiter=',', header=csv_header)
            else:
                if self.sample == 'viva' or self.sample == 'things':
                    csv_header = 'Clipping method = Dame11; rms = ' + str(clip_rms) + ' Jy/b km/s \n \n Spectrum (Jy/b), Velocity (km/s), Velocity offset (km/s), Frequency (GHz)'
                else:
                    csv_header = 'Clipping method = Dame11; rms = ' + str(clip_rms) + ' K km/s \n \n Spectrum (K), Velocity (km/s), Velocity offset (km/s), Frequency (GHz)'
                np.savetxt(self.savepath + 'spectrum.csv',
                           np.column_stack((spectrum, spectrum_velocities, spectrum_vel_offset, spectrum_frequencies)),
                           delimiter=',', header=csv_header)

        # Estimate the rms in the spectrum
        #emis, noise = self.splitCube(cutout, self.galaxy.start, self.galaxy.stop)
        #rms = np.std(np.sum(noise, axis=(1, 2)))
        #np.savetxt(path + 'specnoise.txt', [rms / beamsize])

        #psf = self.makebeam(cube_pbcorr.shape[1], cube_pbcorr.shape[2], cube_pbcorr.header)
        #beamsize = np.sum(psf)
        #spectrum /= beamsize

        return spectrum, spectrum_velocities, spectrum_vel_offset, spectrum_frequencies # / beamsize


    def radial_profile(self, alpha_co=5.4, table_path=None, check_aperture=False, hires=False):

        if hires:
            print('WARNING: using a resolution of 1 pixel')

        _, mom0_hdu_Msun, _, _, _ = self.calc_moms(units='M_Sun/pc^2', alpha_co=alpha_co)
        _, mom0_hdu_K, _, _, _ = self.calc_moms(units='K km/s', alpha_co=alpha_co)

        beam_pix = mom0_hdu_K.header['BMAJ'] / mom0_hdu_K.header['CDELT2']

        # Estimate the rms from the spatial inner part of the cube
        cube_pbcorr, cube_uncorr = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                            savepath=self.savepath,
                                            tosave=self.tosave, sample=self.sample).readfits()
        emiscube, noisecube = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                            savepath=self.savepath,
                                            tosave=self.tosave, sample=self.sample).split_cube(cube_pbcorr)
        inner = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                            savepath=self.savepath,
                                            tosave=self.tosave, sample=self.sample).innersquare(noisecube.data)
        rms = np.nanstd(inner)

        if self.sample == 'viva' or self.sample == 'things':
            coldens_atom_cm = rms * 1.10e24 / (cube_pbcorr.header['BMAJ'] * 3600 * cube_pbcorr.header['BMIN'] * 3600)
            Msol_to_matom = 1.187883838e57
            pc_to_cm = 9.521e36
            rms_Msun = coldens_atom_cm / Msol_to_matom * pc_to_cm
        else:
            rms_Msun = rms * alpha_co

        if hires:
            limit = 0 * rms_Msun
        else:
            limit = 2 * rms_Msun

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
            if inc == [90]:
                inc = [89]
            if inc == [0]:
                inc = [1]
            e = np.sin(np.deg2rad(inc))[0]

        centre = (int(mom0_hdu_K.shape[1] / 2), int(mom0_hdu_K.shape[0] / 2))
        hi_inc = False
        rad_prof_K = []
        rad_prof_Msun = []
        radius = []
        area = []

        b_in = -beam_pix + 0.000000000001
        if hires:
            b_in = -1 + 0.0000000001
        b_out = 0
        theta = np.deg2rad(self.galaxy.angle + 90)

        emission_Msun = 2112
        area_temp = 1
        count = 0

        if check_aperture:
            from matplotlib import pyplot as plt
            plt.figure()
            plt.imshow(mom0_hdu_K.data)

        while (emission_Msun / area_temp > limit) | (count < 10):

            # The second condition was added to deal better with galaxies that have a "hole" in the centre.
            count += 1

            if hires:
                b_in += 1
                b_out += 1
            elif count == 1:
                b_in += beam_pix / 2
                b_out += beam_pix / 2
            else:
                b_in += beam_pix
                b_out += beam_pix

            if emission_Msun == 2112:
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

            emission_K = aperture_photometry(mom0_hdu_K.data, aperture)['aperture_sum'][0]
            emission_Msun = aperture_photometry(mom0_hdu_Msun.data, aperture)['aperture_sum'][0]

            area_temp = aperture.area
            area.append(area_temp)
            rad_prof_K.append(emission_K / area_temp)
            rad_prof_Msun.append(emission_Msun / area_temp)
            radius.append(np.average((a_in, a_out)))

        #if ((len(radius) < 5) & (e > 0.7)) or ((len(np.array(rad_prof_K)[np.log10(np.array(rad_prof_K)) < 0]) > 2) & (e > 0.7)):
        if (e > 0.984807753012208):  # Corresponding to an inclination of 80 degrees

            hi_inc = True

            print('HI_INC')

            rad_prof_K = []
            rad_prof_Msun = []
            radius = []
            area = []

            angle = self.galaxy.angle + 180 + 90

            mom0_K_rot = ndimage.interpolation.rotate(mom0_hdu_K.data, angle, reshape=True)
            mom0_Msun_rot = ndimage.interpolation.rotate(mom0_hdu_Msun.data, angle, reshape=True)

            inner = 0
            centre = (mom0_K_rot.shape[1]) / 2
            emission_Msun = 2112
            count = 0

            while (emission_Msun > limit) | (count < 10):

                count += 1

                if hires:
                    slice1_K = mom0_K_rot[:, int(centre + inner):int(centre + inner + 1)]
                    slice2_K = mom0_K_rot[:, int(centre - inner - 1):int(centre - inner)]
                    emission_K = np.average(np.average(slice1_K) + np.average(slice2_K))

                    slice1_Msun = mom0_Msun_rot[:, int(centre + inner):int(centre + inner + 1)]
                    slice2_Msun = mom0_Msun_rot[:, int(centre - inner - 1):int(centre - inner)]
                    emission_Msun = np.average(np.average(slice1_Msun) + np.average(slice2_Msun))
                else:
                    slice1_K = mom0_K_rot[:, int(centre + inner):int(centre + inner + beam_pix)]
                    slice2_K = mom0_K_rot[:, int(centre - inner - beam_pix):int(centre - inner)]
                    emission_K = np.average(np.average(slice1_K) + np.average(slice2_K))

                    slice1_Msun = mom0_Msun_rot[:, int(centre + inner):int(centre + inner + beam_pix)]
                    slice2_Msun = mom0_Msun_rot[:, int(centre - inner - beam_pix):int(centre - inner)]
                    emission_Msun = np.average(np.average(slice1_Msun) + np.average(slice2_Msun))

                if check_aperture:
                    if hires:
                        mom0_K_rot[:, int(centre + inner):int(centre + inner + 1)] = 10
                        mom0_K_rot[:, int(centre - inner - 1):int(centre - inner)] = 20
                    else:
                        mom0_K_rot[:, int(centre + inner):int(centre + inner + beam_pix)] = 10
                        mom0_K_rot[:, int(centre - inner - beam_pix):int(centre - inner)] = 20
                    from matplotlib import pyplot as plt
                    plt.imshow(mom0_K_rot)

                rad_prof_K.append(emission_K)
                rad_prof_Msun.append(emission_Msun)
                area.append(len(slice1_K[slice1_K > 0]) + len(slice2_K[slice2_K > 0]))

                if hires:
                    radius.append(inner + 0.5)
                    inner += 1
                else:
                    radius.append(inner + 0.5 * beam_pix)
                    inner += beam_pix

        rad_prof_K = rad_prof_K[:-1]
        rad_prof_Msun = rad_prof_Msun[:-1]
        radius = np.array(radius[:-1])
        area = area[:-1]
        radii_deg = radius * mom0_hdu_K.header['CDELT2']
        try:
            radii_kpc = np.deg2rad(radii_deg) * self.galaxy.distance * 1000
        except:
            radii_kpc = np.zeros(len(radii_deg))
        N_beams = np.array(area) / (beam_pix ** 2 * np.pi)
        error_K = np.sqrt(N_beams) * rms
        error_Msun = np.sqrt(N_beams) * rms_Msun

        w = wcs.WCS(mom0_hdu_K.header, naxis=2)
        centre_sky = wcs.utils.pixel_to_skycoord(mom0_hdu_K.shape[0] / 2, mom0_hdu_K.shape[1] / 2, wcs=w)

        clip_rms = self.uncertainty_maps(calc_rms=True)

        if check_aperture:
            plt.savefig('/home/nikki/Documents/Data/VERTICO/apertures/' + self.galaxy.name + '.png')

        if hi_inc:
            if self.sun:
                if self.sample == 'viva' or self.sample == 'things':
                    csv_header = 'Slices parallel to the minor axis centered around (RA; Dec) = (' + str(np.round(centre_sky.ra.deg, 2)) + \
                                 '; ' + str(np.round(centre_sky.dec.deg, 2)) + ') (pixel value = (' + str(mom0_hdu_K.shape[0] / 2) + \
                                 '; ' + str(mom0_hdu_K.shape[1] / 2) + ')). ' \
                                'Radii are equal to one beamsize.  \n ' \
                                 'Clipping method = Sun+18; rms = ' + str(clip_rms) + ' Jy/b km/s \n \n' \
                                 'Integrated intensity (Jy/b km/s), RMS error (Jy/b km/s), ' \
                                 'Surface density (M_Sun pc^-2), RMS error (M_Sun pc^-2), Radii (arcsec), Radii (kpc)'
                else:
                    csv_header = 'Slices parallel to the minor axis centered around (RA; Dec) = (' + str(np.round(centre_sky.ra.deg, 2)) + \
                                 '; ' + str(np.round(centre_sky.dec.deg, 2)) + ') (pixel value = (' + str(mom0_hdu_K.shape[0] / 2) + \
                                 '; ' + str(mom0_hdu_K.shape[1] / 2) + ')). ' \
                                'Radii are equal to one beamsize.  \n ' \
                                 'Clipping method = Sun+18; rms = ' + str(clip_rms) + ' K km/s \n \n' \
                                 'Integrated intensity (K km/s), RMS error (K km/s), ' \
                                 'Surface density (M_Sun pc^-2), RMS error (M_Sun pc^-2), Radii (arcsec), Radii (kpc)'
            else:
                if self.sample == 'viva' or self.sample == 'things':
                    csv_header = 'Slices parallel to the minor axis centered around (RA; Dec) = (' + str(np.round(centre_sky.ra.deg, 2)) + \
                                 '; ' + str(np.round(centre_sky.dec.deg, 2)) + ') (pixel value = (' + str(mom0_hdu_K.shape[0] / 2) + \
                                 '; ' + str(mom0_hdu_K.shape[1] / 2) + ')). ' \
                                'Radii are equal to one beamsize. \n' \
                                 'Clipping method = Dame11; rms = ' + str(clip_rms) + ' Jy/b km/s \n \n' \
                                 'Integrated intensity (Jy/b km/s), RMS error (Jy/b km/s), ' \
                                 'Surface density (M_Sun pc^-2), RMS error (M_Sun pc^-2), Radii (arcsec), Radii (kpc)'
                else:
                    csv_header = 'Slices parallel to the minor axis centered around (RA; Dec) = (' + str(np.round(centre_sky.ra.deg, 2)) + \
                                 '; ' + str(np.round(centre_sky.dec.deg, 2)) + ') (pixel value = (' + str(mom0_hdu_K.shape[0] / 2) + \
                                 '; ' + str(mom0_hdu_K.shape[1] / 2) + ')). ' \
                                'Radii are equal to one beamsize. \n' \
                                 'Clipping method = Dame11; rms = ' + str(clip_rms) + ' K km/s \n \n' \
                                 'Integrated intensity (K km/s), RMS error (K km/s), ' \
                                 'Surface density (M_Sun pc^-2), RMS error (M_Sun pc^-2), Radii (arcsec), Radii (kpc)'
        else:
            if self.sun:
                if self.sample == 'viva' or self.sample == 'things':
                    csv_header = 'Elliptical apertures centered around (RA; Dec) = (' + str(np.round(centre_sky.ra.deg, 2)) + \
                                 '; ' + str(np.round(centre_sky.dec.deg, 2)) + ') (pixel value = (' + str(mom0_hdu_K.shape[0] / 2) + \
                                 '; ' + str(mom0_hdu_K.shape[1] / 2) + ')). ' \
                                'Radii are defined as the semi-major axes of these apertures. \n ' \
                                'Clipping method = Sun+18; rms = ' + str(clip_rms) + ' Jy/b km/s \n \n' \
                                'Integrated intensity (Jy/b km/s), RMS error (Jy/b km/s), ' \
                                 'Integrated intensity (M_Sun pc^-2), RMS error (M_Sun pc^-2), Radii (arcsec), Radii (kpc)'
                else:
                    csv_header = 'Elliptical apertures centered around (RA; Dec) = (' + str(np.round(centre_sky.ra.deg, 2)) + \
                                 '; ' + str(np.round(centre_sky.dec.deg, 2)) + ') (pixel value = (' + str(mom0_hdu_K.shape[0] / 2) + \
                                 '; ' + str(mom0_hdu_K.shape[1] / 2) + ')). ' \
                                'Radii are defined as the semi-major axes of these apertures. \n ' \
                                'Clipping method = Sun+18; rms = ' + str(clip_rms) + ' K km/s \n \n' \
                                'Integrated intensity (K km/s), RMS error (K km/s), ' \
                                 'Surface density (M_Sun pc^-2), RMS error (M_Sun pc^-2), Radii (arcsec), Radii (kpc)'
            else:
                if self.sample == 'viva' or self.sample == 'things':
                    csv_header = 'Elliptical apertures centered around (RA; Dec) = (' + str(np.round(centre_sky.ra.deg, 2)) + \
                                 '; ' + str(np.round(centre_sky.dec.deg, 2)) + ') (pixel value = (' + str(mom0_hdu_K.shape[0] / 2) + \
                                 '; ' + str(mom0_hdu_K.shape[1] / 2) + ')). ' \
                                'Radii are defined as the semi-major axes of these apertures. \n' \
                                'Clipping method = Dame11; rms = ' + str(clip_rms) + ' Jy/b km/s \n \n' \
                                'Integrated intensity (Jy/b km/s), RMS error (Jy/b km/s), ' \
                                 'Surface density (M_Sun pc^-2), RMS error (M_Sun pc^-2), Radii (arcsec), Radii (kpc)'
                else:
                    csv_header = 'Elliptical apertures centered around (RA; Dec) = (' + str(np.round(centre_sky.ra.deg, 2)) + \
                                 '; ' + str(np.round(centre_sky.dec.deg, 2)) + ') (pixel value = (' + str(mom0_hdu_K.shape[0] / 2) + \
                                 '; ' + str(mom0_hdu_K.shape[1] / 2) + ')). ' \
                                'Radii are defined as the semi-major axes of these apertures. \n' \
                                'Clipping method = Dame11; rms = ' + str(clip_rms) + ' K km/s \n \n' \
                                'Integrated intensity (K km/s), RMS error (K km/s), ' \
                                 'Surface density (M_Sun pc^-2), RMS error (M_Sun pc^-2), Radii (arcsec), Radii (kpc)'

        if not hires:
            if self.tosave:
                np.savetxt(self.savepath + 'rad_prof.csv',
                           np.column_stack((rad_prof_K, np.ones(len(rad_prof_K)) * error_K, rad_prof_Msun,
                                            np.ones(len(rad_prof_Msun)) * error_Msun, radii_deg * 3600, radii_kpc)),
                                            delimiter=',', header=csv_header)
        else:
            np.savetxt('/home/nikki/Documents/Data/VERTICO/QuenchMechs/radial_profiles/CO/' + 'radial_profile_' + self.galaxy.name + '.csv',
                       np.column_stack((rad_prof_K, np.ones(len(rad_prof_K)) * error_K, rad_prof_Msun,
                                        np.ones(len(rad_prof_Msun)) * error_Msun, radii_deg * 3600, radii_kpc)),
                       delimiter=',', header=csv_header)

        return rad_prof_K, rms, rad_prof_Msun, rms_Msun, radii_deg * 3600, radii_kpc


    def uncertainty_maps(self, calc_rms=False):

        #rmscube = self.calc_noise_in_cube()

        # Read in masks used to clip
        mask = fits.open(self.savepath + 'mask_cube.fits')[0]

        # Map of the number of channels used to calculate the moments
        if self.redo_clip:
            _, mask_trimmed = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                        savepath=self.savepath, tosave=self.tosave, sample=self.sample).\
                do_clip(clip_also=mask, clip_also_nat='mask')
        elif os.path.exists(self.savepath + 'mask_clipped_trimmed.fits'):
            mask_trimmed = fits.open(self.savepath + 'mask_subcube_slab.fits')[0]
        else:
            _, mask_trimmed = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                       savepath=self.savepath, tosave=self.tosave, sample=self.sample). \
                do_clip(clip_also=mask, clip_also_nat='mask')

        mask_trimmed = fits.open(self.savepath + 'mask_subcube_slab.fits')[0]
        N_map = np.sum(mask_trimmed.data, axis=0)

        # Read in the cubes, split into cubes containing line and line-free channels, take the inner cube and calculate
        # the rms from that.
        cube_pbcorr, cube_uncorr = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                            savepath=self.savepath,
                                            tosave=self.tosave, sample=self.sample).readfits()
        #emiscube, noisecube = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
        #                                    savepath=self.savepath,
        #                                    tosave=self.tosave, sample=self.sample).split_cube(cube_uncorr)
        #noisecube = fits.open(self.savepath + 'noise_subcube_slab.fits')[0]
        #inner = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
        #                                    savepath=self.savepath,
        #                                    tosave=self.tosave, sample=self.sample).innersquare(noisecube.data)

        rms = mask.header['CLIP_RMS']

        if calc_rms:
            return rms

        # Read in moment maps
        cube, mom0_hdu, mom1_hdu, mom2_hdu, sysvel = self.calc_moms(units='K km/s')

        # Read in the PB cube and use the middle channel to estimate the noise
        try:
            pb_hdu = fits.open('/home/nikki/Documents/Data/VERTICO/ReducedData/' + str(self.galaxy.name) + '/' +
                           str(self.galaxy.name) + '_7m_co21_pb_rebin.fits')[0]
            if self.redo_clip:
                _, pb_cube = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                                            savepath=self.savepath, tosave=self.tosave, sample=self.sample).\
                    do_clip(clip_also=pb_hdu, clip_also_nat='pb')
            elif os.path.exists(self.savepath + 'pb_subcube_slab.fits'):
                pb_cube = fits.open(self.savepath + 'pb_subcube_slab.fits')[0]
            else:
                _, pb_cube = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun,
                        savepath=self.savepath, tosave=self.tosave, sample=self.sample).\
                    do_clip(clip_also=pb_hdu, clip_also_nat='pb')

            pb_map = pb_cube.data[int(pb_cube.shape[0] / 2), :, :]

        except:
            pb_map = np.ones(mom0_hdu.shape)

        # Total noise map is the rms divided by the PB map
        noise_map = rms / pb_map

        if cube.header['CTYPE3'] == 'FREQ':
            if self.sample == 'viva' or self.sample == 'things':
                v_val = 299792.458 * (1 - (cube.header['CRVAL3'] / 1e9) / 1.420405752)
                v_shift = 299792.458 * (1 - ((cube.header['CRVAL3'] + cube.header['CDELT3']) / 1e9) / 1.420405752)
                v_step = - (v_val - v_shift)
            else:
                v_val = 299792.458 * (1 - (cube.header['CRVAL3'] / 1e9) / 230.538000)
                v_shift = 299792.458 * (1 - ((cube.header['CRVAL3'] + cube.header['CDELT3']) / 1e9) / 230.538000)
                v_step = - (v_val - v_shift)

            mom0_uncertainty = noise_map * np.sqrt(N_map) * abs(v_step)
            mom0_uncertainty = np.where(mom0_hdu.data > 0, mom0_uncertainty, np.nan)
            mom0_uncertainty = fits.PrimaryHDU(mom0_uncertainty, mom0_hdu.header)
            mom1_uncertainty = (N_map * abs(v_step) / (2 * np.sqrt(3))) * \
                               (mom0_uncertainty.data / mom0_hdu.data)  # Eqn 15 doc. Chris.
            mom2_uncertainty = ((N_map * abs(v_step)) ** 2 / (8 * np.sqrt(5))) * \
                               (mom0_uncertainty.data / mom0_hdu.data) * (mom2_hdu.data) ** -1  # Eqn 30 doc. Chris.

        else:
            mom0_uncertainty = noise_map * np.sqrt(N_map) * abs(cube.header['CDELT3'] / 1000)
            mom0_uncertainty = np.where(mom0_hdu.data > 0, mom0_uncertainty, np.nan)
            mom0_uncertainty = fits.PrimaryHDU(mom0_uncertainty, mom0_hdu.header)
            mom1_uncertainty = (N_map * abs(cube.header['CDELT3'] / 1000) / (2 * np.sqrt(3))) * \
                               (mom0_uncertainty.data / mom0_hdu.data)  # Eqn 15 doc. Chris.

            mom2_uncertainty = ((N_map * abs(cube.header['CDELT3'] / 1000)) ** 2 / (8 * np.sqrt(5))) * \
                               (mom0_uncertainty.data / mom0_hdu.data) * (mom2_hdu.data) ** -1  # Eqn 30 doc. Chris.

        mom1_uncertainty = fits.PrimaryHDU(mom1_uncertainty, mom1_hdu.header)
        mom2_uncertainty = fits.PrimaryHDU(mom2_uncertainty, mom2_hdu.header)

        SN = mom0_hdu.data / mom0_uncertainty.data
        SN_hdu = fits.PrimaryHDU(SN, mom0_hdu.header)
        SN_hdu.header.pop('BUNIT')

        self.add_clipping_keywords(mom0_uncertainty.header)
        self.add_clipping_keywords(SN_hdu.header)
        self.add_clipping_keywords(mom1_uncertainty.header)
        self.add_clipping_keywords(mom2_uncertainty.header)

        if self.tosave:
            mom0_uncertainty.writeto(self.savepath + 'mom0_unc.fits', overwrite=True)
            SN_hdu.writeto(self.savepath + 'mom0_SN.fits', overwrite=True)
            mom1_uncertainty.writeto(self.savepath + 'mom1_unc.fits', overwrite=True)
            mom2_uncertainty.writeto(self.savepath + 'mom2_unc.fits', overwrite=True)

        return mom0_uncertainty, SN_hdu, mom1_uncertainty, mom2_uncertainty


    def peak_temperature(self):
        if self.redo_clip:
            cube, _ = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun, savepath=self.savepath,
                            tosave=self.tosave, sample=self.sample).do_clip()
        elif os.path.exists(self.savepath + 'subcube_slab.fits'):
            cube = fits.open(self.savepath + 'subcube_slab.fits')[0]
        else:
            cube, _ = ClipCube(self.galaxy.name, self.path_pbcorr, self.path_uncorr, sun=self.sun, savepath=self.savepath,
                            tosave=self.tosave, sample=self.sample).do_clip()

        peak_temp = np.amax(cube.data, axis=0)

        peak_temp_hdu = fits.PrimaryHDU(peak_temp, self.new_header(cube.header))
        peak_temp_hdu.header['BTYPE'] = 'Peak temperature'
        peak_temp_hdu.header['BUNIT'] = 'K'; peak_temp_hdu.header.comments['BUNIT'] = ''
        self.add_clipping_keywords(peak_temp_hdu.header)

        if self.tosave:
            peak_temp_hdu.writeto(self.savepath + 'peakT.fits', overwrite=True)

        return peak_temp_hdu