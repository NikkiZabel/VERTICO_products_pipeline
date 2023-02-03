from image_moments import CreateImages

refresh = True
overwrite = True
sun = False
tosave = True

sample = None

path = '/media/nikki/Elements/NGC1436_kana/'
savepath = '/home/nikki/Documents/ForAle/'

galaxy = 'NGC1436'

file_pbcorr = path + 'NGC1436_with_kana.co.image.pbcor.fits'
file_uncorr = path + 'NGC1436.co_robust=2.image.fits'

# Moment maps
CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
             sun=sun, tosave=tosave, sample=sample).moment_zero(units='K km/s')
CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
             sun=sun, tosave=tosave, sample=sample).moment_zero()
CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
             sun=sun, tosave=tosave, sample=sample).moment_zero(peak=True)
CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
             sun=sun, tosave=tosave, sample=sample).moment_1_2()
CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
             sun=sun, tosave=tosave, sample=sample).moment_1_2(moment=2)

# Uncertainty maps
CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
             sun=sun, tosave=tosave).mom0_noise_maps()
CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
             sun=sun, tosave=tosave).mom1_2_noise_maps()

# PVDs
CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
             sun=sun, tosave=tosave).PVD(axis='major', find_angle=False, check_slit=False)
CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
             sun=sun, tosave=tosave).PVD(axis='minor', check_slit=False)

# Spectra
CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
              sun=sun, tosave=tosave).spectrum(x_axis='vel_offset', useclipped=True)
CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
             sun=sun, tosave=tosave).spectrum(x_axis='velocity', useclipped=True)
CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
             sun=sun, tosave=tosave).spectrum(x_axis='frequency', useclipped=True)

# Radial profiles
CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
            sun=sun, tosave=tosave, sample=sample).radial_profile(x_units='arcsec', y_units='M_Sun pc^-2',
                                alpha_co=6.25, table_path=None, check_aperture=False)
CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
            sun=sun, tosave=tosave, sample=sample).radial_profile(x_units='kpc', y_units='M_Sun pc^-2',
                               alpha_co=6.25, table_path=None, check_aperture=False)
CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
            sun=sun, tosave=tosave, sample=sample).radial_profile(y_units='K km/s', x_units='kpc',
                               alpha_co=6.25, table_path=None, check_aperture=False)
CreateImages(galaxy, file_pbcorr, file_uncorr, savepath=savepath, refresh=refresh, overwrite=overwrite,
             sun=sun, tosave=tosave, sample=sample).radial_profile(y_units='K km/s', x_units='arcsec',
                                                                   alpha_co=6.25, table_path=None, check_aperture=False)
