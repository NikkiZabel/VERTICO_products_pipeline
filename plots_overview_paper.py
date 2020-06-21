from matplotlib import pyplot as plt
import os
from glob import glob
import numpy as np
from matplotlib.gridspec import GridSpec
from astropy.io import fits
import aplpy as apl
from matplotlib import rcParams
from sauron_colormap import register_sauron_colormap; register_sauron_colormap()
from gal_params import parameters

distance = 16.5  # Mpc


def set_rc_params():

    params = {
        'text.usetex': True,
        'font.family': 'serif',
        'font.serif': 'Times New Roman',
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Times New Roman',
        'mathtext.it': 'Times New Roman:italic',
        'mathtext.bf': 'Times New Roman: bold',
        'mathtext.sf': 'Times New Roman',
        'mathtext.tt': 'Times New Roman',
        'mathtext.cal': 'Times New Roman: italic',
        'font.size': 13,
        'xtick.labelsize': 12,
        'xtick.minor.width': 1,
        'ytick.minor.width': 1,
        'xtick.direction': 'in',
        'ytick.direction': 'in',

    }

    rcParams.update(params)

    return

def make_square(img):
    shape_diff = np.max(img.shape) - np.min(img.shape)
    square_img = np.ones((np.max(img.shape), np.max(img.shape))) * np.nan

    new_header = img.header.copy()

    if img.shape[0] > img.shape[1]:
        new_header['CRPIX1'] += int(shape_diff / 2)
        square_img[:, int(shape_diff / 2):int(shape_diff / 2 + img.shape[1])] = img.data
    else:
        new_header['CRPIX2'] += int(shape_diff / 2)
        square_img[int(shape_diff / 2):int(shape_diff / 2 + img.shape[0]), :] = img.data

    return fits.PrimaryHDU(square_img, new_header)


def image_mom0(image, units='K km/s', peak=False, show_beam=True, show_scalebar=True):

    if peak:
        subplot = list(gs[1, 3].get_position(f).bounds)
        subplot[1] -= 0.02
        subplot[3] += 0.05
        fig = apl.FITSFigure(image, figure=f, subplot=subplot)
        fig.axis_labels.set_yposition('right')
        fig.tick_labels.set_yposition('right')
    else:
        subplot = list(gs[0, 2].get_position(f).bounds)
        subplot[1] -= 0.01
        subplot[3] += 0.05
        fig = apl.FITSFigure(image, figure=f, subplot=subplot)
        fig.axis_labels.hide_x()
        fig.axis_labels.hide_y()
        fig.tick_labels.hide_x()
        fig.tick_labels.hide_y()

    fig.set_theme('publication')
    fig.add_grid()
    fig.grid.set_color('grey')
    fig.grid.set_alpha(0.5)
    fig.grid.set_linestyle('--')

    fig.show_contour(image, cmap='magma_r', levels=np.linspace(np.nanmax(image.data) * 1e-9, np.nanmax(image.data), 20),
                     filled=True, overlap=True)

    # axes and ticks
    fig.ticks.set_color('black')
    #fig.ticks.set_length(10)
    #fig.ticks.set_linewidth(2)
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.tick_labels.set_yformat('dd:mm:ss')
    fig.ticks.set_minor_frequency(5)

    # add a colourbar
    colors = plt.contourf([[0, 0], [0, 0]],
                          levels=np.linspace(0, np.nanmax(image.data) + np.nanmax(image.data) * 0.05, 20),
                          cmap='magma_r')

    if np.nanmax(image.data) < 0.1:
        ticks = np.arange(0, np.nanmax(image.data) + 0.015, 0.005)
    if np.nanmax(image.data) < 0.5:
        ticks = np.arange(0, np.nanmax(image.data) + 0, 0.05)
    elif np.nanmax(image.data) < 1:
        ticks = np.arange(0, np.nanmax(image.data) + 0.1, 0.1)
    elif np.nanmax(image.data) < 2:
        ticks = np.arange(0, np.nanmax(image.data) + 0.2, 0.2)
    elif np.nanmax(image.data) < 5:
        ticks = np.arange(0, np.nanmax(image.data) + 0.5, 0.5)
    elif np.nanmax(image.data) < 10:
        ticks = np.arange(0, np.nanmax(image.data) + 1, 1)
    elif np.nanmax(image.data) < 20:
        ticks = np.arange(0, np.nanmax(image.data) + 1, 2)
    elif np.nanmax(image.data) < 100:
        ticks = np.arange(0, np.nanmax(image.data) + 5, 10)
    elif np.nanmax(image.data) < 200:
        ticks = np.arange(0, np.nanmax(image.data) + 10, 20)
    elif np.nanmax(image.data) < 1000:
        ticks = np.arange(0, np.nanmax(image.data) + 20, 40)
    else:
        ticks = np.arange(0, np.nanmax(image.data) + 100, 200)

    cbar = f.colorbar(colors, ticks=ticks, location='top', pad=0)
    if peak:
        cbar.set_label('Peak temperature [K]')
    elif units == 'K km/s':
        cbar.set_label(r'Integrated intensity [K km s$^{-1}$]')
    elif units == 'M_Sun/pc^2':
        cbar.set_label(r'Surface density [M$_\odot$ pc$^{-2}$]')
    else:
        raise AttributeError('Please choose between "K km/s" and "M_Sun/pc^2"')

    if show_beam:
        fig.add_beam(frame=False)  # automatically imports BMAJ, BMIN, and BPA
        fig.beam.set_edgecolor('k')
        fig.beam.set_color('k')
        fig.beam.set_borderpad(1)

    if show_scalebar:
        length = np.degrees(1.e-3 / distance)  # length of the scalebar in degrees, corresponding to 1 kpc
        fig.add_scalebar(length=length, label='1 kpc', frame=False)
        fig.scalebar.set_linewidth(5)


def image_mom1_2(image, sysvel, moment=1, show_beam=True, show_scalebar=True):

    name, vrange, vrange2, cliplevel, stokes, start, stop, sysvel_offset, angle, \
    full_width, distance, nchan_low, cliplevel_low, nchan_high, cliplevel_high, prune_by_npix, \
    prune_by_fracbeam, expand_by_fracbeam, expand_by_npix, expand_by_nchan, inclination, eccentricity, \
    figsize = parameters('NGC4216')

    vel_array = np.genfromtxt('/home/nikki/Documents/Data/VERTICO/Products/NGC4216/Smooth_method/NGC4216_7m+tp_co21_pbcorr_round_k_spectrum.csv',
                               delimiter=',', skip_header=1)[:, 1][start:stop]

    sysvel = (sysvel + 5) // 10 * 10

    # show the image in colour
    if moment == 1:
        subplot = list(gs[0, 3].get_position(f).bounds)
        subplot[1] -= 0.01
        subplot[3] += 0.05
        fig = apl.FITSFigure(image, figure=f, subplot=subplot)
        fig.axis_labels.hide_x()
        fig.tick_labels.hide_x()
        fig.axis_labels.set_yposition('right')
        fig.tick_labels.set_yposition('right')
    elif moment == 2:
        subplot = list(gs[1, 2].get_position(f).bounds)
        subplot[1] -= 0.02
        subplot[3] += 0.05
        fig = apl.FITSFigure(image, figure=f, subplot=subplot)
        fig.axis_labels.hide_y()
        fig.tick_labels.hide_y()
    else:
        print('No.')

    # axes and ticks
    fig.ticks.set_color('black')
    #fig.ticks.set_length(10)
    #fig.ticks.set_linewidth(2)
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.tick_labels.set_yformat('dd:mm:ss')
    fig.ticks.show()
    fig.ticks.set_minor_frequency(5)
    fig.add_grid()
    fig.grid.set_color('grey')
    fig.grid.set_alpha(0.5)
    fig.grid.set_linestyle('--')

    #add a colourbar
    if moment == 2:

        if not vrange2:
            vrange2 = np.nanmax(image.data - sysvel)

        fig.show_contour(image, cmap='sauron', levels=np.linspace(0, vrange2, len(vel_array)), vmin=0,
                         vmax=vrange2, extend='both', filled=True, overlap=True)
        colors = plt.contourf([[0, 0], [0, 0]], levels=np.linspace(0, vrange2, len(vel_array)),
                              cmap='sauron')

        if vrange2 < 11:
            ticks = np.arange(0, vrange2 + 1, 1)
        elif vrange2 < 100:
            ticks = np.arange(0, vrange2 + 10, 10)
        else:
            ticks = np.arange(0, vrange2 + 20, 20)
        cbar = f.colorbar(colors, ticks=ticks, location='top', pad=0)
        cbar.set_label(r'Observed $\sigma_v$ [km s$^{-1}$]')

    else:
        if not vrange:
            vrange = 1.5 * (np.nanmax(image.data) - sysvel)

        fig.show_contour(image, cmap='sauron', levels=np.linspace(-vrange, vrange,
            len(vel_array)), vmin=-vrange, vmax=vrange, extend='both', filled=True,
                         overlap=True)
        colors = plt.contourf([[0, 0], [0, 0]], levels=np.linspace(-vrange, vrange,
                                                                   len(vel_array)), cmap='sauron')
        if vrange < 16:
            tickarr = np.arange(-vrange, 0, 3)
        elif vrange < 60:
            tickarr = np.arange(-vrange, 0, 15)
        elif vrange < 130:
            tickarr = np.arange(-vrange, 0, 30)
        else:
            tickarr = np.arange(-vrange, 0, 60)

        ticks = np.concatenate((tickarr, [0], abs(tickarr)))
        cbar = f.colorbar(colors, ticks=ticks, location='top', pad=0)
        cbar.set_label(r'Velocity [km s$^{-1}$]')
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation=45)

    if show_beam:
        fig.add_beam(frame=False)  # automatically imports BMAJ, BMIN, and BPA
        fig.beam.set_edgecolor('k')
        fig.beam.set_color('k')
        fig.beam.set_borderpad(1)

    if show_scalebar:
        length = np.degrees(1.e-3 / distance)  # length of the scalebar in degrees, corresponding to 1 kpc
        fig.add_scalebar(length=length, label='1 kpc', frame=False)
        fig.scalebar.set_linewidth(5)

    return


def contour_plot(image, contours, number=10):

    # Set beam parameters to the ones of the CO image, so that the beam can be plotted
    image.header['BMAJ'] = contours.header['BMAJ']
    image.header['BMIN'] = contours.header['BMIN']
    image.header['BPA'] = contours.header['BPA']

    fig = apl.FITSFigure(image, figure=f, subplot=list(gs[0:2, 0:2].get_position(f).bounds))
    fig.show_grayscale()
    fig.show_contour(contours, levels=np.linspace(0, np.nanmax(contours.data), number),
                     cmap='winter_r')

    # Show the beam
    fig.add_beam(frame=False)
    fig.beam.set_edgecolor('k')
    fig.beam.set_color('white')
    fig.beam.set_borderpad(1)

    # Show a scalebar
    length = np.degrees(1.e-3 / distance)
    fig.add_scalebar(length=length, label='1 kpc', frame=False)
    fig.scalebar.set_linewidth(2)
    fig.scalebar.set_color('white')

    # Add galaxy name in the upper left corner
    fig.add_label(0.15, 0.9, 'NGC4216', color='white', relative=True, size=20)

    # Make some adjustments
    plt.gca().tick_params(which='both', length=10, width=1.5)
    plt.gca().tick_params(which='minor', length=5)
    fig.tick_labels.set_font(size=17)
    fig.axis_labels.set_font(size=20)

    return


# Read in the desired data products CURRENTLY DAME VERSION
path = '/home/nikki/Documents/Data/VERTICO/Products/NGC4216/Smooth_method/'
mom0 = fits.open(path + 'NGC4216_7m+tp_co21_pbcorr_round_k_moment0_K.fits')[0]
mom1 = fits.open(path + 'NGC4216_7m+tp_co21_pbcorr_round_k_moment1.fits')[0]
mom2 = fits.open(path + 'NGC4216_7m+tp_co21_pbcorr_round_k_moment2.fits')[0]
peak = fits.open(path + 'NGC4216_7m+tp_co21_pbcorr_round_k_peak_temperature.fits')[0]

# Set zeros to nans
mom0.data[mom0.data == 0] = np.nan
mom1.data[mom1.data == 0] = np.nan
mom2.data[mom2.data == 0] = np.nan
peak.data[peak.data == 0] = np.nan

# Make those images square
mom0 = make_square(mom0)
mom1 = make_square(mom1)
mom2 = make_square(mom2)
peak = make_square(peak)

# Read in SDSS data
g_band = fits.open('/home/nikki/Documents/Data/VERTICO/SDSS/NGC4216_g.fits')[0]

# Prepare the figure layout
set_rc_params()
f = plt.figure(figsize=(15, 7.2))
gs = GridSpec(2, 4, figure=f)
ax_sdss = f.add_subplot(gs[0:2, 0:2])
ax_mom0 = f.add_subplot(gs[0, 2])
ax_mom1 = f.add_subplot(gs[0, 3])
ax_mom2 = f.add_subplot(gs[1, 2])
ax_peak = f.add_subplot(gs[1, 3])
ax_sdss.axis('off')
ax_mom0.axis('off')
ax_mom1.axis('off')
ax_mom2.axis('off')
ax_peak.axis('off')

# Fill the panels with the plots (position on the gridspec is currently hardcoded)
contour_plot(g_band, mom0, number=8)
image_mom0(mom0, units='K km/s', show_beam=False, show_scalebar=False)
image_mom1_2(mom1, mom1.header['SYSVEL'], moment=1, show_beam=False, show_scalebar=False)
image_mom1_2(mom2, mom1.header['SYSVEL'], moment=2, show_beam=False, show_scalebar=False)
image_mom0(peak, peak=True, show_beam=False, show_scalebar=False)

plt.savefig('/home/nikki/test.pdf', bbox_inches='tight')