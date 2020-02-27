import gal_params

class galaxies:
    def __init__(self, galname):

        name, centre_x, centre_y, size, vrange, vrange_2, cliplevel, stokes, start, stop, sysvel_offset, angle, \
        pv_corr, full_width, distance, nchan_low, cliplevel_low, nchan_high, cliplevel_high, prune_by_npix, \
        prune_by_fracbeam, expand_by_fracbeam, expand_by_npix, expand_by_nchan, inclination, eccentricity, high_inc, \
        figsize, rad_prof_corr = gal_params.parameters(galname)

        self.name = name
        self.centre_x = centre_x
        self.centre_y = centre_y
        self.size = size
        self.vrange = vrange
        self.vrange2 = vrange_2
        self.cliplevel = cliplevel
        self.stokes = stokes
        self.start = start
        self.stop = stop
        self.sysvel_offset = sysvel_offset
        self.angle = angle
        self.pv_corr = pv_corr
        self.full_width = full_width
        self.distance = distance
        self.nchan_low = nchan_low
        self.cliplevel_low = cliplevel_low
        self.nchan_high = nchan_high
        self.cliplevel_high = cliplevel_high
        self.prune_by_npix = prune_by_npix
        self.prune_by_fracbeam = prune_by_fracbeam
        self.expand_by_fracbeam = expand_by_fracbeam
        self.expand_by_npix = expand_by_npix
        self.expand_by_nchan = expand_by_nchan
        self.inclination = inclination
        self.eccentricity = eccentricity
        self.high_inc = high_inc
        self.figsize = figsize
        self.rad_prof_corr = rad_prof_corr
