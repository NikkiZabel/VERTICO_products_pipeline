def parameters(galaxy):

    if galaxy == 'NGC4713':
        centre_y = 25
        centre_x = 32
        size = 70
        vrange = 80
        vrange_2 = 50
        cliplevel = 4
        stokes = False
        start = 24
        stop = 56
        angle = 0
        centre_pvd = 25
        sysvel = 25
        distance = 16.5
        nchan_low = None
        cliplevel_low = None
        nchan_high = None
        cliplevel_high = None
        prune_by_npix = None
        prune_by_fracbeam = None
        expand_by_fracbeam = None
        expand_by_npix = None
        expand_by_nchan = None

    elif galaxy == 'NGC4568':
        centre_y = 23
        centre_x = 27
        size = 60
        vrange = 200
        vrange_2 = 75
        cliplevel = 4
        stokes = False
        start = 13
        stop = 63
        angle = 110
        centre_pvd = 32
        sysvel = 19
        distance = 16.5
        nchan_low = None
        cliplevel_low = None
        nchan_high = None
        cliplevel_high = None
        prune_by_npix = None
        prune_by_fracbeam = None
        expand_by_fracbeam = None
        expand_by_npix = None
        expand_by_nchan = None

    elif galaxy == 'NGC4189':
        centre_y = 32
        centre_x = 20
        size = 60
        vrange = 200
        vrange_2 = 120
        cliplevel = 6
        stokes = False
        start = 19
        stop = 59
        angle = 335
        centre_pvd = 33
        sysvel = 14
        distance = 16.5
        nchan_low = 2
        cliplevel_low = 3
        nchan_high = 3
        cliplevel_high = 5
        prune_by_npix = None
        prune_by_fracbeam = None
        expand_by_fracbeam = None
        expand_by_npix = None
        expand_by_nchan = None

    name = galaxy

    return name, centre_x, centre_y, size, vrange, vrange_2, cliplevel, stokes, start, stop, angle, \
           centre_pvd, sysvel, distance, nchan_low, cliplevel_low, nchan_high, cliplevel_high, prune_by_npix, \
           prune_by_fracbeam, expand_by_fracbeam, expand_by_npix, expand_by_nchan
