def parameters(galaxy):

    if galaxy == 'NGC4713':
        centre_y = 26
        centre_x = 25
        size = 70
        vrange = 100
        vrange_2 = 40
        cliplevel = 4
        stokes = False
        start = 23
        stop = 57
        sysvel_offset = -20
        angle = 176
        vel_centre_pvd = 25
        distance = 16.5
        nchan_low = 2
        cliplevel_low = 2
        nchan_high = 3
        cliplevel_high = 3.5
        prune_by_npix = None
        prune_by_fracbeam = 1
        expand_by_fracbeam = None
        expand_by_npix = None
        expand_by_nchan = 2
        inclination = None
        eccentricity = 1

    elif galaxy == 'NGC4568':
        centre_y = 23
        centre_x = 27
        size = 60
        vrange = 200
        vrange_2 = 75
        cliplevel = 4
        stokes = False
        start = 15
        stop = 66
        sysvel_offset = 70
        angle = 110 + 180
        vel_centre_pvd = 19
        distance = 16.5
        nchan_low = 2
        cliplevel_low = 2
        nchan_high = 3
        cliplevel_high = 3.5
        prune_by_npix = None
        prune_by_fracbeam = 1
        expand_by_fracbeam = None
        expand_by_npix = None
        expand_by_nchan = 2
        inclination = 70
        eccentricity = None

    elif galaxy == 'NGC4189':
        centre_y = 32
        centre_x = 20
        size = 60
        vrange = 200
        vrange_2 = 50
        cliplevel = 1
        stokes = False
        start = 22
        stop = 61
        sysvel_offset = 50
        angle = 335 - 180
        vel_centre_pvd = 14
        distance = 16.5
        nchan_low = 2
        cliplevel_low = 2
        nchan_high = 3
        cliplevel_high = 3.5
        prune_by_npix = None
        prune_by_fracbeam = 1
        expand_by_fracbeam = None
        expand_by_npix = None
        expand_by_nchan = 2
        inclination = None
        eccentricity = 1

    elif galaxy == 'NGC4299':
        centre_y = 32
        centre_x = 20
        size = 60
        vrange = 200
        vrange_2 = 50
        cliplevel = 1
        stokes = False
        start = 22
        stop = 61
        sysvel_offset = 50
        angle = 335 - 180
        vel_centre_pvd = 14
        distance = 16.5
        nchan_low = 2
        cliplevel_low = 2
        nchan_high = 3
        cliplevel_high = 3.5
        prune_by_npix = None
        prune_by_fracbeam = 1
        expand_by_fracbeam = None
        expand_by_npix = None
        expand_by_nchan = 2
        inclination = None
        eccentricity = 1

    name = galaxy

    return name, centre_x, centre_y, size, vrange, vrange_2, cliplevel, stokes, start, stop, sysvel_offset, angle, \
           vel_centre_pvd, distance, nchan_low, cliplevel_low, nchan_high, cliplevel_high, prune_by_npix, \
           prune_by_fracbeam, expand_by_fracbeam, expand_by_npix, expand_by_nchan, inclination, eccentricity
