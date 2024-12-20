def parameters(galaxy, sample=None):

    if sample == 'viva':
        try:
            galnum = galaxy.split('c')[1]
            if galnum == '3392':
                galaxy = 'IC' + galnum
            else:
                galaxy = 'NGC' + galnum
        except:
            pass

    vrange = None
    vrange_2 = None
    cliplevel = 1.5
    stokes = False
    sysvel_offset = 0
    angle = 0
    full_width = False
    distance = 0
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
    eccentricity = None
    figsize = (10, 8.5)

    if galaxy == 'NGC1436':
        vrange = 120
        vrange_2 = 20
        cliplevel = 1.5
        stokes = False
        sysvel_offset = 0
        angle = 345
        full_width = True
        distance = 19.95
        nchan_low = 1
        cliplevel_low = 1
        nchan_high = 2
        cliplevel_high = 3
        prune_by_npix = None
        prune_by_fracbeam = 1
        expand_by_fracbeam = None
        expand_by_npix = None
        expand_by_nchan = 2
        inclination = 47
        eccentricity = None

    elif galaxy == 'NGC4501':
        vrange = 380
        vrange_2 = 100
        cliplevel = 1.5
        stokes = False
        sysvel_offset = 0
        angle = 165
        full_width = True
        distance = 16.5
        nchan_low = 1
        cliplevel_low = 2
        nchan_high = 1
        cliplevel_high = 2.5
        prune_by_npix = None
        prune_by_fracbeam = 1
        expand_by_fracbeam = None
        expand_by_npix = None
        expand_by_nchan = 2
        inclination = None
        eccentricity = None
        figsize = (7, 10)

    elif galaxy == 'NGC4698':
        vrange = 250
        vrange_2 = 20
        cliplevel = 1.5
        stokes = False
        sysvel_offset = 0
        angle = 165
        full_width = True
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
        eccentricity = None
        figsize = (7, 10)

    elif galaxy == 'NGC4694':
        vrange = 50
        vrange_2 = 50
        cliplevel = 1.5
        stokes = False
        sysvel_offset = 0
        angle = 40
        full_width = False
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
        eccentricity = None
        figsize = (8.5, 8)

    elif galaxy == 'NGC4689':
        vrange = 100
        vrange_2 = 30
        cliplevel = 1.5
        stokes = False
        sysvel_offset = 0
        angle = 155
        full_width = False
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
        eccentricity = None
        figsize = (10, 9.5)

    elif galaxy == 'NGC4654':
        vrange = 150
        vrange_2 = 30
        cliplevel = 1.5
        stokes = False
        sysvel_offset = 0
        angle = 122
        full_width = False
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
        eccentricity = None
        figsize = (10, 6.5)

    elif galaxy == 'NGC4579':
        if sample != 'heracles':
            vrange = 200
            vrange_2 = 120
            cliplevel = 1.5
            stokes = False
            sysvel_offset = 0
            angle = 90
            full_width = False
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
            eccentricity = None
            figsize = (10, 7)

    elif galaxy == 'NGC4569':
        if sample != 'heracles':
            vrange = 200
            vrange_2 = 80
            cliplevel = 1.5
            stokes = False
            sysvel_offset = 0
            angle = 15
            full_width = False
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
            eccentricity = None
            figsize = (9.5, 10)

    elif galaxy == 'NGC4548':
        vrange = 150
        vrange_2 = 80
        cliplevel = 1.5
        stokes = False
        sysvel_offset = 0
        angle = 140
        full_width = False
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
        eccentricity = None
        figsize = (10, 7.5)

    elif galaxy == 'NGC4536':
        vrange = 200
        vrange_2 = 100
        cliplevel = 1.5
        stokes = False
        sysvel_offset = 0
        angle = 290
        full_width = False
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
        eccentricity = None
        figsize = (10, 7)

    elif galaxy == 'NGC4535':
        vrange = 150
        vrange_2 = 50
        cliplevel = 1.5
        stokes = False
        sysvel_offset = 0
        angle = 175
        full_width = False
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
        eccentricity = None
        figsize = (10, 10)

    elif galaxy == 'NGC4457':
        vrange = 120
        vrange_2 = 60
        cliplevel = 1.5
        stokes = False
        sysvel_offset = 0
        angle = 75
        full_width = False
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
        eccentricity = None
        figsize = (10, 7)

    elif galaxy == 'NGC4424':
        vrange = 60
        vrange_2 = 35
        cliplevel = 1.5
        stokes = False
        sysvel_offset = 0
        angle = 95
        full_width = False
        distance = 16.5
        nchan_low = 1
        cliplevel_low = 1
        nchan_high = 2
        cliplevel_high = 3
        prune_by_npix = None
        prune_by_fracbeam = 1
        expand_by_fracbeam = None
        expand_by_npix = None
        expand_by_nchan = 2
        inclination = None
        eccentricity = None
        figsize = (10, 7)

    elif galaxy == 'NGC4402':
        vrange = 150
        vrange_2 = 35
        cliplevel = 1.5
        stokes = False
        sysvel_offset = 0
        angle = 87
        full_width = False
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
        eccentricity = None
        figsize = (15, 5.1)

    elif galaxy == 'NGC4321':
        if sample != 'heracles':
            vrange = 150
            vrange_2 = 60
            cliplevel = 1.5
            stokes = False
            sysvel_offset = 0
            angle = 165
            full_width = False
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
            eccentricity = None
            figsize = (10, 7)

    elif galaxy == 'NGC4298':
        vrange = 150
        vrange_2 = 40
        cliplevel = 1.5
        stokes = False
        sysvel_offset = 0
        angle = 315
        full_width = False
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
        eccentricity = None
        figsize = (10, 8.5)

    elif galaxy == 'NGC4293':
        vrange = 120
        vrange_2 = 70
        cliplevel = 1.5
        stokes = False
        sysvel_offset = 0
        angle = 75
        full_width = False
        distance = 16.5
        nchan_low = 2
        cliplevel_low = 2.5
        nchan_high = 3
        cliplevel_high = 3.5
        prune_by_npix = None
        prune_by_fracbeam = 1
        expand_by_fracbeam = None
        expand_by_npix = None
        expand_by_nchan = 2
        inclination = None
        eccentricity = None
        figsize = (10, 5.5)

    elif galaxy == 'NGC4254':
        if sample != 'heracles':
            vrange = 150
            vrange_2 = 40
            cliplevel = 1.5
            stokes = False
            sysvel_offset = 0
            angle = 50
            full_width = False
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
            eccentricity = None
            figsize = (9, 7.5)

    elif galaxy == 'NGC4561':
        vrange = 30
        vrange_2 = 20
        cliplevel = 1.5
        stokes = False
        sysvel_offset = 0
        angle = 0
        full_width = False
        distance = 16.5
        nchan_low = 1
        cliplevel_low = 2
        nchan_high = 2
        cliplevel_high = 2.5
        prune_by_npix = None
        prune_by_fracbeam = 1
        expand_by_fracbeam = None
        expand_by_npix = None
        expand_by_nchan = 2
        inclination = None
        eccentricity = None
        figsize = (7.5, 9)

    elif galaxy == 'NGC4450':
        vrange = 150
        vrange_2 = 75
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 3
        full_width = False
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
        eccentricity = None
        figsize = (9, 10)

    elif galaxy == 'NGC4580':
        vrange = 100
        vrange_2 = 35
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 150
        full_width = False
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
        eccentricity = None
        figsize = (10, 10)

    elif galaxy == 'NGC4772':
        vrange = 400
        vrange_2 = 30
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 130
        full_width = True
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
        eccentricity = None
        figsize = (10, 9.5)

    elif galaxy == 'NGC4567':
        vrange = 120
        vrange_2 = 50
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 100
        full_width = False
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
        eccentricity = None
        figsize = (10, 6)

    elif galaxy == 'NGC4396':
        vrange = 100
        vrange_2 = 30
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 125
        full_width = False
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
        eccentricity = None
        figsize = (10, 7.5)

    elif galaxy == 'NGC4808':
        vrange = 150
        vrange_2 = 40
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 310
        full_width = False
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
        eccentricity = None
        figsize = (10, 7.5)

    elif galaxy == 'NGC4651':
        vrange = 200
        vrange_2 = 75
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 260
        full_width = False
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
        eccentricity = None
        figsize = (10, 6.5)

    elif galaxy == 'NGC4607':
        vrange = 120
        vrange_2 = 30
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 182
        full_width = False
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
        eccentricity = None
        figsize = (8, 10)

    elif galaxy == 'NGC4606':
        vrange = 100
        vrange_2 = 50
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 215
        full_width = False
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
        eccentricity = None
        figsize = (10, 7.5)

    elif galaxy == 'NGC4533':
        vrange = 100
        vrange_2 = 25
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 165
        full_width = False
        distance = 16.5
        nchan_low = 2
        cliplevel_low = 2.5
        nchan_high = 3
        cliplevel_high = 3.5
        prune_by_npix = None
        prune_by_fracbeam = 1
        expand_by_fracbeam = None
        expand_by_npix = None
        expand_by_nchan = 2
        inclination = None
        eccentricity = None
        figsize = (10, 10)

    elif galaxy == 'NGC4532':
        vrange = 100
        vrange_2 = 30
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 330
        full_width = False
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
        eccentricity = None
        figsize = (9, 10)

    elif galaxy == 'NGC4522':
        vrange = 100
        vrange_2 = 30
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 215
        full_width = False
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
        eccentricity = None
        figsize = (9.5, 10)

    elif galaxy == 'NGC4419':
        vrange = 200
        vrange_2 = 75
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 310
        full_width = False
        distance = 16.5
        nchan_low = 2
        cliplevel_low = 2.5
        nchan_high = 3
        cliplevel_high = 3.5
        prune_by_npix = None
        prune_by_fracbeam = 1
        expand_by_fracbeam = None
        expand_by_npix = None
        expand_by_nchan = 2
        inclination = None
        eccentricity = None
        figsize = (10, 7)

    elif galaxy == 'NGC4405':
        vrange = 100
        vrange_2 = 30
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 200
        full_width = 0
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
        eccentricity = None
        figsize = (10, 9)

    elif galaxy == 'NGC4394':
        vrange = 100
        vrange_2 = 50
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 70
        full_width = True
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
        eccentricity = None
        figsize = (10, 10)

    elif galaxy == 'NGC4388':
        vrange = 200
        vrange_2 = 100
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 90
        full_width = False
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
        eccentricity = None
        figsize = (14.5, 6)

    elif galaxy == 'NGC4383':
        vrange = 100
        vrange_2 = 50
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 200
        full_width = False
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
        eccentricity = None
        figsize = (9, 10)

    elif galaxy == 'NGC4380':
        vrange = 200
        vrange_2 = 50
        cliplevel = 3
        stokes = False
        sysvel_offset = 20
        angle = 333
        full_width = False
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
        eccentricity = None
        figsize = (9, 10)

    elif galaxy == 'NGC4351':
        vrange = 50
        vrange_2 = 30
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 65
        full_width = False
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
        eccentricity = None
        figsize = (10, 6.5)

    elif galaxy == 'NGC4330':
        vrange = 150
        vrange_2 = 30
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 58
        full_width = False
        distance = 16.5
        nchan_low = 2
        cliplevel_low = 2.5
        nchan_high = 3
        cliplevel_high = 3.5
        prune_by_npix = None
        prune_by_fracbeam = 1
        expand_by_fracbeam = None
        expand_by_npix = None
        expand_by_nchan = 2
        inclination = None
        eccentricity = None
        figsize = (10, 6.5)

    elif galaxy == 'NGC4302':
        vrange = 200
        vrange_2 = 75
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 180
        full_width = False
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
        eccentricity = None
        figsize = (6, 20)

    elif galaxy == 'NGC4294':
        vrange = 100
        vrange_2 = 40
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 340
        full_width = False
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
        eccentricity = None
        figsize = (9, 10)

    elif galaxy == 'NGC4222':
        vrange = 100
        vrange_2 = 25
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 58
        full_width = False
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
        eccentricity = None
        figsize = (10, 6.5)

    elif galaxy == 'NGC4216':
        vrange = 300
        vrange_2 = 75
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 200
        full_width = True
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
        eccentricity = None
        figsize = (6.5, 10)

    elif galaxy == 'NGC4192':
        vrange = 250
        vrange_2 = 100
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 150
        full_width = False
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
        eccentricity = None
        figsize = (8.5, 10)

    elif galaxy == 'NGC4064':
        vrange = 120
        vrange_2 = 60
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 345
        full_width = False
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
        eccentricity = None
        figsize = (10, 10)

    elif galaxy == 'IC3392':
        vrange = 120
        vrange_2 = 40
        cliplevel = 3
        stokes = False
        sysvel_offset = 10
        angle = 40
        full_width = False
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
        eccentricity = None
        figsize = (10, 8.5)

    elif galaxy == 'NGC4713':
        vrange = 110
        vrange_2 = 40
        cliplevel = 4
        stokes = False
        sysvel_offset = 0
        angle = 274
        full_width = False
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
        eccentricity = None
        figsize = (10, 7)

    elif galaxy == 'NGC4568':
        vrange = 150
        vrange_2 = 75
        cliplevel = 4
        stokes = False
        sysvel_offset = 0
        angle = 20
        full_width = False
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
        eccentricity = None
        figsize = (11.5, 10)

    elif galaxy == 'NGC4189':
        vrange = 150
        vrange_2 = 50
        cliplevel = 3
        stokes = False
        sysvel_offset = 0
        angle = 245
        full_width = False
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
        eccentricity = None
        figsize = (10, 8)

    elif galaxy == 'NGC4299':
        vrange = 100
        vrange_2 = 25
        cliplevel = 1
        stokes = False
        sysvel_offset = 0
        angle = 270
        full_width = True
        distance = 16.5
        nchan_low = 2
        cliplevel_low = 1
        nchan_high = 3
        cliplevel_high = 3.5
        prune_by_npix = None
        prune_by_fracbeam = 1
        expand_by_fracbeam = None
        expand_by_npix = None
        expand_by_nchan = 2
        inclination = None
        eccentricity = None
        figsize = (10, 8)

    if sample == 'heracles' or sample == 'things':
        import numpy as np
        
        if sample == 'heracles':
            if galaxy == 'NGC2841':
                cliplevel_low = 1
                cliplevel_high = 1.5

        def get_inc_pa(galaxy):
            from astropy.io import fits
            table = fits.open('/media/nikki/6719190437AE6F5C/Work/Data/VERTICO/heracles/heracles_basic.fits')[1]
            gal_name_table = table.data['Galaxy']

            try:
                incl = table.data['inclination'][gal_name_table == galaxy]
                pa = table.data['pa'][gal_name_table == galaxy]
                return incl[0], pa[0]
            except:
                gal_num_table = np.array([n.split('C')[1] for n in gal_name_table])
                galaxy_num = galaxy.split('c')[1]
                incl = table.data['inclination'][gal_num_table == galaxy_num]
                pa = table.data['pa'][gal_num_table == galaxy_num]
                return incl[0], pa[0]

        def get_dist(galaxy):
            from astropy.io import fits
            dist_table = fits.open('/media/nikki/6719190437AE6F5C/Work/Data/VERTICO/heracles/heracles_basic.fits')[1]
            try:
                return float(dist_table.data['dist_L19'][dist_table.data['Galaxy'] == galaxy])
            except:
                from astropy.coordinates import Distance
                return Distance(z=dist_table.data['Redshift'][dist_table.data['Galaxy'] == galaxy]).value
            
        inclination, pa = get_inc_pa(galaxy)
        distance = get_dist(galaxy)

        angle = pa

    if sample == 'viva':
        if galaxy == 'NGC4580' or galaxy == 'NGC4654':
            nchan_low = 1
            nchan_high = 2
        elif galaxy == 'NGC4694' or galaxy == 'NGC4606' or galaxy == 'NGC4533' \
            or galaxy == 'NGC4254' or galaxy == 'NGC4561' or galaxy == 'NGC4579' or galaxy == 'NGC4192':
            nchan_low = 1
            nchan_high = 2
            cliplevel_low = 1
            cliplevel_high = 2
        elif galaxy == 'NGC4450' or galaxy == 'NGC4501' or galaxy == 'NGC4569':
            nchan_low = 1
            nchan_high = 2
            cliplevel_low = 1
            cliplevel_high = 2.5
        elif galaxy == 'NGC4713' or galaxy == 'NGC4772' or galaxy == 'NGC4808' \
            or galaxy == 'NGC4222' or galaxy == 'NGC4298' or galaxy == 'NGC4380' \
            or galaxy == 'NGC4522':
            nchan_low = 1
            nchan_high = 3
            cliplevel_low = 1
            cliplevel_high = 2.5
        elif galaxy == 'NGC4189' or galaxy == 'NGC4330':
            cliplevel_low = 1
            cliplevel_high = 2
        elif galaxy == 'NGC4351' or galaxy == 'NGC4419':
            cliplevel_low = 1.5
            cliplevel_high = 2
        elif galaxy == 'NGC4548':
            nchan_low = 2
            nchan_high = 3
            cliplevel_low = 1
            cliplevel_high = 2
            
    name = galaxy
    figsize = (10, 8.5)

    if sample == 'things':
        stokes = False

    # Overwrite PA and inclination with new ones calculated by Toby
    if sample == 'vertico' or sample == 'viva':
        from astropy.io import fits
        import numpy as np
        table = fits.open('/media/nikki/6719190437AE6F5C/Work/Data/VERTICO/VERTICO_basic.fits')[1]
        inclination = int(np.round(table.data['inclination'][table.data['Galaxy'] == galaxy][0]))
        angle = int(np.round(table.data['pa'][table.data['Galaxy'] == galaxy][0]))

    return name, vrange, vrange_2, cliplevel, stokes, sysvel_offset, angle, \
           full_width, distance, nchan_low, cliplevel_low, nchan_high, cliplevel_high, prune_by_npix, \
           prune_by_fracbeam, expand_by_fracbeam, expand_by_npix, expand_by_nchan, inclination, eccentricity, \
           figsize