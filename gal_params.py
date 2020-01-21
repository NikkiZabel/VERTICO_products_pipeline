def parameters(galaxy):

    if galaxy == 'NGC4713':
        centre_y = 25
        centre_x = 32
        size = 70
        vrange = 80
        vrange_2 = 50
        pbcor = False
        cliplevel = 4
        stokes = False
        start = 24
        stop = 56
        dv = 10
        angle = 0
        centre_pvd = 25
        sysvel = 25
        distance = 14.8

    elif galaxy == 'NGC4568':
        centre_y = 23
        centre_x = 27
        size = 60
        vrange = 200
        vrange_2 = 75
        pbcor = False
        cliplevel = 4
        stokes = False
        start = 13
        stop = 63
        dv = 10
        angle = 110
        centre_pvd = 32
        sysvel = 19
        distance = 17.4

    elif galaxy == 'NGC4189':
        centre_y = 32
        centre_x = 20
        size = 60
        vrange = 200
        vrange_2 = 120
        pbcor = False
        cliplevel = 4
        stokes = False
        start = 19
        stop = 59
        dv = 10
        angle = 335
        centre_pvd = 33
        sysvel = 14
        distance = 32.1

    name = galaxy

    return name, centre_x, centre_y, size, vrange, vrange_2, pbcor, cliplevel, stokes, start, stop, dv, angle, \
           centre_pvd, sysvel, distance
