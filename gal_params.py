def parameters(galaxy):

    if galaxy == 'NGC4713':
        centre_y = 25
        centre_x = 32
        size = 150
        vrange = 80
        vrange_2 = 50
        pbcor = False
        cliplevel = 4
        stokes = False
        start = 38
        stop = 60
        dv = 10

    elif galaxy == 'NGC4568':
        centre_y = 23
        centre_x = 27
        size = 60
        vrange = 120
        vrange_2 = 50
        pbcor = False
        cliplevel = 2
        stokes = False
        start = 21
        stop = 55
        dv = 10

    elif galaxy == 'NGC4189':
        centre_y = 32
        centre_x = 20
        size = 60
        vrange = 120
        vrange_2 = 120
        pbcor = False
        cliplevel = 2
        stokes = False
        start = 27
        stop = 56
        dv = 10


    return centre_x, centre_y, size, vrange, vrange_2, pbcor, cliplevel, stokes, start, stop, dv