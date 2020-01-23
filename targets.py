import gal_params

class galaxies:
    def __init__(self, galname):

        name, centre_x, centre_y, size, vrange, vrange_2, cliplevel, stokes, start, stop, angle, centre_pvd, \
        sysvel, distance = gal_params.parameters(galname)

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
        self.angle = angle
        self.centre_pvd = centre_pvd
        self.sysvel = sysvel
        self.distance = distance
