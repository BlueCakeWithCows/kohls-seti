from astropy.coordinates import EarthLocation, Longitude, Latitude
import astropy.units as u


def get_earth_grid(n):
    grid = Grid(n)
    for i in range(0, n):
        for j in range(0, n):
            earth_loc = EarthLocation.from_geodetic(Longitude(j * 360 / n - 180, u.degree),
                                                    Latitude(i * 180 / n - 90, u.degree))
            grid.set(earth_loc, i, j)
    return grid


class Grid:
    def __init__(self, n):
        self.n = n
        self.lst = [0] * (n ** 2)

    def getIndex(self, x, y):
        n = self.n
        x = (x % n + n) % n
        y = (y % n + n) % n
        return x * n + y

    def get(self, x, y=None):
        if y is not None:
            return self.lst[self.getIndex(x, y)]
        return self.lst[x]

    def set(self, val, x, y=None):
        if y is not None:
            self.lst[self.getIndex(x, y)] = val
        self.lst[x] = val

    def getAdjacent(self, x, y):
        return [self.getIndex(v_x, v_y) for v_x, v_y in [(x-1, y), (x+1, y), (x, y-1), (x, y+1)]]
