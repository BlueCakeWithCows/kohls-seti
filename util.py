from astropy.coordinates import EarthLocation, Longitude, Latitude
import astropy.units as u
import numpy as np


def get_earth_grid(n):
    grid = Grid(n)
    for i in range(0, n):
        for j in range(0, n):
            earth_loc = EarthLocation.from_geodetic(Longitude(j * 360 / n - 180, u.degree),
                                                    Latitude(i * 180 / n - 90 + 180 / n / 2, u.degree))
            grid.set(earth_loc, i, j)
    return grid


def nearest(trial_sols, observation):
    def center(array):
        tmp = sum(~np.isnan(array))
        if tmp:
            average = np.nansum(array) / sum(~np.isnan(array))
            return array-average
        return array

    obs = center(observation)
    n = 20
    errors = [0] * len(trial_sols)

    for i in range(0, n ** 2):
        result = center(trial_sols[i]) - obs
        tmp = sum(~np.isnan(result))
        if tmp:
            result = (result ** 2) / tmp
        err = np.nansum(result).value
        not_matched = sum(np.isnan(obs) * (~np.isnan(trial_sols[i]))) + sum((~np.isnan(obs)) * np.isnan(trial_sols[i]))
        if not_matched >= 2:
            err = np.inf
        errors[i] = err
    indices = np.argsort(errors)
    return indices


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
        else:
            self.lst[x] = val

    def getAdjacent(self, x, y):
        return [self.getIndex(v_x, v_y) for v_x, v_y in [(x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)]]
