from astropy.coordinates import EarthLocation, Longitude, Latitude
import astropy.units as u

def get_earth_grid(n):
    emitters = []
    for i in range(0, n):
        tmp = []
        for j in range(0, n):
            earth_loc = EarthLocation.from_geodetic(Longitude(i * 360 / n, u.degree),
                                                    Latitude(j * 180 / n - 90, u.degree))
            tmp.append(earth_loc)
        emitters.append(tmp)
    return [item for sublist in emitters for item in sublist], emitters
