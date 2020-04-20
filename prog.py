#
#
#
import numpy as np
from astropy.coordinates import get_body_barycentric,get_body_barycentric_posvel, get_body, get_moon
from astropy.coordinates import CartesianRepresentation, CartesianDifferential, SkyCoord, GCRS

CONST_C = 299792458 #m/s


def relativistic_doppler(v, costheta):
    return 1/((1-v**2/CONST_C**2)**(1/2) *  (1 + v/CONST_C * costheta))


class DopplerSystemEME:
    """
        Requires JPL emphemeris. pip install jplephem
        from astropy.coordinates import solar_system_ephemeris; solar_system_ephemeris.set('jpl')
        
        
    """
    
    def __init__(self, emitters, reciever, doppler = relativistic_doppler):
        self.emitters = emitters
        self.reciever = reciever
        self.doppler = doppler
        
    def _get_dopplers(self, time):
    
        #get current moon
        moon_loc, moon_vel = get_body_barycentric_posvel('moon', time)
        moon_loc = SkyCoord(
        representation_type = 'cartesian', differential_type=CartesianDifferential,unit='m',frame='icrs', 
        obstime =time, x = moon_loc.x, y= moon_loc.y, z = moon_loc.z, v_x = moon_vel.x, v_y = moon_vel.y, v_z = moon_vel.z)
        moon = moon_loc.transform_to(GCRS(obstime=time))
        
        #get current emitters
        emitters = []
        for e in self.emitters:
            emitter_loc, emitter_vel = e.get_gcrs_posvel(time)
            emitter = SkyCoord(
                representation_type = 'cartesian', differential_type=CartesianDifferential, unit = 'm',frame = 'gcrs', 
                obstime = time, x = emitter_loc.x, y= emitter_loc.y, z = emitter_loc.z, 
                v_x = emitter_vel.x, v_y = emitter_vel.y, v_z = emitter_vel.z)
            emitters.append(emitter)
        
        #get current reciever
        reciever_loc, reciever_vel = self.reciever.get_gcrs_posvel(time)
        reciever = SkyCoord(
        representation_type = 'cartesian', differential_type=CartesianDifferential,unit='m',frame='gcrs', 
        obstime =time, x = reciever_loc.x, y= reciever_loc.y, z = reciever_loc.z, 
        v_x = reciever_vel.x, v_y = reciever_vel.y, v_z = reciever_vel.z)

        #calculate moon to earth doppler shift
        delta_pos_m_e = CartesianRepresentation(x=(moon.cartesian.xyz-reciever.cartesian.xyz), unit='m')
        delta_v_m_e = CartesianDifferential(d_x = (moon.velocity-reciever.velocity).get_d_xyz(), unit='m/s').to_cartesian()
        theta_m_e = delta_pos_m_e.dot(delta_v_m_e)/(delta_pos_m_e.norm() * delta_v_m_e.norm())
        delta_v_m_e_norm= delta_v_m_e.norm().value
        shifted_m_e = self.doppler(delta_v_m_e_norm, theta_m_e)
        
        #calculate earth to moon doppler shifts 
        shifts = np.zeros(len(emitters))
        for i, emitter in enumerate(emitters):
            delta_pos_e_m = CartesianRepresentation(x=(emitter.cartesian.xyz-moon.cartesian.xyz), unit='m')
            delta_v_e_m = CartesianDifferential(d_x = (emitter.velocity-moon.velocity).get_d_xyz(), unit='m/s').to_cartesian()
            theta_e_m = delta_pos_e_m.dot(delta_v_e_m)/(delta_pos_e_m.norm() * delta_v_e_m.norm())
            delta_v_e_m_norm= delta_v_e_m.norm().value
            shifts[i] = self.doppler(delta_v_e_m_norm, theta_e_m)
       
        #calculate total shift
        shifts = shifts * shifted_m_e
        return shifts
    
    def get_dopplers(self, times):
        result = np.zeros(shape=(len(self.emitters), len(times)))
        for i, time in enumerate(times['time']):
            result[:, i] = self._get_dopplers(time)
        return result
    
    