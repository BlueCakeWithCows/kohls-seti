#
#
#
import numpy as np
from astropy.coordinates import get_body_barycentric,get_body_barycentric_posvel, get_body, get_moon
from astropy.coordinates import CartesianRepresentation, CartesianDifferential, SkyCoord, GCRS, AltAz
import astropy.units as u
CONST_C = 299792458. #m/s

from multiprocessing.dummy import Pool as ThreadPool

def relativistic_doppler(v, costheta):
    return 1/((1-v**2/CONST_C**2)**(1/2) *  (1 + v/CONST_C * costheta))


class DopplerSystemEME:
    """
        Requires JPL emphemeris. pip install jplephem
        from astropy.coordinates import solar_system_ephemeris; solar_system_ephemeris.set('jpl')
        
        
    """
    
    def __init__(self, emitters, reciever, doppler = relativistic_doppler, only_visible = True):
        self.emitters = emitters
        self.reciever = reciever
        self.only_visible = only_visible
        self.doppler = doppler
        
    def _get_dopplers(self, time):
        #print(len(self.emitters.shape))
        shifts = np.empty(shape=(len(self.emitters), len(time)))
        shifts[:,:] = 1
            
        #get current moon
        moon_loc, moon_vel = get_body_barycentric_posvel('moon', time)
        moon_loc = SkyCoord(
        representation_type = 'cartesian', differential_type=CartesianDifferential,unit='m',frame='icrs', 
        obstime =time, x = moon_loc.x, y= moon_loc.y, z = moon_loc.z, v_x = moon_vel.x, v_y = moon_vel.y, v_z = moon_vel.z)
        moon = moon_loc.transform_to(GCRS(obstime=time))
        
        #get current reciever
        reciever_loc, reciever_vel = self.reciever.get_gcrs_posvel(time)
        reciever = SkyCoord(
        representation_type = 'cartesian', differential_type=CartesianDifferential,unit='m',frame='gcrs', 
        obstime =time, x = reciever_loc.x, y= reciever_loc.y, z = reciever_loc.z, 
        v_x = reciever_vel.x, v_y = reciever_vel.y, v_z = reciever_vel.z)
        
        #get moon altitude
        moon_altaz = moon.transform_to(AltAz(obstime=time, location = reciever))
        if self.only_visible:
            shifts[:,] = np.greater(moon_altaz.alt, 0, dtype=np.float64)
            
            
        
        #get current emitters
        emitters = []
        for e in self.emitters:
            emitter_loc, emitter_vel = e.get_gcrs_posvel(time)
            emitter = SkyCoord(
                representation_type = 'cartesian', differential_type=CartesianDifferential, unit = 'm',frame = 'gcrs', 
                obstime = time, x = emitter_loc.x, y= emitter_loc.y, z = emitter_loc.z, 
                v_x = emitter_vel.x, v_y = emitter_vel.y, v_z = emitter_vel.z)
            emitters.append(emitter)

        #calculate moon to earth doppler shift
        delta_pos_m_e = CartesianRepresentation(x=(moon.cartesian.xyz-reciever.cartesian.xyz), unit='m')
        delta_v_m_e = CartesianDifferential(d_x = (moon.velocity-reciever.velocity).get_d_xyz(), unit='m/s').to_cartesian()
        theta_m_e = delta_pos_m_e.dot(delta_v_m_e)/(delta_pos_m_e.norm() * delta_v_m_e.norm())
        delta_v_m_e_norm= delta_v_m_e.norm().value
        shifted_m_e = self.doppler(delta_v_m_e_norm, theta_m_e)
        
        #calculate earth to moon doppler shifts 
        
        temp = np.zeros(shape=(len(time)))
        for i, emitter in enumerate(emitters):
            delta_pos_e_m = CartesianRepresentation(x=(emitter.cartesian.xyz-moon.cartesian.xyz), unit='m')
            delta_v_e_m = CartesianDifferential(d_x = (emitter.velocity-moon.velocity).get_d_xyz(), unit='m/s').to_cartesian()
            theta_e_m = delta_pos_e_m.dot(delta_v_e_m)/(delta_pos_e_m.norm() * delta_v_e_m.norm())
            delta_v_e_m_norm= delta_v_e_m.norm().value
            moon_altaz_2 = moon.transform_to(AltAz(obstime=time, location = emitter))
            temp[:]=  self.doppler(delta_v_e_m_norm, theta_e_m)
            if self.only_visible:
                shifts[i] = temp[i] * np.greater(moon_altaz.alt, 0, dtype=np.float64)
            else: 
                shifts[i] = temp[i]
       
        shifts = shifts * shifted_m_e
        shifts[shifts==0] = np.nan
        return shifts
    
    def get_dopplers(self, times, threads = 1):
        result = np.zeros(shape=(len(self.emitters), len(times)))
        return self._get_dopplers(times)
    
    