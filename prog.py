#
#
#
import numpy as np
from astropy.coordinates import get_body_barycentric,get_body_barycentric_posvel, get_body, get_moon
from astropy.coordinates import CartesianRepresentation, CartesianDifferential, SkyCoord, GCRS, AltAz, EarthLocation, Longitude, Latitude   
from astropy.timeseries import TimeSeries
import astropy.units as u
import os
import matplotlib.pyplot as plt
from math import ceil
CONST_C = 299792458. #m/s

from multiprocessing.dummy import Pool as ThreadPool

def relativistic_doppler(v, costheta):
    return 1/((1-v**2/CONST_C**2)**(1/2) *  (1 + v/CONST_C * costheta))


def get_earth_location_grid(n):
    emitters = []
    for i in range(0, n):
        tmp = []
        for j in range(0,n+1):
            earth_loc = EarthLocation.from_geodetic(Longitude(i * 360 / n, u.degree), Latitude(j * 180 /n -90, u.degree))
            tmp.append(earth_loc)
        emitters.append(tmp)
    return [item for sublist in emitters for item in sublist], emitters


class DopplerSystemEME:
    """
        Requires JPL emphemeris. pip install jplephem
        from astropy.coordinates import solar_system_ephemeris; solar_system_ephemeris.set('jpl')
        
        
    """
    
    def __init__(self, emitters, reciever, start_time, timestep, timespan, doppler = relativistic_doppler, only_visible = True, signal = None):
        self.emitters = emitters
        self.reciever = reciever
        self.only_visible = only_visible
        self.doppler = doppler
        
        self.start_time = start_time
        self.timestep = timestep
        self.timespan = timespan
        self.samples = ceil(timespan.to('hour')/timestep.to('hour'))
        self.times = TimeSeries(time_start= self.start_time, time_delta = self.timestep, n_samples= self.samples)['time']
        
        self._last_dopplers = None
        self._last_moon_emitters_altitudes = None
        self._last_moon_reciever_altitudes = None
        self.signal = signal
        
    def update(self):
        self.get_dopplers()
        self.get_moon_altitudes_emitters()
        self.get_moon_altitudes_reciever()
        
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
        
        for i, emitter in enumerate(emitters):
            delta_pos_e_m = CartesianRepresentation(x=(emitter.cartesian.xyz-moon.cartesian.xyz), unit='m')
            delta_v_e_m = CartesianDifferential(d_x = (emitter.velocity-moon.velocity).get_d_xyz(), unit='m/s').to_cartesian()
            theta_e_m = delta_pos_e_m.dot(delta_v_e_m)/(delta_pos_e_m.norm() * delta_v_e_m.norm())
            delta_v_e_m_norm= delta_v_e_m.norm().value
            if self.only_visible:
                moon_altaz_2 = moon.transform_to(AltAz(obstime=time, location = emitter))
                shifts[i] = shifts[i] * self.doppler(delta_v_e_m_norm, theta_e_m) * np.greater(moon_altaz_2.alt, 0, dtype=np.float64)
            else: 
                shifts[i] = self.doppler(delta_v_e_m_norm, theta_e_m)
       
        shifts = shifts * shifted_m_e
        shifts[shifts==0] = np.nan
        
        if self.signal:
            shifts[:]  = shifts[:]*self.signal-self.signal
        
        return shifts
    
    def get_dopplers(self):
        if self._last_dopplers is None:
            self._last_dopplers = self._get_dopplers(self.times)
        return self._last_dopplers
        
    def _get_moon_altitudes(self, time, pos):
        #get current moon
        moon_loc, moon_vel = get_body_barycentric_posvel('moon', time)
        moon_loc = SkyCoord(
        representation_type = 'cartesian', differential_type=CartesianDifferential,unit='m',frame='icrs', 
        obstime =time, x = moon_loc.x, y= moon_loc.y, z = moon_loc.z, v_x = moon_vel.x, v_y = moon_vel.y, v_z = moon_vel.z)
        moon = moon_loc.transform_to(GCRS(obstime=time))
        
        #get current reciever
        reciever_loc, reciever_vel = pos.get_gcrs_posvel(time)
        reciever = SkyCoord(
        representation_type = 'cartesian', differential_type=CartesianDifferential,unit='m',frame='gcrs', 
        obstime =time, x = reciever_loc.x, y= reciever_loc.y, z = reciever_loc.z, 
        v_x = reciever_vel.x, v_y = reciever_vel.y, v_z = reciever_vel.z)
        moon_altaz = moon.transform_to(AltAz(obstime=time, location = reciever))
        return moon_altaz.alt
    
    def get_moon_altitudes_reciever(self):
        if self._last_moon_reciever_altitudes is None:
            self._last_moon_reciever_altitudes = self._get_moon_altitudes(self.times, self.reciever)
        return self._last_moon_reciever_altitudes
        
    def get_moon_altitudes_emitters(self):
        if self._last_moon_emitters_altitudes is None:
            altitudes = []
            for emitter in self.emitters:
                altitudes.append(self._get_moon_altitudes(self.times, emitter))
            self._last_moon_emitters_altitudes = altitudes
        return self._last_moon_emitters_altitudes 
        
    def save_plots(self, path):
        os.makedirs(path, exist_ok=True)
        reciever = self.reciever
        x = range(0, self.samples)* self.timestep.to('hour').value
        start_string = str(self.start_time)
        reciever_moon_altitudes = self._last_moon_reciever_altitudes
        emitters_moon_altitudes = self._last_moon_emitters_altitudes
        for idx, emitter in enumerate(self.emitters):
            fig, (ax, ax2) = plt.subplots(2)
            title = '(%f, %f) to (%f, %f)'%(emitter.lat.value, emitter.lon.value, reciever.lat.value, reciever.lon.value)
            ax.set_title(title, pad = '20')
            ax.set_xlim(0,x[-1])
            ax2.set_xlabel('Hours Since %s'%start_string)
            ax.set_ylim(1 - 4 * 1e-6, 1 + 4* 1e-6)
            ax.plot(x, self._last_dopplers[idx], label = 'Doppler Factor')
            ax.set_ylabel('Doppler Factor')
            #ax2= fig.gca()
            ax2.plot(x, reciever_moon_altitudes, color ='g', label = 'Moon Alt, Reciever')
            ax2.plot(x, emitters_moon_altitudes[idx], color ='r', label = 'Moon Alt, Emitter')
            ax2.set_ylabel('Altitude (Deg)')
            ax2.set_xlim(0, x[-1])
            ax2.set_ylim(0,90)
            lgd = fig.legend(bbox_to_anchor=(.95,.8), loc="upper left", borderpad= 1.2)
            fig.savefig(path + title+'.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
            plt.close(fig)
            
    def get_plot(self, idx):
        reciever = self.reciever
        x = range(0, self.samples)* self.timestep.to('hour').value
        start_string = str(self.start_time)
        reciever_moon_altitudes = self._last_moon_reciever_altitudes
        emitters_moon_altitudes = self._last_moon_emitters_altitudes
        emitter = self.emitters[idx]
        fig, (ax, ax2) = plt.subplots(2)
        title = '(%f, %f) to (%f, %f)'%(emitter.lat.value, emitter.lon.value, reciever.lat.value, reciever.lon.value)
        if self.signal:
            title += '\n{:3e} Hz'.format(self.signal)
        ax.set_title(title, pad = '20')
        ax.set_xlim(0,x[-1])
        ax2.set_xlabel('Hours Since %s'%start_string)
        ax.set_ylim(1 - 4 * 1e-6, 1 + 4* 1e-6)
        if self.signal:
                    
                    tmp = self.signal/ (10e6/40)
                    ax.set_ylim(-tmp, tmp)
                    
        label = 'Doppler Shift (Hz)' if self.signal else 'Doppler Factor'
        ax.plot(x, self._last_dopplers[idx], label = label)
        ax.set_ylabel(label)
        #ax2= fig.gca()
        ax2.plot(x, reciever_moon_altitudes, color ='g', label = 'Moon Alt, Reciever')
        ax2.plot(x, emitters_moon_altitudes[idx], color ='r', label = 'Moon Alt, Emitter')
        ax2.set_ylabel('Altitude (Deg)')
        ax2.set_xlim(0, x[-1])
        ax2.set_ylim(0,90)
        lgd = fig.legend(bbox_to_anchor=(.95,.8), loc="upper left", borderpad= 1.2)
    
    