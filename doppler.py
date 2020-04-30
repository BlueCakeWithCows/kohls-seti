#
#
#
import os
import pickle, dill
from math import ceil

import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import CartesianRepresentation, CartesianDifferential, SkyCoord, GCRS, AltAz
from astropy.coordinates import get_body_barycentric_posvel
from astropy.timeseries import TimeSeries

CONST_C = 299792458.  # m/s


def relativistic_doppler(v, costheta):
    return 1 / ((1 - v ** 2 / CONST_C ** 2) ** (1 / 2) * (1 + v / CONST_C * costheta))


class DopplerSystemEME:
    """
        Requires JPL emphemeris. pip install jplephem
        from astropy.coordinates import solar_system_ephemeris; solar_system_ephemeris.set('jpl')

        Calculates and plots doppler shift for EME system.
        
    """

    def __init__(self, emitters, receiver, start_time, timestep, timespan, doppler=relativistic_doppler,
                 only_visible=True, signal=None):
        """

        @param emitters: List of EarthLocation
        @param receiver: EarthLocation
        @param start_time: Time
        @param timestep: Quantity with time unit
        @param timespan: Quantity with time unit
        @param doppler: Function to calculate doppler shift; use default
        @param only_visible: Bool, filter doppler with moon visibility
        @param signal: None or Number, if none plots doppler factor else plots doppler shift.
        """
        self.emitters = emitters
        self.receiver = receiver
        self.only_visible = only_visible
        self.doppler = doppler

        self.start_time = start_time
        self.timestep = timestep
        self.timespan = timespan
        self.samples = ceil(timespan.to('hour') / timestep.to('hour'))
        self.times = TimeSeries(time_start=self.start_time, time_delta=self.timestep, n_samples=self.samples)['time']

        self._last_dopplers = None
        self._last_moon_emitters_altitudes = None
        self._last_moon_receiver_altitudes = None
        self.signal = signal

    def update(self):
        self.get_dopplers()
        self.get_moon_altitudes_emitters()
        self.get_moon_altitudes_receiver()

    def _get_dopplers(self, time):
        shifts = np.empty(shape=(len(self.emitters), len(time)))
        shifts[:, :] = 1

        # get current moon
        moon_loc, moon_vel = get_body_barycentric_posvel('moon', time)
        moon_loc = SkyCoord(
            representation_type='cartesian', differential_type=CartesianDifferential, unit='m', frame='icrs',
            obstime=time, x=moon_loc.x, y=moon_loc.y, z=moon_loc.z, v_x=moon_vel.x, v_y=moon_vel.y, v_z=moon_vel.z)
        moon = moon_loc.transform_to(GCRS(obstime=time))

        # get current receiver
        receiver_loc, receiver_vel = self.receiver.get_gcrs_posvel(time)
        receiver = SkyCoord(
            representation_type='cartesian', differential_type=CartesianDifferential, unit='m', frame='gcrs',
            obstime=time, x=receiver_loc.x, y=receiver_loc.y, z=receiver_loc.z,
            v_x=receiver_vel.x, v_y=receiver_vel.y, v_z=receiver_vel.z)

        # get moon altitude
        moon_altaz = moon.transform_to(AltAz(obstime=time, location=receiver))
        if self.only_visible:
            shifts[:, ] = np.greater(moon_altaz.alt, 0, dtype=np.float64)

        # get current emitters
        emitters = []
        for e in self.emitters:
            emitter_loc, emitter_vel = e.get_gcrs_posvel(time)
            emitter = SkyCoord(
                representation_type='cartesian', differential_type=CartesianDifferential, unit='m', frame='gcrs',
                obstime=time, x=emitter_loc.x, y=emitter_loc.y, z=emitter_loc.z,
                v_x=emitter_vel.x, v_y=emitter_vel.y, v_z=emitter_vel.z)
            emitters.append(emitter)

        # calculate moon to earth doppler shift
        delta_pos_m_e = CartesianRepresentation(x=(moon.cartesian.xyz - receiver.cartesian.xyz), unit='m')
        delta_v_m_e = CartesianDifferential(d_x=(moon.velocity - receiver.velocity).get_d_xyz(),
                                            unit='m/s').to_cartesian()
        theta_m_e = delta_pos_m_e.dot(delta_v_m_e) / (delta_pos_m_e.norm() * delta_v_m_e.norm())
        delta_v_m_e_norm = delta_v_m_e.norm().value
        shifted_m_e = self.doppler(delta_v_m_e_norm, theta_m_e)

        # calculate earth to moon doppler shifts

        for i, emitter in enumerate(emitters):
            delta_pos_e_m = CartesianRepresentation(x=(emitter.cartesian.xyz - moon.cartesian.xyz), unit='m')
            delta_v_e_m = CartesianDifferential(d_x=(emitter.velocity - moon.velocity).get_d_xyz(),
                                                unit='m/s').to_cartesian()
            theta_e_m = delta_pos_e_m.dot(delta_v_e_m) / (delta_pos_e_m.norm() * delta_v_e_m.norm())
            delta_v_e_m_norm = delta_v_e_m.norm().value
            if self.only_visible:
                moon_altaz_2 = moon.transform_to(AltAz(obstime=time, location=emitter))
                shifts[i] = shifts[i] * self.doppler(delta_v_e_m_norm, theta_e_m) * np.greater(moon_altaz_2.alt, 0,
                                                                                               dtype=np.float64)
            else:
                shifts[i] = self.doppler(delta_v_e_m_norm, theta_e_m)

        shifts = shifts * shifted_m_e
        shifts[shifts == 0] = np.nan

        if self.signal:
            shifts[:] = shifts[:] * self.signal - self.signal

        return shifts

    def get_dopplers(self):
        if self._last_dopplers is None:
            self._last_dopplers = self._get_dopplers(self.times)
        return self._last_dopplers

    def _get_moon_altitudes(self, time, pos):
        # get current moon
        moon_loc, moon_vel = get_body_barycentric_posvel('moon', time)
        moon_loc = SkyCoord(
            representation_type='cartesian', differential_type=CartesianDifferential, unit='m', frame='icrs',
            obstime=time, x=moon_loc.x, y=moon_loc.y, z=moon_loc.z, v_x=moon_vel.x, v_y=moon_vel.y, v_z=moon_vel.z)
        moon = moon_loc.transform_to(GCRS(obstime=time))

        # get current receiver
        receiver_loc, receiver_vel = pos.get_gcrs_posvel(time)
        receiver = SkyCoord(
            representation_type='cartesian', differential_type=CartesianDifferential, unit='m', frame='gcrs',
            obstime=time, x=receiver_loc.x, y=receiver_loc.y, z=receiver_loc.z,
            v_x=receiver_vel.x, v_y=receiver_vel.y, v_z=receiver_vel.z)
        moon_altaz = moon.transform_to(AltAz(obstime=time, location=receiver))
        return moon_altaz.alt

    def get_moon_altitudes_receiver(self):
        if self._last_moon_receiver_altitudes is None:
            self._last_moon_receiver_altitudes = self._get_moon_altitudes(self.times, self.receiver)
        return self._last_moon_receiver_altitudes

    def get_moon_altitudes_emitters(self):
        if self._last_moon_emitters_altitudes is None:
            altitudes = []
            for emitter in self.emitters:
                altitudes.append(self._get_moon_altitudes(self.times, emitter))
            self._last_moon_emitters_altitudes = altitudes
        return self._last_moon_emitters_altitudes

    def save(self, path):
        save(self, path)

    def save_plots(self, path):
        os.makedirs(path, exist_ok=True)
        receiver = self.receiver
        x = range(0, self.samples) * self.timestep.to('hour').value
        start_string = str(self.start_time)
        receiver_moon_altitudes = self._last_moon_receiver_altitudes
        emitters_moon_altitudes = self._last_moon_emitters_altitudes
        for idx, emitter in enumerate(self.emitters):
            fig, (ax, ax2) = plt.subplots(2)
            titlel = '(%f, %f) to (%f, %f)' % (
                emitter.lat.value, emitter.lon.value, receiver.lat.value, receiver.lon.value)
            if self.signal:
                title = titlel + '\n{:3e} Hz'.format(self.signal)
            ax.set_title(title, pad='20')
            ax.set_xlim(0, x[-1])
            ax2.set_xlabel('Hours Since %s' % start_string)
            ax.set_ylim(1 - 4 * 1e-6, 1 + 4 * 1e-6)
            if self.signal:
                tmp = self.signal / (10e6 / 40)
                ax.set_ylim(-tmp, tmp)

            label = 'Doppler Shift (Hz)' if self.signal else 'Doppler Factor'
            ax.plot(x, self._last_dopplers[idx], label=label)
            ax.set_ylabel(label)
            # ax2= fig.gca()
            ax2.plot(x, receiver_moon_altitudes, color='g', label='Moon Alt, receiver')
            ax2.plot(x, emitters_moon_altitudes[idx], color='r', label='Moon Alt, Emitter')
            ax2.set_ylabel('Altitude (Deg)')
            ax2.set_xlim(0, x[-1])
            ax2.set_ylim(0, 90)
            lgd = fig.legend(bbox_to_anchor=(.95, .8), loc="upper left", borderpad=1.2)
            fig.savefig(path + titlel + '.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
            plt.close(fig)

    def get_plot(self, idx):
        receiver = self.receiver
        x = range(0, self.samples) * self.timestep.to('hour').value
        start_string = str(self.start_time)
        receiver_moon_altitudes = self._last_moon_receiver_altitudes
        emitters_moon_altitudes = self._last_moon_emitters_altitudes
        emitter = self.emitters[idx]
        fig, (ax, ax2) = plt.subplots(2)
        title = '(%f, %f) to (%f, %f)' % (emitter.lat.value, emitter.lon.value, receiver.lat.value, receiver.lon.value)
        if self.signal:
            title += '\n{:3e} Hz'.format(self.signal)
        ax.set_title(title, pad='20')
        ax.set_xlim(0, x[-1])
        ax2.set_xlabel('Hours Since %s' % start_string)
        ax.set_ylim(1 - 4 * 1e-6, 1 + 4 * 1e-6)
        if self.signal:
            tmp = self.signal / (10e6 / 40)
            ax.set_ylim(-tmp, tmp)

        label = 'Doppler Shift (Hz)' if self.signal else 'Doppler Factor'
        ax.plot(x, self._last_dopplers[idx], label=label)
        ax.set_ylabel(label)
        # ax2= fig.gca()
        ax2.plot(x, receiver_moon_altitudes, color='g', label='Moon Alt, receiver')
        ax2.plot(x, emitters_moon_altitudes[idx], color='r', label='Moon Alt, Emitter')
        ax2.set_ylabel('Altitude (Deg)')
        ax2.set_xlim(0, x[-1])
        ax2.set_ylim(0, 90)
        lgd = fig.legend(bbox_to_anchor=(.95, .8), loc="upper left", borderpad=1.2)


def save(obj, path):
    with open(path, 'wb') as pickle_file:
        pickle.dump(obj, pickle_file)


def load(filename):
    return pickle.load(filename)
