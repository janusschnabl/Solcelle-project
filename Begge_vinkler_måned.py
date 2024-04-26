import pandas as pd
import pvlib
import numpy as np
from pvlib.location import Location
from scipy import integrate
import matplotlib.pyplot as plt
from datetime import date, timedelta

PANEL_EFFICIENCY = 0.214
A0 = 0.5
S0 = 1100
L = 2.384
B = 1.303 
PANEL_COUNT = 38

# Projectionen af solens stråler på et panel
def solar_panel_projection(theta_s, phi_s, theta_p, phi_p):
    proj = np.sin(theta_p) * np.sin(theta_s) * np.cos(phi_p - phi_s) + np.cos(theta_p) * np.cos(theta_s)
    return np.where(proj < 0, 0,proj)

# Solens position
def get_solar_position(start_dato, slut_dato,latitude, longtitude, altitude, tz, dt = "H"):
    """dt = 'H' or 'Min'"""
    tidszone = "Europe/Copenhagen"
    
    site = Location(latitude, longtitude, tz, altitude)

    # Definition of a time range of simulation
    times = pd.date_range (start_dato + " 00:00:00", slut_dato + " 23:59:00", inclusive="left", freq=dt, tz=tidszone)
    
    # Estimate Solar Position with the 'Location' object
    solpos = site.get_solarposition(times)
    return solpos

# Udregner fluxen ud fra solens position vinkel og meget mere:
def solar_flux(Længde, bredde, S_0, A_0, PANEL_EFFICIENCY, theta_panel, phi_panel, solpos):
    """Returns time,x,y,z as np arrays"""
    # Convert zenith and azimuth from degrees to radians
    theta = np.deg2rad(np.array(solpos["zenith"]))
    phi = np.deg2rad(np.array(solpos["azimuth"]))
    proj = solar_panel_projection(theta, phi, np.deg2rad(theta_panel), np.deg2rad(phi_panel))

    valid_mask = (theta >= 0) & (theta <= np.pi / 2)
    proj = np.where(valid_mask, proj, 0)
    
    flux = Længde * bredde * S_0 * A_0 * proj * PANEL_EFFICIENCY
    
    return flux

# Precompute solar positions
def monthly_angles(longtitude, altitude, elevation, tz, year):
    pdStartDates = pd.date_range(start='2024-01-01', end='2024-12-01', freq='MS')
    pdEndDates = pdStartDates + pd.offsets.MonthEnd(1)
    startDates = np.array([d.strftime('%Y-%m-%d') for d in pdStartDates])
    endDates = np.array([d.strftime('%Y-%m-%d') for d in pdEndDates])
    result = []
    for i in range(12):
        solpos = get_solar_position(startDates[i], endDates[i], longtitude, altitude, elevation, tz)
        energy = np.zeros((91, 360))
        for theta in range(91):
            for phi in range(360):
                flux = solar_flux(L, B, 1100, 0.5, PANEL_EFFICIENCY,  theta, phi, solpos)
                energy[theta][phi] = integrate.simps(flux, dx=3600)
        max_energy_index = np.unravel_index(np.argmax(energy), energy.shape)
        result.append(max_energy_index)

        print(f"Maximum energy produced at Theta = {max_energy_index[0]} degrees and Phi = {max_energy_index[1]} degrees")
    return result

angles = monthly_angles(55.7861, 12.5234, 10, "Europe/Copenhagen", 2024)

# Extracting theta and phi angles
thetas, phis = zip(*angles)

# Plotting the optimal angles
months = range(1, 13)
plt.figure(figsize=(10, 5))
plt.plot(months, thetas, label='Theta')
plt.plot(months, phis, label='Phi')
plt.xlabel('Month')
plt.ylabel('Angle (degrees)')
plt.title('Optimal Solar Panel Angles by Month')
plt.legend()
plt.grid(True)
plt.xticks(months)
plt.show()

