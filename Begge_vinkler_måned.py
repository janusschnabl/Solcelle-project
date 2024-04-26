import pandas as pd
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from helpers import *

def get_monthly_prices():
    prices = pd.read_json("Elspotprices.json")
    prices = prices[prices.PriceArea == "DK2"]
    prices["HourDK"] = pd.to_datetime(prices["HourDK"])
    prices["HourUTC"] = pd.to_datetime(prices["HourUTC"])
    monthly_prices = []
    for i in range(1,13):
        month = prices[prices["HourDK"].dt.month == i]
        month.sort_values(by=["HourDK","HourUTC"])
        monthly_prices.append(np.array(month["SpotPriceDKK"] / 3600000000))
    return monthly_prices

# Precompute solar positions
def monthly_angles(longtitude, altitude, elevation, tz, optimize_price=True):
    pdStartDates = pd.date_range(start='2023-01-01', end='2023-12-01', freq='MS')
    pdEndDates = pdStartDates + pd.offsets.MonthEnd(1)
    startDates = np.array([d.strftime('%Y-%m-%d') for d in pdStartDates])
    endDates = np.array([d.strftime('%Y-%m-%d') for d in pdEndDates])
    monthly_prices = get_monthly_prices()
    result = []
    for i in range(12):
        solpos = get_solar_position(startDates[i], endDates[i], longtitude, altitude, elevation, tz)
        energy = np.zeros((91, 360))
        for theta in range(91):
            for phi in range(360):
                flux = solar_flux(L, B, 1100, 0.5,  theta, phi, solpos)
                price = monthly_prices[i] * flux
                energy[theta][phi] = integrate.simps(price if optimize_price else flux, dx=3600)
        max_energy_index = np.unravel_index(np.argmax(energy), energy.shape)
        result.append(max_energy_index)

        print(f"Maximum energy produced at Theta = {max_energy_index[0]} degrees and Phi = {max_energy_index[1]} degrees")
    return result

angles = monthly_angles(55.7861, 12.5234, 10, "Europe/Copenhagen",optimize_price=False)

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

