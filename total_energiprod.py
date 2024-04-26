PANEL_COUNT = 38

def daterange(start_date, end_date):
    for n in range(int((end_date - start_date).days)):
        yield start_date + timedelta(n)

def get_daily_energy(start_date, end_date, lat, lon, alt, theta, phi):
    result = []
    for single_date in daterange(start_date, end_date): 
        date_string = single_date.strftime("%Y-%m-%d")
        solpos = get_solar_position(date_string, date_string, lat, lon, alt, "Europe/Copenhagen", dt = "Min")
        flux = PANEL_COUNT * solar_flux(L,B,S0,A0, PANEL_EFFICIENCY, theta, phi, solpos )
        joule = integrate.simps(flux, dx=60)
        result.append(joule)
    return result
    
    
pdStartDates = pd.date_range(date(2023,1,1), date(2023,12,31), freq='MS')
pdEndDates = pdStartDates + pd.offsets.MonthEnd(1)
startDates = np.array([d.strftime('%Y-%m-%d') for d in pdStartDates])
endDates = np.array([d.strftime('%Y-%m-%d') for d in pdEndDates])

print(endDates)

energy_optimized = []
for i in range(12):
    energy_optimized += get_daily_energy(pdStartDates[i], pdEndDates[i], 55.781, 12.5234, 10, thetas[i], phis[i])
    
plt.plot(energy_optimized)
sum(energy_optimized)
