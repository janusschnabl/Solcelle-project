from helpers import *
#from Begge_vinkler_måned import monthly_angles

def daterange(start_date, end_date):
    for n in range(int((end_date - start_date).days)):
        yield start_date + timedelta(n)

def get_daily_prices():
    prices = pd.read_json("Elspotprices.json")
    prices = prices[prices.PriceArea == "DK2"]
    prices["HourDK"] = pd.to_datetime(prices["HourDK"])
    prices["HourUTC"] = pd.to_datetime(prices["HourUTC"])
    monthly_prices = []
    for i in range(1,13):
        month_prices = []
        month = prices[prices["HourDK"].dt.month == i]
        for j in range(1,max(month["HourDK"].dt.day)+1):
            day = month[month["HourDK"].dt.day == j]
            day.sort_values(by=["HourDK","HourUTC"])
            month_prices.append(np.array(day["SpotPriceDKK"] / 3600000000))
        monthly_prices.append(month_prices)
    return monthly_prices


def get_daily_energy(start_date, end_date, lat, lon, alt, theta, phi):
    daily_prices = get_daily_prices()
    result = []
    for single_date in daterange(start_date, end_date):
        date_string = single_date.strftime("%Y-%m-%d")
        solpos = get_solar_position(date_string, date_string, lat, lon, alt, "Europe/Copenhagen", dt = "H")
        flux = PANEL_COUNT * solar_flux(L,B,S0,A0, theta, phi, solpos )
        price = flux * daily_prices[single_date.month-1][single_date.day-1]
        joule = integrate.simps(price, dx=3600)
        result.append(joule)
    return result


pdStartDates = pd.date_range(date(2023,1,1), date(2023,12,31), freq='MS')
pdEndDates = pdStartDates + pd.offsets.MonthEnd(1)
startDates = np.array([d.strftime('%Y-%m-%d') for d in pdStartDates])
endDates = np.array([d.strftime('%Y-%m-%d') for d in pdEndDates])

print(endDates)

angles = monthly_angles(55.7861, 12.5234, 10, "Europe/Copenhagen",optimize_price=True)
thetas, phis = zip(*angles)

energy_optimized = []
for i in range(12):
    energy_optimized += get_daily_energy(pdStartDates[i], pdEndDates[i], 55.781, 12.5234, 10, thetas[i], phis[i])

plt.plot(energy_optimized)
print(energy_optimized)
