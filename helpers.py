#indeholder alle de funktioner, som der står vi skal skrive.

#funktionsværdier i `f` i intervallet $[Low, Up]$ og angiver de tilhørende `t`-værdier.
def interval(f, t, lower, upper):
    for i in range(0, len(f)):
        if f[i] > lower and f[i] < upper:
            print("Data(f): ", f[i], "Time(t): ", t[i])
    return None
        

# Fortegn skift
def fortegnskift(f, t):
    for i in range(0, len(f) - 1):
        if f[i] * f[i+1] < 0:
            print("Skift ved: ", t[i])
    if f[len(f) - 1] * f[len(f)- 2] < 0:
        print("Skift ved: ", t[len(f)-1])
    return None


#Solens højeste punkt indefor interval af dato
def max_sol(start_dato, slut_dato, latitude, longtitude, place):
    delta_tid = "H" #"M"
    tidszone = "Europe/Copenhagen"
    
    site = Location(
        latitude, longtitude, "Europe/Copenhagen", 10, place
    )

    # Definition of a time range of simulation
    times = pd.date_range(
        start_dato + " 00:00:00", slut_dato + " 23:59:00", inclusive="left", freq=delta_tid, tz=tidszone
    )
    
    # Estimate Solar Position with the 'Location' object
    solpos = site.get_solarposition(times)
    
    f = np.array(solpos.loc[start_dato].elevation)

    return f.max()

#Solens x,y,z koordinator
def solar_position_to_xyz(start_dato, latitude, longtitude, altitude, tz):
    """Returns time,x,y,z as np arrays"""

    delta_tid = "H" #"M"
    tidszone = "Europe/Copenhagen"
    
    site = Location(latitude, longtitude, tz, altitude)

    # Definition of a time range of simulation
    times = pd.date_range (start_dato + " 00:00:00", start_dato + " 23:59:00", inclusive="left", freq=delta_tid, tz=tidszone)
    
    # Estimate Solar Position with the 'Location' object
    solpos = site.get_solarposition(times)
    
    # Get the Earth-Sun distance and convert to meters
    r = np.array(nrel_earthsun_distance(times) * 149597870700)  # 1 AU in meters
    
    # Convert zenith and azimuth from degrees to radians
    theta = np.deg2rad(np.array(solpos.loc[start_dato].zenith))
    phi = np.deg2rad(np.array(solpos.loc[start_dato].azimuth))
    
    # Calculate Cartesian coordinates
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    return  times, x,  y,  z


#xyz til sfærisk coordinator
def solar_position_to_spherical(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)
    
    return r, theta, phi


# Theta til alfa koordinater
def solar_elevation_angle(theta):
    return np.pi/2 - theta

# solar panel elevation projection funktioner:
def solar_panel_projection(theta_s, phi_s, theta_p, phi_p):
    proj = np.sin(theta_p) * np.sin(theta_s) * np.cos(phi_p - phi_s) + np.cos(theta_p + theta_s)
    return max(proj,0)

def solar_panel_projection_arrays(theta_s, phi_s, theta_p, phi_p):
    if (len(theta_s) == len(phi_s)):
        return np.array([solar_panel_projection(theta_s[i], phi_s[i], theta_p, phi_p) for i in range(len(theta_s))])
    else: 
        raise ValueError("Arrays do not have equal number of entries")

def solar_flux(Længde, bredde, S_0, A_0, theta_panel, phi_panel, start_dato, latitude, longtitude, altitude, tz):
    """Returns time,x,y,z as np arrays"""

    delta_tid = "H" #"Min"
    tidszone = "Europe/Copenhagen"
    
    site = Location(latitude, longtitude, tz, altitude)

    # Definition of a time range of simulation
    times = pd.date_range (start_dato + " 00:00:00", "2024-06-03" + " 23:59:00", inclusive="left", freq=delta_tid, tz=tidszone)
    
    # Estimate Solar Position with the 'Location' object
    solpos = site.get_solarposition(times)
    
    # Convert zenith and azimuth from degrees to radians
    theta = np.deg2rad(np.array(solpos.loc[start_dato].zenith))
    phi = np.deg2rad(np.array(solpos.loc[start_dato].azimuth))
    print(np.deg2rad(phi_panel))
 
    proj = solar_panel_projection_arrays(theta, phi, np.deg2rad(theta_panel), np.deg2rad(phi_panel))
    
    flux = Længde * bredde * S_0 * A_0 * proj
    
    return flux
