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

# Theta til alfa koordinater
def solar_elevation_angle(theta):
    return np.pi/2 - theta

# solar panel elevation projection funktioner:
def solar_panel_projection(theta_s, phi_s, theta_p, phi_p):
    proj = sin(theta_p) * sin(theta_s) * cos(phi_p - phi_s) + cos(theta_p + theta_s)
    return max(proj,0)

def solar_panel_projection_arrays(theta_s, phi_s, theta_p, phi_p):
    if (len(theta_s) == len(phi_s) == len(theta_p) == len(phi_p)):
        return np.array([solar_panel_projection(theta_s[i], phi_s[i], theta_p[i], phi_p[i]) for i in range(len(theta_s))])
    else: 
        raise ValueError("Arrays do not have equal number of entries")
