#indeholder alle de funktioner, som der står vi skal skrive.

#funktionsværdier i `f` i intervallet $[Low, Up]$ og angiver de tilhørende `t`-værdier.
def interval(f, t, low, up):
    for i in range(0, len(f)):
        if f[i] > low and f[i] < up:
            print("Data(f): ", f[i], "Time(t): ", t[i])
    return None
        

# Fortegn skift
def fortegnskift(f, t):
    for i in range(0, len(f) - 1):
        if f[i] * f[i+1] < 0:
            print("Skift ved: ", t[i])
    if f[len(f) - 1] * f[len(f)- 2] < 0:
        print("Skift ved: ", t[i])
    return None
