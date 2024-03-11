import numpy as np

# Constants
pi = np.pi
R = 8.314462618
T = 300
P = 1.0
h = 6.62607015e-34
kB = 1.380649e-23
NA = R/kB
# He, Ne, Al, Ba, Ar, Be, Br, C
molar_masses = [4.00260128, 20.18004638, 26.98153841, 137.32667172, 39.94779856, 9.01218306, 79.90432616, 12.01063556]

def get_mass_of_an_atom (molarmass: float):
    return molarmass / NA / 1000  # Se retorna en kg/atomo

# Translational partition function
def translational_partition(m):
    return ((2 * pi * m * kB * T) / (h**2))**(3/2) * (kB * T / P)

def entropia(q):
    return 3/2 * kB * np.log(q)

mass_he = get_mass_of_an_atom(molar_masses[0])
print(mass_he)
tras_He = translational_partition(mass_he)
print(tras_He)
print(entropia(tras_He))

