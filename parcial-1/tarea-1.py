import numpy as np

# Constants
pi = np.pi
R = 8.314462618
T = 298
P = 1.0
h = 6.62607015e-34
kB = 1.380649e-23
NA = R/kB
V = 4.11e-26
# He, Ne, Al, Ba, Ar, Be, Br, C
molar_masses = [4.00260128, 20.18004638, 26.98153841, 137.32667172, 39.94779856, 9.01218306, 79.90432616, 12.01063556]
Qe=[1.00, 1.00, 2.00, 1.00, 1.00, 1.00, 2.00, 1.00] # Vector para estado de ge1 de cada elemento de la tabla

def get_mass_of_an_atom (molarmass: float):
    return molarmass / NA / 1000  # Se retorna en kg/atomo

# Translational partition function
def translational_partition(m):
    return ((2 * pi * m * kB * T) / (h**2))**(3/2) * V

def entropia(q):
    return NA*(3/2 * kB + kB * np.log(q))

def energia_interna():
    return 3/2 * R * T

def gibbs(q):
    return kB*T*np.log(np.e/q)

mass_he = get_mass_of_an_atom(molar_masses[7])
print(mass_he)
tras_He = translational_partition(mass_he)
print(tras_He)
print(entropia(tras_He))
print(energia_interna())
print(gibbs(tras_He))



