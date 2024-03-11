import numpy as np

# Constants
pi = np.pi                  # Pi
T = 298                     # Standard Temperature (K)
P = 10**5                   # Standard pressure (Pa)
h = 6.62607015e-34          # Plank's constant (J s)
k_B = 1.380649e-23          # Boltzmann constant (J / K)
N_A = 6.02214076e-23        # Avogadro's Number (1 / mol)

R = k_B * N_A               # Ideal gas constant (J / mol K)
V = k_B * T / P             # Standard volume (m^3)
U = 3/2 * R * T /1000       # Standard internal energy U_298 (kJ)

print('Standard Volume: ' + str(V))


# List of molar masses of the elements listed (g / mol)
# He, Ne, Al, Ba, Ar, Be, Br, C
molar_masses = [4.00260128, 20.18004638, 26.98153841, 137.32667172, 39.94779856, 9.01218306, 79.90432616, 12.01063556]

# List of degeneracies (ge1)
degeneracies = [1, 1, 2, 1, 1, 1, 2, 1]


def mass_of_an_atom (m: float):
    return m / N_A / 1000  # (kg / atom)

def translational_partition(m):
    return ((2 * pi * m * k_B * T) / (h**2))**(3/2) * V

def entropy(q):
    return N_A*(3/2 * k_B + k_B * np.log(q))

def gibbs(q):
    return k_B*T*np.log(np.e/q)


