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
print('Standard internal energy: ' + str(U) + '\n')


# He, Ne, Al, Ba, Ar, Be, Br, C
elements = ['He', 'Ne', 'Al', 'Ba', 'Ar', 'Be', 'Br', 'C']

# List of molar masses of the elements listed (g / mol)
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
    return N_A*k_B*T*np.log(np.e/q)

for i in range(len(elements)):
    print('Element:\t' + elements[i])
    mass = mass_of_an_atom(molar_masses[i])
    print('Mass of an atom =\t' + str(mass) + ' kg')
    Q = translational_partition(mass)
    print('q_tras = \t' + str(Q))
    print('S_298 = \t' + str(entropy(Q)))
    print('G_298 = \t' + str(gibbs(Q)) + '\n')

