import numpy as np

# Constants
pi = np.pi
R = 8.314462618  # J/(mol*K) - Ideal gas constant
T = 298.0        # K - Temperature
P = 1.0          # atm - Pressure
h = 6.62607015e-34
kB = 1.380649e-23

# Translational partition function
def translational_partition(m):
    return ((2.0 * pi * m * kB * T) / (h**2))**(3.0 / 2.0) * (kB * T / P)

# Mass (in kg/mol)
mass = 28.97e-3  # Molar mass of nitrogen in kg/mol

# Calculate the translational partition function for the given mass
translational_partition_value = translational_partition(mass)

# Print the result
print("Molar mass:", mass, "kg/mol")
print("Translational partition function:", translational_partition_value)
