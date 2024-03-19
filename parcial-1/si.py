import numpy as np

# Constantes
pi = np.pi                  # Pi
h = 6.62607015e-34          # Plank's constant (J s)
k_B = 1.380649e-23          # Boltzmann constant (J / K)
N_A = 6.02214076e+23        # Avogadro's Number (1 / mol)
R = k_B * N_A               # Ideal gas constant (J / mol K)


def mass_of_an_atom(m: float):
    return m / N_A / 1000  # (kg / atom)


def translational_partition(m, T, V):
    return ((2 * pi * m * k_B * T) / (h ** 2)) ** (3 / 2) * V


def internal_energy(T):
    return 3 / 2 * R * T


def entropy(q, N):
    return N * (3 / 2 * k_B + k_B * np.log(q))


def gibbs(q, N, T):
    return N * k_B * T * (1 - np.log(q))


def electronic_partition(ge):
    return ge


def entropy_ele(q_ele):
    return R * np.log(q_ele)


def entropy_both(q_tras, q_ele):
    return 5 / 2 * R + R * np.log(q_tras * q_ele / N_A)


def main():
    T = 298  # Standard Temperature (K)
    P = 10 ** 5  # Standard pressure (Pa)
    V = k_B * T / P  # Standard volume (m^3)
    N = N_A  # Number of particles

    print('Con este programa calculas U_298, S_298, G_298.')
    atom_in = input('Ingrese el átomo: ')
    m_in = float(input('Ingrese la masa: '))  # Convert to float
    ge_in = float(input('Ingrese la degenerancia del estado electrónico: '))  # Convert to float

    print('--------------------------------------------------------')
    print('Función de partición traslacional')
    m = mass_of_an_atom(m_in)
    q_tras = translational_partition(m, T, V)
    print('q_tras = ' + str(q_tras))
    print('\nU°_298 = ' + str(internal_energy(T) / 1000) + 'kJ')
    print('S°_298 = ' + str(entropy(q_tras, N)) + 'unidad')
    print('G°_298 = ' + str(gibbs(q_tras, N, T)))


main()
