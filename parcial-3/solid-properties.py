import numpy as np

# Constants
hbar = 1.054571817e-34  # Planck's constant over 2 pi in J·s
kB = 1.380649e-23  # Boltzmann constant in J/K
me = 9.1093837e-31  # Mass of electron in kg
eV = 1.60218e-19  # 1 eV in joules

def calculate_volume_and_density(bravais, Z):
    a = float(input('Enter the lattice parameter a (in Ångströms): '))
    a_meters = a * 1e-10  # Convert a from Ångströms to meters

    if bravais == 'cubic':
        N = 1
        V = a_meters**3

    elif bravais == 'fcc':
        N = 4
        V = a_meters**3

    elif bravais == 'bcc':
        N = 2
        V = a_meters**3

    elif bravais == 'tetra':
        c = float(input('Enter the lattice parameter c (in Ångströms): '))
        c_meters = c * 1e-10
        N = 1
        V = a_meters**2 * c_meters

    elif bravais == 'bc-tetra':
        c = float(input('Enter the lattice parameter c (in Ångströms): '))
        c_meters = c * 1e-10
        N = 2
        V = a_meters**2 * c_meters / 2

    elif bravais == 'ortho':
        b = float(input('Enter the lattice parameter b (in Ångströms): '))
        c = float(input('Enter the lattice parameter c (in Ångströms): '))
        b_meters, c_meters = b * 1e-10, c * 10e-10
        N = 1
        V = a_meters * b_meters * c_meters

    elif bravais == 'bc-ortho':
        b = float(input('Enter the lattice parameter b (in Ångströms): '))
        c = float(input('Enter the lattice parameter c (in Ångströms): '))
        b_meters, c_meters = b * 1e-10, c * 10e-10
        N = 2
        V = a_meters * b_meters * c_meters / 2

    elif bravais == 'fc-ortho':
        b = float(input('Enter the lattice parameter b (in Ångströms): '))
        c = float(input('Enter the lattice parameter c (in Ångströms): '))
        b_meters, c_meters = b * 1e-10, c * 10e-10
        N = 4
        V = a_meters * b_meters * c_meters / 4

    elif bravais == 'hex':
        c = float(input('Enter the lattice parameter c (in Ångströms): '))
        c_meters = c * 1e-10
        N = 1
        V = (3 * np.sqrt(3) / 2) * a_meters**2 * c_meters

    elif bravais == 'rhombo':
        alpha = float(input('Enter the angle alpha (in degrees): '))
        alpha_radians = np.radians(alpha)
        N = 1
        V = a_meters**3 * np.sqrt(1 - 3 * np.cos(alpha_radians)**2 + 2 * np.cos(alpha_radians)**3)

    elif bravais == 'mono':
        b = float(input('Enter the lattice parameter b (in Ångströms): '))
        c = float(input('Enter the lattice parameter c (in Ångströms): '))
        beta = float(input('Enter the angle beta (in degrees): '))
        b_meters, c_meters = b * 1e-10, c * 10e-10
        beta_radians = np.radians(beta)
        N = 1
        V = a_meters * b_meters * c_meters * np.sin(beta_radians)

    elif bravais == 'bc-mono':
        b = float(input('Enter the lattice parameter b (in Ångströms): '))
        c = float(input('Enter the lattice parameter c (in Ångströms): '))
        beta = float(input('Enter the angle beta (in degrees): '))
        b_meters, c_meters = b * 1e-10, c * 10e-10
        beta_radians = np.radians(beta)
        N = 2
        V = a_meters * b_meters * c_meters * np.sin(beta_radians) / 2

    elif bravais == 'tri':
        b = float(input('Enter the lattice parameter b (in Ångströms): '))
        c = float(input('Enter the lattice parameter c (in Ångströms): '))
        alpha = float(input('Enter the angle alpha (in degrees): '))
        beta = float(input('Enter the angle beta (in degrees): '))
        gamma = float(input('Enter the angle gamma (in degrees): '))
        b_meters, c_meters = b * 1e-10, c * 10e-10
        alpha_radians = np.radians(alpha)
        beta_radians = np.radians(beta)
        gamma_radians = np.radians(gamma)
        N = 1
        V = a_meters * b_meters * c_meters * np.sqrt(1 - np.cos(alpha_radians)**2 - np.cos(beta_radians)**2 - np.cos(gamma_radians)**2 + 2 * np.cos(alpha_radians) * np.cos(beta_radians) * np.cos(gamma_radians))

    else:
        raise ValueError("Unsupported Bravais lattice type")

    rho = N * Z / V  # Electron density in electrons/m^3
    return V, rho


def calculate_fermi_energy(rho):
    return hbar**2 / (2 * me) * (3 * np.pi**2 * rho)**(2/3)


def calculate_fermi_velocity(rho):
    return hbar / me * (3 * np.pi**2 * rho)**(1/3)


def calculate_chemical_potential(T, Ef):
    return Ef * (1 - np.pi**2 / 12 * (kB * T / Ef)**2)


def calculate_energy(T, Ef):
    return Ef * (1 + 5 * np.pi**2 / 12 * (kB * T / Ef)**2)


def main():
    # Prompt user for inputs
    bravais_cell = input("\nEnter the type of Bravais lattice cell: ")
    valence_electrons = float(input("Enter the number of valence electrons Z: "))

    # Call the function with user inputs
    V, rho = calculate_volume_and_density(bravais_cell, valence_electrons)
    T = float(input('Enter the temperature in K: '))
    print('')

    Ef_J = calculate_fermi_energy(rho)  # Fermi energy in joules
    Vf = calculate_fermi_velocity(rho)  # Fermi velocity in m/s
    mu_J = calculate_chemical_potential(T, Ef_J)  # Chemical potential in joules
    E_J = calculate_energy(T, Ef_J)  # Energy in joules

    # Convert energies from joules to eV
    Ef_eV = Ef_J / eV
    mu_eV = mu_J / eV
    E_eV = E_J / eV

    # Print results in two columns with formatted alignment
    print('-------------------------------------------------------------------------\n')
    print(f"{'Volume:':<40}{V:.4e} m^3")
    print(f"{'Electron Density:':<40}{rho:.4e} electrons/m^3")
    print(f"{'Fermi Energy:':<40}{Ef_eV:.4e} eV")
    print(f"{'Fermi Velocity:':<40}{Vf:.4e} m/s")
    print(f"{'Chemical Potential at T = ' + str(T) + ' K:':<40}{mu_eV:.4e} eV")
    print(f"{'Energy at T = ' + str(T) + ' K:':<40}{E_eV:.4e} eV")
    print(f"{'Chemical Potential (eV/mol):':<40}{mu_eV:.4e} eV/mol\n")
    print('-------------------------------------------------------------------------\n')

main()