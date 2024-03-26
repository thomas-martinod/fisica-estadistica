import numpy as np

# Constantes
pi = np.pi                  # Pi
h = 6.62607015e-34          # Plank's constant (J s)
k_B = 1.380649e-23          # Boltzmann constant (J / K)
N_A = 6.02214076e+23        #s Avogadro's Number (1 / mol)

class Tras:
    def __init__(self, m, N, T, P) -> None:
        self.m = m / 1000 / N_A  # En kg
        self.N = N
        self.T = T
        self.P = P
        self.V = self.N * k_B * self.T / self.P
        self.q = (2 * pi * self.m * k_B * self.T / h**2)**(3/2) * self.V
        self.U = 1.5 * self.N * k_B * self.T
        self.S = self.N * k_B * (2.5 + np.log(self.q / self.N))
        self.A = -self.N * k_B * self.T * (np.log(self.q / self.N) + 1)
        self.G = -self.N * k_B * self.T * np.log(self.q / self.N)
        self.H = 2.5 * self.N * k_B * self.T
        self.mu = - k_B * self.T * np.log(self.q / self.N)
    def get_q(self):
        return self.q

    def get_U(self):
        return self.U

    def get_S(self):
        return self.S

    def get_A(self):
        return self.A

    def get_G(self):
        return self.G

    def get_H(self):
        return self.H

    def get_mu(self):
        return self.mu


class Tras_and_ele:
    def __init__(self, m, N, T, P, ge1) -> None:
        self.m = m
        self.N = N
        self.T = T
        self.P = P
        self.V = self.N * k_B * self.T / self.P
        self.q_tras = (2 * pi * self.m * k_B * self.T / h**2)**(3/2) * self.V
        self.q_ele = ge1
        self.U = 1.5 * self.N * k_B * self.T
        self.S = self.N * k_B * (2.5 + np.log(self.q_tras * self.q_ele / self.N))
        self.A = -self.N * k_B * self.T * (np.log(self.q_tras * self.q_ele / self.N) + 1)
        self.G = -self.N * k_B * self.T * np.log(self.q_tras * self.q_ele / self.N)
        self.H = 2.5 * self.N * k_B * self.T
        self.mu = - k_B * self.T * np.log(self.q_tras * self.q_ele / self.N)

    def set_m(self, m):
        self.m = m

    def set_N(self, N):
        self.N = N

    def set_T(self, T):
        self.T = T

    def set_P(self, P):
        self.P = P

    def get_q(self):
        return self.q

    def get_U(self):
        return self.U

    def get_S(self):
        return self.S

    def get_A(self):
        return self.A

    def get_G(self):
        return self.G

    def get_H(self):
        return self.H

    def get_mu(self):
        return self.mu

def main():
    # Ingreso de data
    aux = input('Ingrese las condiciones de presión, temperatura y número de moles separadas por comas (si ingresa \'st\' se asume 1 bar, 298K y un mol):')
    if aux == 'st':
        P = 10000
        T = 298.15
        N = N_A
    else:
        pass
    elem = input('Ingrese el elemento a considerar: ')
    mass = float(input('Ingrese su masa: '))
    config = input('Ingrese su terminación de configuración electrónica (ej: \'3p^4\'): ')
    # Create an instance of the tras class
    sys_tras = Tras(mass, N, T, P)

    # Print out the values of q, U, S, A, G, H, and mu
    print("q =", sys_tras.get_q())
    print("U =", sys_tras.get_U())
    print("S =", sys_tras.get_S())
    print("A =", sys_tras.get_A())
    print("G =", sys_tras.get_G())
    print("H =", sys_tras.get_H())
    print("μ =", sys_tras.get_mu())

if __name__ == "__main__":
    main()

